suppressPackageStartupMessages({
library(tidyverse)
library(glue)
library(data.table)
})

real_path <- "data/admixture_array_nea.tsv"
p0 <- 0.03
m <- 1e-4
t <- 15000
d <- 15000
output_dir <- NULL

real_samples <- read_tsv(real_path, col_types="cidci") %>%
  filter(name != "Oase1") %>%
  group_by(pop) %>% 
  mutate(name=ifelse(pop == "EMH", paste0("emh_", 1:n()), paste0("eur_", 1:n())),
         t_admix=55000 - age) %>% 
  ungroup

scale_t <- function(t) { t / gen_time / (4 * Ne0) }
scale_Ne <- function(Ne) { Ne / Ne0 }
scale_m <- function(m) { 4 * m * Ne0 }

# generation time
gen_time <- 25

# effective population size used for scaling below
Ne0 <- 20000

# how old are the sampled Neanderthal haplotypes?
T_nea_age <- 70000

# how many haplotypes to sample from each population?
n_afr <- 200
n_asn <- 0 # 50
n_nea <- 4 # 6
n_chimp <- 0 # 1

n_haplotypes <- 20000
hap_length <- 5001

# split between Homo and Chimp
T_chimp_split <- 6000000
# split between MH and Neanderthals
T_nea_mh_split <- 600000
# African/non-African split and bottleneck
T_afr_nonafr_split <- 60000

# Ne during the Out of Africa bottleneck
Ne_bottleneck <- 2000

# size of the Neanderthal population
Ne_nea <- 1000

# size of the ancestral population
Ne_anc <- 10000

# EUR <-> AFR migration rates
m_afr_eur <- m
T_m_afr_eur_start <- t - d
T_m_afr_eur_end <- t

# Neanderthal introgression parameters
T_nea_admix_start <- 55000
T_nea_admix_end <- T_nea_admix_start - 1000
m_nea_admix <- p0 / (1000 / gen_time)

# mutation rate per site per generation
mut_rate <- 2.5e-8
# probability of cross-over between adjacent bases per generation 
recomb_rate <- 1e-8

# population scaled mutation and recombination rates
theta <- 4 * Ne0 * mut_rate * hap_length
rho <- 4 * Ne0 * recomb_rate * hap_length



n_eur <- filter(real_samples, pop != "EMH") %>% nrow
emh_ages <- filter(real_samples, pop == "EMH")$age


eI <- c(
  glue("-eI {scale_t(T_nea_age)}         0       0       0       {n_nea}   0"),
  glue("-eI {map_dbl(emh_ages, scale_t)} 0       1       0       0         0"),
  glue("-eI 0                            {n_afr} 0       0       0         0"),
  glue("-eI 0                            0       {n_eur} 0       0         0"),
  glue("-eI 0                            0       0       {n_asn} 0         0"),
  glue("-eI 0                            0       0       0       0 {n_chimp}")
) %>%
  paste(collapse=" ")



scrm_cmd <- glue("

scrm \\

{n_nea + n_afr + n_asn + nrow(real_samples) + n_chimp} {n_haplotypes} \\
-I 5 0 0 0 0 0 \\
-r {rho} {hap_length} \\
-t {theta} \\
-n 4 {scale_Ne(Ne_nea)} \\

-em {scale_t(T_nea_admix_end)} 2 4 {scale_m(m_nea_admix)} \\
-em {scale_t(T_nea_admix_start)} 2 4 0 \\
-en {scale_t(T_nea_admix_start)} 2 {scale_Ne(Ne_bottleneck)} \\

-em {scale_t(T_m_afr_eur_start)} 1 2 {scale_m(m_afr_eur)} \\
-em {scale_t(T_m_afr_eur_end)} 1 2 0 \\

-ej {scale_t(30000)} 3 2 \\
-ej {scale_t(T_afr_nonafr_split)} 2 1 \\
-ej {scale_t(T_nea_mh_split)} 4 1 \\
-ej {scale_t(T_chimp_split)} 5 1 \\

{eI} \\

--transpose-segsites | tail -n+4 | grep '^[0-9]' | cut -d' ' -f3-

") %>% str_replace_all("\n", " ") %>% str_replace_all("  +", " ")

# generate column names for the output df
col_names <- c(paste0("nea_", 1 : n_nea),
               filter(real_samples, pop == "EMH")$name,
               paste0("afr_", 1 : n_afr),
               filter(real_samples, pop != "EMH")$name)

# run scrm and parse all simulated segregating sites
all_sites <- fread(scrm_cmd, col.names=col_names, showProgress=FALSE)

# calculate Neanderthal and African allele frequencies
pop_freq <- function(dt, pop) {
  dt[, rowSums(.SD) / length(.SD), .SDcols=colnames(dt)[str_detect(colnames(dt), pop)]]
}
afr_freq <- pop_freq(all_sites, "afr")
nea_freq <- pop_freq(all_sites, "nea")

# process the sites based on the archaic admixture array conditioning
afr_cutoff <- 1.0
admix_array <- all_sites[(afr_freq == 0 | afr_freq >= afr_cutoff) &
                         (nea_freq == 0 | nea_freq == 1) &
                         (abs(afr_freq - nea_freq) > 0.5)]
rm(all_sites); gc()

# detect all Neanderthal alleles in the non-African population (TRUE/FALSE at
# each site converted to integers - 1 is Nea-like allele, 0 is an African-like
# MH allele)
sim_samples <- map_df(admix_array[, real_samples$name, with=FALSE],
                  ~ mean(. == admix_array$nea_1)) %>%
  gather(name, nea) %>% inner_join(select(real_samples, -nea), ., by="name")




nea_est %>%
  ggplot(aes(age, nea)) +
  geom_point(size=5) + xlim(55000, 0) + ylim(0, 0.05) +
  geom_smooth(method="lm", linetype=2, color="red", alpha=1/2, fullrange=TRUE) +
  geom_hline(yintercept=m_nea_admix * (T_nea_admix_start - T_nea_admix_end) / gen_time, linetype=2, alpha=1/2) +
  xlab("years after admixture") + ylab("Nea ancestry proportion") +
  geom_vline(xintercept=c(T_m_afr_eur_start, T_m_afr_eur_end), linetype=2, color="blue", alpha=1/2) +
  theme(text=element_text(size=15))
