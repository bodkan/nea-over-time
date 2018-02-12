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

gen_time <- 25
Ne0 <- 10000

plot_nea <- function(df, real_lm=NULL, sim_lm=NULL) {
  p <- ggplot(df, aes(t_admix, nea)) + geom_point(size=5) + xlim(0, 55000) + ylim(0, 0.05)
  if (!is.null(real_lm)) {
    p <- p + geom_abline(intercept=coef(real_lm)[1], slope=coef(real_lm)[2], linetype=2, color="blue")
  }
  if (!is.null(sim_lm)) {
    p <- p + geom_abline(intercept=coef(sim_lm)[1], slope=coef(sim_lm)[2], linetype=2, color="red")
  }
  p
}

# read admixture array estimates and ages
real_data <- function(path) {
  read_tsv(path, col_types="cidci") %>%
  filter(name != "Oase1") %>%
  group_by(pop) %>% 
  mutate(name=ifelse(pop == "EMH", paste0("emh_", 1:n()), paste0("eur_", 1:n())),
         t_admix=55000 - age) %>% 
  ungroup
}

# calculate Neanderthal and African allele frequencies
pop_freq <- function(dt, pop) {
  dt[, rowSums(.SD) / length(.SD), .SDcols=colnames(dt)[str_detect(colnames(dt), pop)]]
}

admixture_array <- function(all_sites) {
  afr_freq <- pop_freq(all_sites, "afr")
  nea_freq <- pop_freq(all_sites, "nea")
  
  # process the sites based on the archaic admixture array conditioning
  afr_cutoff <- 1.0
  array_sites <- all_sites[(afr_freq == 0 | afr_freq >= afr_cutoff) &
                           (nea_freq == 0 | nea_freq == 1) &
                           (abs(afr_freq - nea_freq) > 0.5)]
  
  array_sites
}

big_yoruba_array <- function(all_sites) {
  all_sites[(nea_1 != nea_2) | (afr_1 != afr_2) | (afr_3 != afr_4)]
}

real_eur <- real_data("data/admixture_array_nea.tsv")

scale_t <- function(t) { t / gen_time / (4 * Ne0) }
scale_Ne <- function(Ne) { Ne / Ne0 }
scale_m <- function(m) { 4 * m * Ne0 }


simulate_sites <- function(p0, m, t, d,
                           n_afr=0, n_nea=0, n_eur=0, n_asn=0, n_chimp=0,
                           n_haps=20000, hap_length=5001,
                           emh_ages=NULL)  {
  if (n_afr + n_nea + n_eur + n_asn + n_chimp == 0)
    stop("No haplotypes sampled")
  
  # mutation rate per site per generation
  mut_rate <- 2.5e-8
  # probability of cross-over between adjacent bases per generation 
  recomb_rate <- 1e-8
  
  gen_time <- 25

  # effective population size used for scaling below
  Ne0 <- 20000

  # how old are the sampled Neanderthal haplotypes?
  T_nea_age <- 70000
  
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
  nea_admix_duration <- 1000
  T_nea_admix_end <- T_nea_admix_start - nea_admix_duration
  m_nea_admix <- p0 / (nea_admix_duration / gen_time)
  
  # population scaled mutation and recombination rates
  theta <- 4 * Ne0 * mut_rate * hap_length
  rho <- 4 * Ne0 * recomb_rate * hap_length
  
  T_nea_gens <- scale_t(T_nea_age)
  emh_gens <- map_dbl(emh_ages, scale_t)

  eI <- c(
    glue("-eI {T_nea_gens} 0       0       0       {n_nea} 0"),
    glue("-eI {emh_gens}   0       1       0       0       0"),
    glue("-eI 0            {n_afr} 0       0       0       0"),
    glue("-eI 0            0       {n_eur} 0       0       0"),
    glue("-eI 0            0       0       {n_asn} 0       0"),
    glue("-eI 0            0       0       0       0       {n_chimp}")
  ) %>%
    paste(collapse=" ")
  
  scrm_cmd <- glue("
  
  scrm \\
  
  {n_nea + n_afr + n_asn + n_eur + length(emh_ages) + n_chimp} {n_haps} \\
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
                 paste0("emh_", 1 : length(emh_ages)),
                 paste0("afr_", 1 : n_afr),
                 paste0("eur_", 1 : n_eur))
  
  # run scrm and parse all simulated segregating sites
  all_sites <- fread(scrm_cmd, col.names=col_names, showProgress=FALSE)
  
  all_sites
}


# detect all Neanderthal alleles in the non-African population (TRUE/FALSE at
# each site converted to integers - 1 is Nea-like allele, 0 is an African-like
# MH allele)
sim_samples <- map_df(admix_array[, real_eur$name, with=FALSE],
                      ~ mean(. == admix_array$nea_1)) %>%
  gather(name, nea) %>% inner_join(select(real_eur, -nea), ., by="name")

sim_samples
