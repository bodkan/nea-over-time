library(tidyverse)
library(admixstats)

# Plot Neanderthal ancestry over time.
plot_nea <- function(df, est_lm = FALSE, real_lm = NULL, sim_lm = NULL) {
  p <- ggplot(df, aes(t_admix, nea)) + geom_point(alpha = 1/2, size = 5) + xlim(0, 55000) + ylim(0, 0.08)
  if (!is.null(real_lm)) {
    p <- p + geom_abline(intercept=coef(real_lm)[1], slope=coef(real_lm)[2], linetype=2, color="blue")
  }
  if (!is.null(sim_lm)) {
    p <- p + geom_abline(intercept=coef(sim_lm)[1], slope=coef(sim_lm)[2], linetype=2, color="red")
  }
  if (est_lm) {
    p <- p + geom_smooth(method = "lm", linetype = 2, color = "red")
  }
  p
}

# Read empirical Neanderthal ancestry estimates.
load_estimates <- function(path) {
  read_tsv(path, col_types="cidci") %>%
    filter(name != "Oase1") %>%
    group_by(pop) %>% 
    mutate(name=ifelse(pop == "EMH", paste0("emh_", 1:n()), paste0("eur_", 1:n())),
           t_admix=55000 - age) %>% 
    ungroup
}

nea_lm <- function(df) {
  lm(nea ~ t_admix, data=df)
}

direct_props <- function(sites, sample_names) {
  direct_prop(sites, sample_names) %>% inner_join(select(real_est, -nea), by = "name")
}

f4_ratios <- function(sites, sample_names) {
  f4_ratio(sites,
           x = sample_names,
           a = c("nea_1", "nea_2"), b = c("nea_3", "nea_4"), c = paste0("afr_", 1:10),
           o = "chimp_1") %>% inner_join(select(real_est, -nea), by = "name")
}

real_est <- load_estimates("data/admixture_array_nea.tsv")

library(parallel)
sims_df <- data.frame()
for (Ne in c(10000, 20000, 30000, 40000, 50000)) {
  for (m in c(0, 0.0001, 0.001)) {
    for (t in c(20000, 15000, 10000, 5000)) {

reps <- mclapply(1:20, mc.cores = 10, function(rep_i) {

sites <- nea_admix_scrm(0.03, 0, m, t, t, Ne_eur = Ne,
                        n_afr = 200,
                        n_eur = nrow(filter(real_est, pop != "EMH")),
                        n_asn = 2, n_nea = 4, n_chimp = 1,
                        n_haps = 100000, hap_length = 5001,
                        emh_ages = filter(real_est, pop == "EMH")$age)
 
admix_array <- archaic_admixture_array(sites)
bigyri_array <- big_yoruba_array(sites)
ho_array <- human_origins_array(sites)

sample_names <- filter(real_est, !name %in% c("eur_1", "eur_2", "asn_1", "asn_2"))$name

admix <- direct_props(admix_array, sample_names) %>% mutate(sites = "admix")
bigyri <- f4_ratios(bigyri_array, sample_names) %>% mutate(sites = "bigyri")
ho <- f4_ratios(ho_array, sample_names) %>% mutate(sites = "ho")
all <- f4_ratios(sites, sample_names) %>% mutate(sites = "all")

bind_rows(admix, bigyri, ho, all) %>% mutate(Ne = Ne, m = m, t = t, rep = rep_i)

}) %>% bind_rows

sims_df <- bind_rows(sims_df, reps)

    }
  }
}








# plot_nea(admix, est_lm = T)
# plot_nea(bigyri, est_lm = T)
# plot_nea(ho, est_lm = T)
# plot_nea(all, est_lm = T)



# nea_admix_scrm(0.03, 1e-4, 1e-4, 15000, 15000, Ne_eur = 20000,
# n_afr = 10,
# n_eur = nrow(filter(real_est, pop != "EMH")),
# n_asn = 6, n_nea = 4, n_chimp = 1,
# n_haps = 100000, hap_length = 5001,
# emh_ages = filter(real_est, pop == "EMH")$age, dump_cmd = TRUE) %>% PlotMS(type = "scrm", time.scale = "log10year")











# # Calculate the ratio between magnitudes of slopes in real vs simulated data.
# slope_ratio <- function(real_lm, sim_lm) {
#   as.vector(coef(real_lm)[2] / coef(sim_lm)[2])
# }
# 
# load_abc <- function(path, prefix) {
#   list.files(path, paste0(prefix, ".*"), full.names=TRUE) %>%
#     map_df(readRDS)
# }



# real_eur <- load_estimates("data/admixture_array_nea.tsv")
# 
# real_lm <- lm(nea ~ t_admix, data=real_eur, weights=snp_count)
# 
# MIN_M <- 0.00001; MAX_M <- 0.001
# MIN_T <- 1; MAX_T <- 18000
# 
# prior_p0 <- function() { as.vector(coef(real_lm)[1]) }
# prior_m <- function() { runif(1, MIN_M, MAX_M) }
# prior_t <- function() { sample(MIN_T : MAX_T, 1) }
# prior_d <- function(t) { sample(1 : t, 1) }
# 
# 
# all_sites <- simulate_sites(
#   p0 = prior_p0(), m = 1e-4, t = 15000, d = 15000,
#   n_afr = 200, n_eur = nrow(filter(real_eur, pop != "EMH")),
#   n_nea = 4, n_chimp = 1, n_haps = 10, hap_length = 5001,
#   emh_ages = filter(real_eur, pop == "EMH")$age
# )
# 
# admix_array <- archaic_admixture_array(all_sites)
# bigyri_array <- big_yoruba_array(all_sites)
# 
# 
# real_lm <- lm(nea ~ t_admix, data=real_eur, weights=snp_count)
# sim_lm <- lm(nea ~ t_admix, data=sim_eur)





















# 
# run_iter <- function() {
#   p0 <- prior_p0()
#   m <- prior_m()
#   t <- prior_t()
#   d <- prior_d(t)
#   
#   sim <- simulate_data(p0, m, samples$gen)
#   
#   sim_lm <- lm(sim ~ samples$t_admix)
#   ratio <- slope_ratio(real_lm, sim_lm)
# 
#   res <- data.frame(p0, m, t, d, ratio)
#   
#   if (ratio > 0.98 & ratio < 1.02) {
#     #plot_nea(samples, real_lm, sim_lm
#     saveRDS(res, tempfile(tmpdir="tmp/abc", pattern="abc_", fileext=".rds"))
#   }
# }
# 



# unlink("tmp/abc/*")
# mclapply(1:5000, function(i) run_sim(), mc.cores=detectCores()) -> .
# 
# 
# x <- load_abc("tmp/abc", "abc_")

# filter(x, ratio > 0.99, ratio < 1.01) %>% sample_n(1) %>%
# {
#   sim <- simulate_data(.[["p0"]], .[["m"]], samples$gen)
#   sim_lm <- lm(sim ~ samples$t_admix)
#   print(.)
#   plot_nea(samples, real_lm, sim_lm)
# }
