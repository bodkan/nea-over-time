library(tidyverse)
library(parallel)

# load("data/RData/admixtools_nea_ancestry.RData")
# filter(nea_estimates, method == "direct") %>%
#   select(name, age, alpha, pop, snp_count) %>%
#   rename(nea=alpha) %>%
#   write_tsv("data/admixture_array_nea.tsv")

simulate_decrease <- function(p0, m, total_time) {
  p_trajectory <- c(p0, rep(0, total_time - 1))
  
  for (gen in seq_len(total_time)[-1]) {
    p_trajectory[gen] <- p_trajectory[gen - 1] * (1 - m)
  }
  
  p_trajectory
}

sample_datapoints <- function(sample_times, traj, snp_count) {
  sapply(sample_times, function(gen) { mean(rbinom(snp_count, 1, traj[gen]))  })
}

simulate_data <- function(p0, m, sample_times) {
  total_time <- as.integer(55000 / 25) # generations
  
  # simulate the dilution of Nea. alleles with
  traj <- simulate_decrease(p0, m, total_time)
  
  # sample aDNA samples
  sample_datapoints(sample_times, traj, snp_count=100000)
}

plot_nea <- function(df, real_lm=NULL, sim_lm=NULL) {
  p <- ggplot(df, aes(t_admix, nea)) + geom_point(size=5) + xlim(0, 55000) + ylim(0, 0.05)
  if (!is.null(real_lm)) {
    p <- p + geom_abline(intercept=coef(real_lm)[1], slope=coef(real_lm)[2], color="red")
  }
  if (!is.null(sim_lm)) {
    p <- p + geom_abline(intercept=coef(sim_lm)[1], slope=coef(sim_lm)[2], color="blue")
  }
  p
}

real_lm <- lm(nea ~ t_admix, data=real_samples)
sim_lm <- lm(nea ~ t_admix, data=sim_samples)

slope_ratio <- function(real_lm, sim_lm) {
  as.vector(coef(real_lm)[2] / coef(sim_lm)[2])
}

prior_p0 <- function() { runif(1, coef(real_lm)[1] - 0.0025, coef(real_lm)[1] + 0.0025) }
prior_m <- function() { runif(1, 0.00001, 0.001) }
prior_t <- function() { as.integer(runif(1, 2000, 18000)) }
prior_d <- function(t) { as.integer(runif(1, 2000, 20000 - t)) }



samples <- read_tsv("data/admixture_array_nea.tsv") %>%
  mutate(gen=as.integer(55000 / 25) - as.integer(age / 25),
         t_admix=55000 - age)
real_lm <- lm(nea ~ t_admix, data=samples, weights=snp_count)

# sim <- simulate_data(0.031, 0.000230525548254992, samples$gen)
# sim <- simulate_data(0.02041546, 0.0004839762)
# sim_lm <- lm(sim ~ samples$age)
# plot_nea(samples, real_lm, sim_lm)


run_sim <- function() {
  p0 <- prior_p0()
  m <- prior_m()
  t <- prior_t()
  d <- prior_d(t)
  
  sim <- simulate_data(p0, m, samples$gen)
  
  sim_lm <- lm(sim ~ samples$t_admix)
  ratio <- slope_ratio(real_lm, sim_lm)

  res <- data.frame(p0, m, t, d, ratio)
  
  if (ratio > 0.98 & ratio < 1.02) {
    #plot_nea(samples, real_lm, sim_lm
    saveRDS(res, tempfile(tmpdir="tmp/abc", pattern="abc_", fileext=".rds"))
  }
}

load_abc <- function(path, prefix) {
  list.files(path, paste0(prefix, ".*"), full.names=TRUE) %>%
    map_df(readRDS)
}


unlink("tmp/abc/*")
mclapply(1:5000, function(i) run_sim(), mc.cores=detectCores()) -> .


x <- load_abc("tmp/abc", "abc_")

# filter(x, ratio > 0.99, ratio < 1.01) %>% sample_n(1) %>%
# {
#   sim <- simulate_data(.[["p0"]], .[["m"]], samples$gen)
#   sim_lm <- lm(sim ~ samples$t_admix)
#   print(.)
#   plot_nea(samples, real_lm, sim_lm)
# }
