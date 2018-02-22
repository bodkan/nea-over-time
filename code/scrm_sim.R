suppressPackageStartupMessages({
library(tidyverse)
library(glue)
library(data.table)
library(purrr)
})

gen_time <- 25

scale_t_fun <- function(t, Ne0) { t / gen_time / (4 * Ne0) }
scale_Ne_fun <- function(Ne, Ne0) { Ne / Ne0 }
scale_m_fun <- function(m, Ne0) { 4 * m * Ne0 }

make_names <- function(prefix, num) {
  if (!num) return(NULL)
  paste0(prefix, "_", 1 : num)
}

simulate_sites <- function(p0, m_afr_eur, m_eur_afr, t, d, Ne0,
                           n_afr=0, n_nea=0, n_eur=0, n_asn=0, n_chimp=0,
                           n_haps=20000, hap_length=5001,
                           emh_ages=NULL)  {
  if (n_afr + n_nea + n_eur + n_asn + n_chimp == 0)
    stop("No haplotypes sampled")
  
  scale_t <- partial(scale_t_fun, Ne0 = Ne0)
  scale_Ne <- partial(scale_Ne_fun, Ne0 = Ne0)
  scale_m <- partial(scale_m_fun, Ne0 = Ne0)

  # mutation rate per site per generation
  mut_rate <- 2.5e-8
  # probability of cross-over between adjacent bases per generation 
  recomb_rate <- 1e-8
  
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
  
  -em {scale_t(T_m_afr_eur_start)} 1 2 {scale_m(m_eur_afr)} \\
  -em {scale_t(T_m_afr_eur_end)} 1 2 0 \\

  -em {scale_t(T_m_afr_eur_start)} 2 1 {scale_m(m_afr_eur)} \\
  -em {scale_t(T_m_afr_eur_end)} 2 1 0 \\
  
  -ej {scale_t(30000)} 3 2 \\
  -ej {scale_t(T_afr_nonafr_split)} 2 1 \\
  -ej {scale_t(T_nea_mh_split)} 4 1 \\
  -ej {scale_t(T_chimp_split)} 5 1 \\
  
  {eI} \\
  
  --transpose-segsites | tail -n+4 | grep '^[0-9]' | cut -d' ' -f3-
  
  ") %>% str_replace_all("\n", " ") %>% str_replace_all("  +", " ")
  
  # generate column names for the output df
  col_names <- c(make_names("nea",   n_nea),
                 make_names("emh",   length(emh_ages)),
                 make_names("afr",   n_afr),
                 make_names("eur",   n_eur),
                 make_names("asn",   n_asn),
                 make_names("chimp", n_chimp))
  
  # run scrm and parse all simulated segregating sites
  all_sites <- fread(scrm_cmd, col.names=col_names, showProgress=FALSE)
  
  all_sites
}

