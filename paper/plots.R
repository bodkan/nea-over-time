# install.packages("devtools")
# devtools::install_github("thomasp85/patchwork")
# devtools::install_github('cttobin/ggthemr')

library(tidyverse)
library(patchwork)

source("my_theme.R")

plot_nea_time <- function(df, snp_cutoff = 0, ymax = 0.06, title="") {
  filter(df,
         pop %in% c("EMH", "WestEurasia"),
         !X %in% near_east,
         snp_count > snp_cutoff) %>% 
    ggplot(aes(age, alpha)) +
    geom_point(aes(size = snp_count), alpha = 1/2) +
    geom_errorbar(aes(ymin = alpha - 2 * stderr, ymax = alpha + 2 * stderr), alpha = 1/5) +
    geom_smooth(aes(weight = snp_count), method = "lm", linetype = 2, fullrange = TRUE, size = 0.5, se = FALSE) +
    xlim(40000, 0) + ylim(0, ymax) +
    xlab("years before present") + ylab("Neandertal ancestry proportion") +
    ggtitle(title)
}

rdata_dir <- "rds"

# Near East populations to exclude from West Eurasians in plots over time
near_east <- c("BedouinB", "Druze", "Iranian", "Iraqi_Jew", "Jordanian",
               "Palestinian", "Samaritan", "Turkish", "Yemenite_Jew")


# Neandertal ancestry over time -------------------------------------------

# new direct f4 ratio estimates using Altai and Vindija
direct_f4_all <- readRDS(file.path(rdata_dir, "direct_f4_combined.rds")) %>% filter(X != "UstIshim")
direct_f4_bigyri <- readRDS(file.path(rdata_dir, "direct_f4_BigYRI.rds")) %>% filter(X != "UstIshim")
direct_f4_ho <- readRDS(file.path(rdata_dir, "direct_f4_HO.rds")) %>% filter(X != "UstIshim")

# indirect f4 ratio estimates based on assumptions about African populations
indirect_f4_all <- readRDS(file.path(rdata_dir, "indirect_f4_combined.rds")) %>% filter(X != "UstIshim")
indirect_f4_bigyri <- readRDS(file.path(rdata_dir, "indirect_f4_BigYRI.rds")) %>% filter(X != "UstIshim")
indirect_f4_ho <- readRDS(file.path(rdata_dir, "indirect_f4_HO.rds")) %>% filter(X != "UstIshim")

# admixture array estimates
admix_prop <- readRDS(file.path(rdata_dir, "admixture_array_prop.rds")) %>% filter(X != "UstIshim")

indirect_f4_all %>%
  plot_nea_time(snp_cutoff = 200000, title = "Ice Age indirect f4 ratio") +
  my_theme(legend.position = "none")

admix_prop %>%
  plot_nea_time(title = "Proportion of Neandertal-like alleles") +
  my_theme(legend.position = "none")

filter(direct_f4_all, C == "Dinka") %>%
  plot_nea_time(snp_cutoff = 200000, title = "Direct f4 ratio using Altai and Vindija") +
  my_theme(legend.position = "none")

filter(direct_f4_bigyri, C == "Dinka") %>%
  plot_nea_time(snp_cutoff = 50000, title = "Direct f4 ratio - Big Yoruba sites") +
  my_theme(legend.position = "none")

filter(direct_f4_ho, C == "Dinka") %>%
  plot_nea_time(snp_cutoff = 150000, title = "Direct f4 ratio - Human Origins sites") +
  my_theme(legend.position = "none")


# Affinities over time to different SGDP populations ----------------------

affinities <- readRDS(file.path(rdata_dir, "sgdp_affinities.rds"))

affinities %>% filter(str_detect(Y, "Africa"), Y != "NorthAfrica", n_snps > 50000) %>%
  mutate(Y = str_replace(Y, "Africa", " Africa")) %>% 
  ggplot(aes(age, D)) + geom_point(aes(shape = abs(Zscore) > 3), size = 5, alpha = 1/5) +
  geom_hline(yintercept = 0, linetype = 2, color = "blue") +
  xlim(40000, -2000) + ylim(-0.05, 0.05) +
  xlab("years before present") + ylab("D(Ust'-Ishim, West Eurasian; African population, Chimp)") +
  facet_grid(. ~ Y) +
  scale_shape_manual(values = c(21, 19)) +
  my_theme(legend.position = "none", axis.text.x = element_text(hjust = 1, angle = 30))

affinities %>% filter(Y == "Oceania") %>%
  ggplot(aes(age, Zscore)) + geom_point(aes(shape = abs(Zscore) > 3), size = 5, alpha = 1/2) +
  geom_hline(yintercept = 0, linetype = 2, color = "blue") +
  geom_hline(yintercept = c(-3, 3), linetype = 2, color = "red") +
  xlim(40000, -2000) + ylim(-20, 20) +
  xlab("years before present") + ylab("Z score for D(Ust'-Ishim, West Eurasian; Papuans, Chimp)") +
  facet_grid(. ~ Y) +
  scale_shape_manual(values = c(21, 19)) +
  my_theme(legend.position = "none", axis.text.x = element_text(hjust = 1, angle = 30))


# Deltas of allele frequencies --------------------------------------------

deltas <- readRDS(file.path(rdata_dir, "deltas.rds"))


# Deleterious delta plots -------------------------------------------------

s_breaks <- c(0, -3, -5, -Inf)
s_labels <- as.character(10^s_breaks)
s_bins <- rev(paste0("s = ", s_labels[-length(s_labels)], paste0(" - ", s_labels[-1])))
deltas %>% mutate(logS_bin=cut(logS,
                               breaks=s_breaks,
                               labels=s_bins)) %>%
  replace_na(list(logS_bin="s = 1e-05 - 0")) %>% 
  mutate(logS_bin=factor(logS_bin, levels=rev(s_bins))) %>%
  group_by(g, mut_type, rep, logS_bin) %>%    
  summarise(avg_delta = mean(delta_f), sd_delta = sd(delta_f), n_delta = n()) %>%
  mutate(se_delta = sd_delta / sqrt(n_delta),
         lower_ci = avg_delta - qt(1 - (0.05 / 2), n_delta - 1) * se_delta,
         upper_ci = avg_delta + qt(1 - (0.05 / 2), n_delta - 1) * se_delta) %>% 
  filter(g >= 100) %>% 
  filter(mut_type %in% c("MH_del", "Nea_del")) %>% 
  ungroup %>% 
  mutate(mut_type = str_replace(mut_type, "MH_del", "modern human") %>%
                    str_replace("Nea_del", "introgressed")) %>% 
  ggplot(aes(g, avg_delta)) +
  geom_line(aes(color = rep), linetype = 2) +
  geom_smooth(color="black", alpha = 1/2, size = 0.5) + 
  geom_hline(yintercept = 0, linetype = 2, size = 0.25) + facet_grid(mut_type ~ logS_bin) +
  xlab("generation since admixture") + ylab("average frequency change per generation") +
  my_theme(legend.position = "none")

# Neutral marker delta plots ----------------------------------------------

group_by(deltas, g, mut_type, rep) %>%    
  summarise(avg_delta=mean(delta_f), sd_delta=sd(delta_f), n_delta=n()) %>%
  mutate(se_delta=sd_delta / sqrt(n_delta),
         lower_ci=avg_delta - qt(1 - (0.05 / 2), n_delta - 1) * se_delta,
         upper_ci=avg_delta + qt(1 - (0.05 / 2), n_delta - 1) * se_delta) %>% 
  filter(g >= 100) %>% 
  ungroup %>% 
  filter(mut_type %in% c("gap_marker", "region_marker")) %>% 
  mutate(mut_type = str_replace(mut_type, "gap_marker", "neutral markers outside selected regions") %>% 
                    str_replace("region_marker", "neutral markers within selected regions")) %>% 
  ggplot(aes(g, avg_delta)) +
  geom_line(aes(color = rep), linetype = 2) +
  geom_smooth(color = "black", alpha = 1/2, size = 0.5) + 
  geom_hline(yintercept = 0, color = "blue", linetype = 2) +
  facet_grid(mut_type ~ .) + 
  xlab("generation since admixture") + ylab("average frequency change per generation") +
  my_theme(legend.position = "none")



muts <- readRDS("")