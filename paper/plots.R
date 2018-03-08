# install.packages("devtools")
# devtools::install_github("thomasp85/patchwork")
# devtools::install_github('cttobin/ggthemr')

library(tidyverse)
library(patchwork)

source("paper_theme.R")









# Neandertal ancestry over time -------------------------------------------

plot_nea_time <- function(df, snp_cutoff = 0, ymax = 0.06, title="") {
  filter(df,
         pop %in% c("EMH", "WestEurasia"),
         !X %in% near_east,
         snp_count > snp_cutoff) %>% 
    ggplot(aes(age, alpha)) +
    geom_point(aes(size = snp_count), alpha = 1/2) +
    geom_errorbar(aes(ymin = alpha - 2 * stderr, ymax = alpha + 2 * stderr), alpha = 1/5) +
    geom_smooth(aes(weight = snp_count), method = "lm", linetype = 2, fullrange = TRUE, size = 0.5) +
    xlim(40000, 0) + ylim(0, ymax) +
    xlab("years before present") + ylab("Neandertal ancestry proportion") +
    ggtitle(title)
}

# Near East populations to exclude from West Eurasians in plots over time
near_east <- c("BedouinB", "Druze", "Iranian", "Iraqi_Jew", "Jordanian",
               "Palestinian", "Samaritan", "Turkish", "Yemenite_Jew")

ignore_samples <- c(near_east, "UstIshim")

# new direct f4 ratio estimates using Altai and Vindija
direct_f4_all <- readRDS("rds/direct_f4_combined.rds") %>% filter(!X %in% ignore_samples)
direct_f4_bigyri <- readRDS("rds/direct_f4_BigYRI.rds") %>% filter(!X %in% ignore_samples)
direct_f4_ho <- readRDS("rds/direct_f4_HO.rds") %>% filter(!X %in% ignore_samples)

# indirect f4 ratio estimates based on assumptions about African populations
indirect_f4_all <- readRDS("rds/indirect_f4_combined.rds") %>% filter(!X %in% ignore_samples)
indirect_f4_bigyri <- readRDS("rds/indirect_f4_BigYRI.rds") %>% filter(!X %in% ignore_samples)
indirect_f4_ho <- readRDS("rds/indirect_f4_HO.rds") %>% filter(!X %in% ignore_samples)

# admixture array estimates
admix_prop <- readRDS("rds/admixture_array_prop.rds") %>% filter(!X %in% ignore_samples)

# plots

indirect_f4_all %>%
  plot_nea_time(snp_cutoff = 200000, title = "Ice Age indirect f4 ratio") +
  paper_theme(legend.position = "none")
ggsave("fig/indirect_f4.png")

admix_prop %>%
  plot_nea_time(title = "Proportion of Neandertal-like alleles") +
  paper_theme(legend.position = "none")
ggsave("fig/admix_prop.png")

filter(direct_f4_all, C == "Dinka") %>%
  plot_nea_time(snp_cutoff = 200000, title = "Direct f4 ratio") +
  paper_theme(legend.position = "none")
ggsave("fig/direct_f4_all.png")

filter(direct_f4_bigyri, C == "Dinka") %>%
  plot_nea_time(snp_cutoff = 50000, title = "Direct f4 ratio (Big Yoruba sites)") +
  paper_theme(legend.position = "none")
ggsave("fig/direct_f4_bigyri.png")

filter(direct_f4_ho, C == "Dinka") %>%
  plot_nea_time(snp_cutoff = 150000, title = "Direct f4 ratio (Human Origins sites)") +
  paper_theme(legend.position = "none")
ggsave("fig/direct_f4_ho.png")









# Affinities over time to different SGDP populations ----------------------

afr_affinities <- readRDS("rds/african_affinities.rds")

afr_affinities %>% filter(Y %in% c("Mbuti", "Dinka", "Yoruba"), n_snps > 50000) %>%
  mutate(Y = paste0("Y = ", Y)) %>% 
  mutate(Y = str_replace(Y, "Africa", " Africa")) %>% 
  ggplot(aes(age, D)) + geom_point(aes(shape = abs(Zscore) > 3), size = 5, alpha = 1/5) +
  geom_hline(yintercept = 0, linetype = 2, color = "blue") +
  xlim(40000, -2000) + ylim(-0.05, 0.02) +
  xlab("years before present") + ylab("D statistic") +
  ggtitle("D(Ust'-Ishim, West Eurasian; African Y, Chimp)") +
  facet_grid(. ~ Y) +
  scale_shape_manual(values = c(21, 19)) +
  paper_theme(legend.position = "none", axis.text.x = element_text(hjust = 1, angle = 30))
ggsave("fig/affinities_africa.png")

afr_affinities %>% filter(Y == "Oceania") %>%
  ggplot(aes(age, Zscore)) + geom_point(aes(shape = abs(Zscore) > 3), size = 5, alpha = 1/2) +
  geom_hline(yintercept = 0, linetype = 2, color = "blue") +
  geom_hline(yintercept = c(-3, 3), linetype = 2, color = "red") +
  xlim(40000, -2000) + ylim(-20, 20) +
  xlab("years before present") + ylab("Z score") +
  ggtitle("D(Ust'-Ishim, West Eurasian; Oceanians, Chimp)") +
  scale_shape_manual(values = c(21, 19)) +
  paper_theme(legend.position = "none", axis.text.x = element_text(hjust = 1, angle = 30))
ggsave("fig/affinities_oceania.png")










# Deltas of allele frequencies --------------------------------------------

deltas <- readRDS("rds/deltas.rds")


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
  geom_line(aes(color = rep), size = 1, alpha = 1/2) +
  geom_smooth(color = "black", alpha = 1/2) + 
  geom_hline(yintercept = 0, linetype = 2, size = 0.25) + facet_grid(mut_type ~ logS_bin) +
  xlab("generation since admixture") + ylab("average frequency change") +
  paper_theme(legend.position = "none", axis.text.x = element_text(hjust = 1, angle = 30))
ggsave("fig/deltas_deleterious.png")

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
  geom_line(aes(color = rep), size = 1, alpha = 1/2) +
  geom_smooth(color = "black", alpha = 1/2) + 
  geom_hline(yintercept = 0, color = "blue", linetype = 2) +
  facet_grid(. ~ mut_type) + 
  xlab("generation since admixture") +
  ylab("average frequency change per generation") +
  paper_theme(legend.position = "none", axis.text.x = element_text(hjust = 1, angle = 30))
ggsave("fig/deltas_neutral.png")









# Frequency vs selection coefficient --------------------------------------

# drift-selection "boundary"
# http://people.bu.edu/msoren/BI5152013/PGLecture27.pdf
# https://math.la.asu.edu/~jtaylor/teaching/Fall09/APM541/selection.pdf

muts <- readRDS("rds/first_gen_mutations.rds")

mutate(muts, mut_type = str_replace(mut_type, "MH_del", "modern human") %>%
                    str_replace("Nea_del", "introgressed")) %>% 
  group_by(mut_type) %>%
  mutate(freq = freq / max(freq)) %>%
  ggplot(aes(log10(-S), freq, color = (freq == 1), group = mut_type)) +
  geom_point(size = 3, alpha = 1/2) +
  ylab("allele frequency") +
  facet_grid(mut_type ~ .) +
  geom_vline(xintercept = log10(1 / 1000), linetype = 2) +
  geom_vline(xintercept = log10(1 / 10000), linetype = 2) +
  scale_x_reverse(name = "selection coefficient",
                  breaks = log10(1 / c(1000, 10000)),
                  labels = c(1/1000, 1/10000)) +
  scale_color_manual(values = c("gray", "black")) +
  paper_theme(legend.position = "none",
              axis.text.x=element_text(angle = 45, hjust = 1))
