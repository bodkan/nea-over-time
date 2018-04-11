# install.packages("devtools")
# devtools::install_github("thomasp85/patchwork")
# devtools::install_github('cttobin/ggthemr')

library(tidyverse)
library(patchwork)

source("paper_theme.R")




geneflow <- bind_rows(
  readRDS("rds/grid_scrm_afr-eur.rds") %>% mutate(geneflow = "AFR-EUR"),
  readRDS("rds/grid_scrm_eur-afr.rds") %>% mutate(geneflow = "EUR-AFR"),
  readRDS("rds/grid_scrm_bidirect.rds") %>% mutate(geneflow = "both")
)

filter(geneflow, Ne == 10000, m == 0.0001, t == 15000, sites != "ho") %>%
  mutate(sites = case_when(sites == "admix" ~ "admixture array",
                           sites == "all" ~ "all SNPs",
                           sites == "bigyri" ~ "Big Yoruba")) %>% 
  ggplot(aes(age, nea)) +
  geom_point(alpha = 1/2) +
  geom_smooth(method = "lm", linetype = 2, color = "red") +
  xlim(40000, 0) + ylim(0, 0.05) +
  xlab("time before present [years]") + ylab("Neanderthal proportion") +
  facet_grid(geneflow ~ sites) +
  paper_theme(axis.text.x = element_text(hjust = 1, angle = 30))




# Neandertal ancestry over time -------------------------------------------

plot_nea_time <- function(df, statistic, array, snp_cutoff = 0, ymax = 0.06, title="") {
  filter(df,
         stat == statistic,
         sites == array,
         pop %in% c("EMH", "WestEurasia"),
         snp_count > snp_cutoff) %>% 
    ggplot(aes(age, alpha)) +
    geom_point(aes(size = snp_count), alpha = 1/2) +
    geom_errorbar(aes(ymin = alpha - 2 * stderr, ymax = alpha + 2 * stderr), alpha = 1/5) +
    geom_smooth(aes(weight = snp_count), method = "lm", linetype = 2, fullrange = TRUE, size = 0.5) +
    xlim(40000, 0) + ylim(0, ymax) +
    xlab("years before present") + ylab("Neandertal ancestry proportion") +
    ggtitle(title) +
    paper_theme(legend.position = "none")
}

# Near East populations to exclude from West Eurasians in plots over time
near_east <- c("BedouinB", "Druze", "Iranian", "Iraqi_Jew", "Jordanian",
               "Palestinian", "Samaritan", "Turkish", "Yemenite_Jew")
ignore_samples <- c(near_east)

# load the estimates
admix_array <- readRDS("rds/admixture_array_prop.rds") %>% filter(!X %in% ignore_samples) %>% mutate(stat = "admix_prop")
f4_nea <- readRDS("rds/f4_nea.rds") %>% filter(!X %in% ignore_samples)

nea_est <- bind_rows(admix_array, f4_nea) %>%
  select(X, stat, alpha, Zscore, snp_count, age, pop, stderr, sites)


plot_nea_time(nea_est, statistic = "indirect_f4", array = "all",
              snp_cutoff = 1000, title = "Direct f4 ratio")
ggsave("fig/indirect_f4.png")

plot_nea_time(nea_est, statistic = "direct_f4", array = "all",
              snp_cutoff = 1000, title = "Direct f4 ratio")
ggsave("fig/direct_f4.png")

plot_nea_time(nea_est, statistic = "admix_prop", array = "admixture_array",
              snp_cutoff = 1000, title = "Admixture array proportion")
ggsave("fig/admix_array_prop.png")







# Affinities over time to different SGDP populations ----------------------




# Grouped African populations ---------------------------------------------

sgdp_affinities <- readRDS("rds/sgdp_affinities_ui_chimp.rds")

sgdp_affinities %>% filter(Y %in% c("EastAfrica", "CentralAfrica", "WestAfrica"), n_snps > 50000) %>%
  mutate(Y = paste0("Y = ", Y)) %>% 
  mutate(Y = str_replace(Y, "Africa", " Afr.")) %>% 
  ggplot(aes(age, D)) + geom_point(aes(shape = abs(Zscore) > 3), size = 5, alpha = 1/5) +
  geom_hline(yintercept = 0, linetype = 2, color = "blue") +
  xlim(40000, -2000) + ylim(-0.05, 0.02) +
  xlab("years before present") + ylab("D statistic") +
  ggtitle("D(Ust'-Ishim, West Eurasian; African Y, Chimp)") +
  facet_grid(. ~ Y) +
  scale_shape_manual(values = c(21, 19)) +
  paper_theme(legend.position = "none", axis.text.x = element_text(hjust = 1, angle = 30))
ggsave("fig/sgdp_affinities_ui_chimp.png")

sgdp_affinities %>% filter(Y == "Oceania") %>%
  ggplot(aes(age, Zscore)) + geom_point(aes(shape = abs(Zscore) > 3), size = 5, alpha = 1/2) +
  geom_hline(yintercept = 0, linetype = 2, color = "blue") +
  geom_hline(yintercept = c(-3, 3), linetype = 2, color = "red") +
  xlim(40000, -2000) + ylim(-20, 20) +
  xlab("years before present") + ylab("Z score") +
  ggtitle("D(Ust'-Ishim, West Eurasian; Oceanians, Chimp)") +
  scale_shape_manual(values = c(21, 19)) +
  paper_theme(legend.position = "none", axis.text.x = element_text(hjust = 1, angle = 30))
ggsave("fig/african_affinities_ui_chimp_oceania.png")



# Individual Yoruba, Dinka and Mbuti --------------------------------------

afr_affinities <- readRDS("rds/african_affinities_ui_chimp.rds")

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
ggsave("fig/african_affinities_ui_chimp.png")





near_east_proportion <- readRDS("rds/nea_east_proportions.rds")

near_east_proportion %>% filter(B == "Levant_N") %>%
  filter(snp_count > 200000) %>% 
  ggplot(aes(age, alpha)) +
  geom_point(size=3, alpha=1/2) +
  geom_errorbar(aes(ymin=alpha - 2 * stderr, ymax=alpha + 2 * stderr)) +
  geom_hline(yintercept=c(0, 1), linetype=2, color="red") +
  facet_wrap(~ A) +
  xlim(40000, 0) + labs(y="proportion of ancient Near East ancestry") +
  theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1)) +
  ggtitle("Near East ancestry in Ice Age samples and SGDP Europeans over time",
          subtitle="Calculated on the 2.2M sites array, as a ratio of f4 statistics from the equation above")


near_east_proportion %>% filter(B == "Levant_N", A == "BedouinB") %>%
  filter(snp_count > 200000) %>%
  ggplot(aes(age, alpha)) +
  geom_point(size=3, alpha=1/2) +
  geom_errorbar(aes(ymin=alpha - 3 * stderr, ymax=alpha + 3 * stderr)) +
  geom_hline(yintercept=c(0, 1), linetype=2, color="red") +
  xlim(40000, 0) + labs(y="proportion of ancient Near East ancestry") +
  ggtitle("Proportion of ancient Near East ancestry in Europeans over time") +
  theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))




bin_props <- readRDS("rds/depletion_near_genes.rds")

bin_props %>% filter(gen == 2200) %>% 
  ggplot(aes(factor(dist_bin), nea_prop)) +
  geom_jitter(alpha = 1/5, size = 0.2) +
  geom_boxplot(alpha = 0.8) +
  geom_smooth(aes(group = 1), method = "lm", linetype = 2, color = "red") +
  paper_theme()






# Deltas of allele frequencies --------------------------------------------

deltas <- readRDS("rds/deltas.rds")


mutations <- readRDS("rds/mutations.rds")

mutations %>%
  filter(mut_type %in% c("gap_marker", "region_marker")) %>% 
  group_by(gen, rep, mut_type) %>%
  summarise(avg_nea=mean(freq)) %>%
  group_by(gen, mut_type) %>%
  summarise(mean_rep=mean(avg_nea), sd_rep=sd(avg_nea), n_rep=n()) %>%
  mutate(se_rep=sd_rep / sqrt(n_rep),
         lower_ci=mean_rep - qt(1 - (0.05 / 2), n_rep - 1) * se_rep,
         upper_ci=mean_rep + qt(1 - (0.05 / 2), n_rep - 1) * se_rep,
         mut_type=ifelse(mut_type == "gap_marker",
                         "markers outside selected regions",
                         "markers within selected regions")) %>%
ggplot(aes(gen, mean_rep, color=mut_type, group=mut_type)) +
  geom_line() +
  geom_ribbon(aes(ymin=lower_ci, ymax=upper_ci, fill=mut_type), alpha=1/2, color = NA) +
  geom_smooth(method = "lm",
              data=filter(nea_est, stat == "indirect_f4",sites == "all", pop %in% c("EMH", "WestEurasia")) %>%
                mutate(gen = (55000 - age) / 25),
              aes(gen, alpha), alpha = 1/5, inherit.aes = FALSE, se = FALSE, linetype = 2, color = "black") +
  geom_smooth(method = "lm",
              data=filter(nea_est, stat == "direct_f4",sites == "all", pop %in% c("EMH", "WestEurasia")) %>%
                mutate(gen = (55000 - age) / 25),
              aes(gen, alpha), inherit.aes = FALSE, se = FALSE, linetype = 4, color = "black") +
  xlab("generations after admixture") + ylab("mean frequency") +
  coord_cartesian(y=c(0, 0.1)) + labs(fill="", color="") +
  theme(legend.position = "bottom")






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
