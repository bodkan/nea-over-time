# idea - previous evidence suggests that there has been a potential
# gene flow between some middle eastern population after at least ~14 years ago,
# possibly related to the original farmers
# however, this evidence is based on D statistics, which can reject a hypothesis
# of no admixture and indicate a direction of gene-flow, but it cannot tell
# us anything about the admixture proportions
# here we want to see how much potential admixture from a mysterious basal
# eurasian lineage was there
library(tidyverse)
library(stringr)

source("~/projects/slim-neanderthal/R/utils.R")
source("admixr.R")
source("../R/utils.R")

fix_name <- function(str) {
    str_replace_all(str, "-|\\.", "_") %>%
        str_replace_all("^S_|_[1-9]+", "")
}




######################################################################
# prepare the tables with the SGDP non-African populations (ages
# geographical coordinates, etc.)
sgdp <- load_sgdp_info("../raw_data/10_24_2014_SGDP_metainformation_update.txt") %>%
    filter(! Region %in% c("Africa", "Oceania")) %>%
    select(-Country, pop=Region) %>%
    mutate(age=0, name=fix_name(name)) %>%
    group_by(name, age, pop) %>%
    summarise(Latitude=mean(Latitude), Longitude=mean(Longitude)) %>%
    ungroup
    
emhs <- read_delim("../emh_ages.txt", delim=" ", col_names=c("name", "age")) %>%
    mutate(pop="EMH", Latitude=NA, Longitude=NA)

samples <- bind_rows(emhs, sgdp)


############################################################
# calculate Nea. ancestry proportions using f4
create_qpF4ratio_pops(X=samples$name,
                      A="Yoruba", B="Dinka", C="Altai", O="Chimp",
                      file="nea_estimate.pop")

# generate a table of SNPs to filter out (transitions)
read_fwf("UPA.K.P.V1.3.2.snp",
         fwf_widths(c(20, 6, 16, 16, 2, 2),
                    col_names=c("id", "chrom", "gen", "pos", "alt", "ref")),
         progress=FALSE) %>%
    keep_transitions %>%
    write_tsv("UPA.K.P.V1.3.2.transitions.snp", col_names=FALSE)

create_param_file("nea_estimate.par",
                  "nea_estimate.pop",
                  eigenstrat_prefix="UPA.K.P.V1.3.2",
                  badsnp_file="UPA.K.P.V1.3.2.transitions.snp")

run_cmd("qpF4ratio", param_file="nea_estimate.par", log_file="nea_estimate.log")

nea_f4_ratios <- read_qpF4ratio("nea_estimate.log")
nea_f4_df <- select(nea_f4_ratios, name=X, alpha) %>%
    mutate(alpha=1 - alpha,
           method="f4")


############################################################
# calculate "direct" estimate of Nea. ancestry
array_snps <- load_dataset("../clean_data/ice_age.tsv",
                           "../clean_data/sgdp.tsv",
                           "../clean_data/archaics.tsv",
                           filter_damage=TRUE,
                           random_sample=TRUE)

direct_nea <- select(array_snps, -c(chrom, pos, ref, alt, contains("archaic"))) %>%
    summarise_all(function(ind) { mean(ind, na.rm=TRUE) / 2 }) %>%
    gather(name, alpha) %>%
    mutate(name=fix_name(name)) %>%
    group_by(name) %>%
    summarise(alpha=mean(alpha)) %>%
    filter(name %in% nea_f4_df$name) %>%
    mutate(method="direct")


nea_long <- bind_rows(inner_join(samples, nea_f4_df), inner_join(samples, direct_nea)) %>% filter(name !=  "Oase1", alpha > 0, alpha < 0.08)
nea_wide <- spread(nea_long, method, alpha) %>% filter(!is.na(direct))

save.image("../RData/nea_estimate.RData")

############################################################
# plot the comparison of f4 and direct Nea. estimates

load("../RData/nea_estimate.RData")

pdf("f4_vs_direct_props.pdf")

# comparison of f4 vs direct Nea. estimates
ggplot(nea_wide, aes(f4, direct)) +
    geom_point() +
    geom_smooth(method="lm", linetype=2, size=0.5) +
    ggtitle("Correlation of Nea% estimates using f4 ratios or admixture array SNPs",
            "f4 ratio = 1 - f4(YRI, Chimp; X, Altai) / f4(YRI, Chimp; Dinka, Altai)")

ggplot(nea_wide, aes(f4, direct, color=pop, group=pop)) +
    geom_point() +
    geom_smooth(method="lm", linetype=2, size=0.5) +
    ggtitle("Correlation of Nea% estimates using f4 ratios or admixture array SNPs",
            "f4 ratio = 1 - f4(YRI, Chimp; X, Altai) / f4(YRI, Chimp; Dinka, Altai)")

ggplot(filter(nea_wide, pop != "EMH"), aes(f4, direct, color=pop, group=pop)) +
    geom_point() +
    geom_smooth(method="lm", linetype=2, size=0.5) +
    ggtitle("Correlation of Nea% estimates using f4 ratios or admixture array SNPs",
            "f4 ratio = 1 - f4(YRI, Chimp; X, Altai) / f4(YRI, Chimp; Dinka, Altai)")

# Nea over time in West Eurasians
filter(nea_long, pop %in% c("EMH", "WestEurasia")) %>%
ggplot(aes(age, alpha)) +
    geom_point() +
    geom_smooth(method="lm", linetype=2, fullrange=TRUE, size=0.5) +
    facet_grid(method ~ .) +
    xlim(55000, 0) + ylim(0, 0.1) +
    ggtitle("Nea% (alpha) over time using Ice Age paper EMHs and SGDP West Eurasians",
            "Upper panel - Nea. estimates using admixture array SNPs, bottom panel - ratio of f4 statistics")

dev.off()





######################################################################
# test for the presence of basal Eurasian ancestry, following the
# approach outlined in Lazaridis et al. 2016 (supplementary
# material page 18)





############################################################
# estimating the basal Eurasian ancestry, following the
# approach outlined in Lazaridis et al. 2016 (supplementary
# material page 23)

