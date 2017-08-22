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
    mutate(pop="EMH", Latitude=NA, Longitude=NA) %>%
    filter(name != "Oase1")

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

############################################################
# plot the comparison of f4 and direct Nea. estimates


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
## test for the presence of basal Eurasian ancestry, following the
## approach outlined in Lazaridis et al. 2016 (supplementary
## material page 18)

## To test for Basal Eurasian ancestry, we computed the statistic
## f4(Test, Han; Ust’- Ishim, Chimp) (Supplementary Information,
## section 4), which measures the excess of allele sharing of
## Ust’-Ishim with a variety of Test popula- tions compared to Han as
## a baseline. This statistic is significantly negative (Z < −3.7) for
## all ancient Near Easterners as well as Neolithic and later
## Europeans, consistent with them having ancestry from a deeply
## divergent Eurasian lineage that separated from the ancestors of
## most Eurasians before the separation of Han and Ust’-Ishim.
## (main text of the paper)

## basal eurasian ancestry evidence in ancient near east + ice age samples
create_Dstats_pops(
    W=c("Levant_N", "Natufian", "Seh_Gabi_Iran_Chalcolithic", "Iran_N", "Anatolia_N", emhs$name),
    X=c("Han", "Onge", "Karitiana"),
    Y="UstIshim",
    Z="Chimp",
    file="basal_ancient_f4.pop"
)
create_param_file(
    param_file="basal_ancient_f4.par",
    pops_file="basal_ancient_f4.pop",
    eigenstrat_prefix="UPA.K.P.V1.3.2",
    badsnp_file="UPA.K.P.V1.3.2.transitions.snp",
    f4mode=TRUE
)
run_cmd("qpDstat", param_file="basal_ancient_f4.par", log_file="basal_ancient_f4.log")

# basal eurasian in modern West Eurasians
create_Dstats_pops(
    W=filter(sgdp, pop == "WestEurasia")$name,
    X=c("Han", "Onge", "Karitiana"),
    Y="UstIshim",
    Z="Chimp",
    file="basal_modern.pop"
)
create_param_file(
    param_file="basal_modern.par",
    pops_file="basal_modern.pop",
    eigenstrat_prefix="UPA.K.P.V1.3.2",
    badsnp_file="UPA.K.P.V1.3.2.transitions.snp",
    f4mode=TRUE
)
run_cmd("qpDstat", param_file="basal_modern.par", log_file="basal_modern.log")


pdf("f4_basal_eurasian_zscores.pdf", width=10, height=4)

basal_ancient_f4 <- read_qpDstat("basal_ancient.log")

# samples by zscore
group_by(basal_ancient_f4, W) %>%
    summarise(Dstat=mean(Dstat), Zscore=mean(Zscore)) %>%
    mutate(W=factor(W, levels=W[order(Zscore)])) %>%
    ggplot(aes(W, Zscore)) +
    geom_point() +
    geom_hline(yintercept=c(1, -1) * qnorm(0.975, 0, 1), linetype=2) +
    coord_cartesian(ylim=c(2, -10)) +
    ggtitle("Ancient samples - Z-scores for f4 detecting basal Eurasian ancestry",
            "Dashed line indicates 5% significance cutoff (what is bellow the line shows evidence of basal Eurasian ancestry)") +
    theme(axis.text.x=element_text(angle=35, hjust=1, size=8))

basal_modern_f4 <- read_qpDstat("basal_modern.log")

group_by(basal_modern_f4, W) %>%
    summarise(Dstat=mean(Dstat), Zscore=mean(Zscore)) %>%
    mutate(W=factor(W, levels=W[order(Zscore)])) %>%
    ggplot(aes(W, Zscore)) +
    geom_point() +
    geom_hline(yintercept=c(1, -1) * qnorm(0.975, 0, 1), linetype=2) +
    coord_cartesian(ylim=c(2, -10)) +
    ggtitle("Modern samples - Z-scores for f4 detecting basal Eurasian ancestry",
            "Dashed line indicates 5% significance cutoff (what is bellow the line shows evidence of basal Eurasian ancestry)") +
    theme(axis.text.x=element_text(angle=35, hjust=1, size=8))


inner_join(
    group_by(filter(bind_rows(basal_ancient_f4, basal_modern_f4), W %in% samples$name),
             W) %>% rename(name=W) %>% summarise(Zscore=mean(Zscore)),
    samples
) %>%
    ggplot(aes(age, Zscore)) +
    geom_point() +
    geom_hline(yintercept=c(1, -1) * qnorm(0.975, 0, 1), linetype=2) +
    coord_cartesian(ylim=c(2, -10)) +
    ggtitle("Basal Eurasian ancestry significance level vs sample ages",
            "Dashed line indicates 5% significance cutoff (what is bellow the line shows evidence of basal Eurasian ancestry)") +
    xlim(55000, 0)


dev.off()













############################################################
# estimating the basal Eurasian ancestry, following the
# approach outlined in Lazaridis et al. 2016 (supplementary
# material page 23)


pdf("f4_basal_eurasian_ratios.pdf", width=8, height=5)


create_qpF4ratio_pops(X=filter(samples, pop %in% c("EMH", "WestEurasia"),
                               !name %in% c("UstIshim", "Kostenki14"))$name,
                      A="UstIshim", B="Kostenki14", C="Mbuti", O="Mbuti",
                      file="basal_ratios.pop")

create_param_file("basal_ratios.par",
                  "basal_ratios.pop",
                  eigenstrat_prefix="UPA.K.P.V1.3.2",
                  badsnp_file="UPA.K.P.V1.3.2.transitions.snp")

run_cmd("qpF4ratio", param_file="basal_ratios.par", log_file="basal_ratios.log")

nea_f4_ratios <- read_qpF4ratio("basal_ratios.log")
nea_f4_ratios

nea_f4_df <- select(nea_f4_ratios, name=X, alpha) %>%
    mutate(alpha=1 - alpha,
           method="f4")



## F4 ratio estimates using raw Dstats/f4

## denominator
create_Dstats_pops(
    W="Mbuti",
    X="Loschbour",
    Y="UstIshim",
    Z="Kostenki14",
    file="blah1.pop"
)
create_param_file(
    param_file="blah1.par",
    pops_file="blah1.pop",
    eigenstrat_prefix="UPA.K.P.V1.3.2",
    badsnp_file="UPA.K.P.V1.3.2.transitions.snp",
    f4mode=TRUE
)
run_cmd("qpDstat", param_file="blah1.par", log_file="blah1.log")

# numerator - ancient pops
create_Dstats_pops(
    W=c("Levant_N", "Natufian", "Iran_N", "Anatolia_N",
        filter(samples, pop == "EMH",
               !name %in% c("UstIshim", "Kostenki14", "Loschbour", "Oase1"))$name),
    X="Loschbour",
    Y="UstIshim",
    Z="Kostenki14",
    file="blah2_ancient.pop"
)
create_param_file(
    param_file="blah2_ancient.par",
    pops_file="blah2_ancient.pop",
    eigenstrat_prefix="UPA.K.P.V1.3.2",
    badsnp_file="UPA.K.P.V1.3.2.transitions.snp",
    f4mode=TRUE
)
run_cmd("qpDstat", param_file="blah2_ancient.par", log_file="blah2_ancient.log")

# numerator - modern pops
create_Dstats_pops(
    W=filter(samples, pop == "WestEurasia")$name,
    X="Loschbour",
    Y="UstIshim",
    Z="Kostenki14",
    file="blah2_modern.pop"
)
create_param_file(
    param_file="blah2_modern.par",
    pops_file="blah2_modern.pop",
    eigenstrat_prefix="UPA.K.P.V1.3.2",
    badsnp_file="UPA.K.P.V1.3.2.transitions.snp",
    f4mode=TRUE
)
run_cmd("qpDstat", param_file="blah2_modern.par", log_file="blah2_modern.log")


blah1 <- read_qpDstat("blah1.log")
blah2 <- bind_rows(read_qpDstat("blah2_ancient.log"),
                   read_qpDstat("blah2_modern.log"))

## calculate ratios blah2 / blah1
basal_ratios <- mutate(blah2, alpha=Dstat/blah1$Dstat) %>%
    select(name=W, alpha) %>%
    left_join(samples) %>%
    mutate(name=factor(name, levels=name[order(alpha, decreasing=TRUE)]))

## basal Eurasian proportion over time
ggplot(basal_ratios, aes(age, alpha)) +
    geom_point() +
    ggtitle("Basal Eurasian ancestry proportion over time") +
    geom_smooth(method="lm") +
    xlim(55000, 0)

## basal Eurasian proportions in West Eurasia
filter(basal_ratios, pop == "WestEurasia") %>%
    ggplot(aes(name, alpha)) +
    geom_bar(stat="identity") +
    theme(axis.text.x=element_text(angle=45, hjust=1))


basal_vs_nea <- inner_join(
    rename(nea_long, alpha_nea=alpha),
    rename(basal_ratios, alpha_basal=alpha)
)
    

ggplot(filter(basal_vs_nea, pop == "WestEurasia", method == "direct"),
       aes(alpha_nea, alpha_basal)) +
    geom_point() +
    geom_smooth(method="lm")


dev.off()



