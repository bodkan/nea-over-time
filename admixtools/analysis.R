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
source("admixtools.R")
source("../R/utils.R")

fix_name <- function(str) {
    str_replace_all(str, "-|\\.", "_") %>%
        str_replace_all("^S_|_[1-9]*", "")
}


# some small-scale testing
create_f4_poplist(X=c("French", "Sardinian", "Han", "Dai", "Stuttgart", "Oase1", "UstIshim"),
                  A="Yoruba", B="Dinka", C="Altai", O="Chimp",
                  filename="test.pop")
create_f4_parfile("test.par", "UPA.K.P.V1.3.2", "test.pop")




# load names of SGDP West Eurasian subpopulations
west_eurasians <-
    load_sgdp_info("../raw_data/10_24_2014_SGDP_metainformation_update.txt") %>%
    filter(Region == "WestEurasia") %>%
    .[["name"]] %>%
    str_replace("^S_", "") %>%
    str_replace("_[0-9]+$", "") %>%
    unique

# load names of the Ice Age upper-paleolithic humans
emhs <- read.table("../raw_data/Fu/archaic.ind",
                   stringsAsFactors=FALSE) %>%
    filter(V3 != "Oase1") %>%
    .[["V3"]]


############################################################
# load sample ages
emh_ages <-
    read_delim("../clean_data/ages.txt", delim=" ") %>%
    mutate(name=fix_name(name))
sgdp_ages <- tibble(name=west_eurasians, age=0L)
ages <- bind_rows(emh_ages, sgdp_ages)


############################################################
# calculate Nea. ancestry proportions using f4
create_f4_poplist_file(X=c(west_eurasians, emhs),
                       A="Yoruba", B="Dinka", C="Altai", O="Chimp",
                       filename="nea_in_west_eurasians.pop")

create_f4_param_file("nea_in_west_eurasians.par",
                     "nea_in_west_eurasians.pop",
                     eigenstrat_prefix="UPA.K.P.V1.3.2")

run_f4("nea_in_west_eurasians.par", "nea_in_west_eurasians.log")

f4_res <- read_f4_ratios("nea_in_west_eurasians.log")
f4_nea <- select(f4_res, name=X, alpha) %>%
    mutate(alpha=1 - alpha,
           name=fix_name(name),
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
    filter(name %in% f4_nea$name) %>%
    mutate(method="direct")


nea_estimates <- bind_rows(inner_join(ages, f4_nea), inner_join(ages, direct_nea))

nea_cor <- spread(nea_estimates, method, alpha)
plot(nea_cor$direct, nea_cor$f4, xlim=c(0, 0.1), ylim=c(0, 0.1))
summary(lm(direct ~ f4, data=nea_cor))
