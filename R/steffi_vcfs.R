#!/usr/bin/env Rscript
library(tidyverse)
library(stringr)
library(admixr)
devtools::reload(devtools::inst("admixr"))

source("R/utils.R")

setwd("raw_data/merged_vcfs")

snp_all <- data.frame()
geno_all <- c()

for (chrom in as.character(1:22)) {
    vcf_file <- paste0("merged_arch_SGDP_chr", chrom, ".almostvcf.gz")

    cat("Processing file: ", vcf_file, "...")

    # load the pseudo VCF file
    chrom_vcf <- read_tsv(vcf_file) %>%
        mutate(PanTro4=toupper(PanTro4), PonAbe2=toupper(PonAbe2)) %>%
        mutate(Chimp=as.integer(PanTro4 != REF), Orang=as.integer(PonAbe2 != REF)) %>%
        rename(chrom=`#CHROM`, pos=POS, ref=REF, alt=ALT, new_Altai=AltaiNeandertal,
               new_Vindija=Vindija) %>%
        select(-c(ID, QUAL, FILTER, INFO, FORMAT, PanTro4, PonAbe2))

    # generate dataframes with the 3 EIGENSTRAT info tables
    snp <- select(chrom_vcf, chrom, pos, ref, alt) %>%
        mutate(snp_id=paste(chrom, pos, sep="_"), gen_dist="0.0") %>%
        select(snp_id, chrom, gen_dist, pos, ref, alt)
    ind <- tibble(
        sample_id=select(chrom_vcf, -c(chrom, pos, ref, alt)) %>% names,
        sex="U",
        label=fix_name(sample_id)
    )
    geno <- select(chrom_vcf, -c(chrom, pos, ref, alt)) %>%
        mutate_all(gt_to_eigenstrat)

    # merge with the data processed so far
    snp_all <- bind_rows(snp_all, snp)
    geno_all <- bind_rows(geno_all, geno)

    cat("done.\n")
}


write_tsv(snp_all, "steffi.snp", col_names=FALSE)
write_tsv(ind, "steffi.ind", col_names=FALSE)
writeLines(apply(geno_all, 1, paste, collapse=""), "steffi.geno")







## i <- read_ind("../eigenstrat_all/UPA_merged_all.ind")
## g <- read_geno("../eigenstrat_all/UPA_merged_all.geno")
## names(g) <- i$id
## s <- read_snp("../eigenstrat_all/UPA_merged_all.snp")


## geno <- read_geno("steffi.geno")
## ind <- read_ind("steffi.ind")
## names(geno) <- ind$id
## snp <- read_snp("steffi.snp")

## qiaomei <- bind_cols(s, g) %>% select(chrom, pos, ref, alt, Ust_Ishim, new_Altai, new_Vindija, Chimp)
## steffi <- bind_cols(snp, geno) %>% select(chrom, pos, ref, alt, Ust_Ishim, new_Altai, new_Vindija, Chimp)


## x <- inner_join(qiaomei, steffi, by=c("chrom", "pos"))


## select(x, contains("Altai")) %>% table
## select(x, contains("Vindija")) %>% table
## select(x, contains("Ust_Ishim")) %>% table

## select(x, contains("Chimp")) %>% table
