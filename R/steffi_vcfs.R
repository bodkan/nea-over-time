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
    vcf_file <- paste0("merged_arch_SGDP_chr", chrom, ".vcf.gz")

    cat("Processing file: ", vcf_file, "...")

    # load the pseudo VCF file
    chrom_vcf <- read_tsv(vcf_file) %>%
        rename(chrom=`#CHROM`, pos=POS, ref=REF, alt=ALT, new_Altai=AltaiNeandertal,
               new_Vindija=Vindija) %>%
        select(-c(ID, QUAL, FILTER, INFO, FORMAT))

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







## qiaomei_ind <- read_ind("../eigenstrat_all/UPA_merged_all.ind")
## qiaomei_geno <- read_geno("../eigenstrat_all/UPA_merged_all.geno")
## names(qiaomei_geno) <- qiaomei_ind$id
## qiaomei_snp <- read_snp("../eigenstrat_all/UPA_merged_all.snp") %>% select(chrom, pos)
## qiaomei_pos <- select(qiaomei_snp, chrom, pos)

## steffi_geno <- read_geno("steffi.geno")
## steffi_ind <- read_ind("steffi.ind")
## names(steffi_geno) <- steffi_ind$id
## steffi_snp <- read_snp("steffi.snp") %>% select(chrom, pos)
## steffi_pos <- select(steffi_snp, chrom, pos)

## shared_samples <- intersect(names(qiaomei_geno), names(steffi_geno))

## steffi_geno <- select(steffi_geno, shared_samples)
## qiaomei_geno <- select(qiaomei_geno, shared_samples)

## steffi <- bind_cols(steffi_snp, steffi_geno) %>% inner_join(qiaomei_snp)
## qiaomei <- bind_cols(qiaomei_snp, qiaomei_geno) %>% inner_join(select(steffi, chrom, pos))


## table(steffi$new_Altai)
## table(qiaomei$new_Altai)

## table(steffi$new_Vindija)
## table(qiaomei$new_Vindija)

## table(steffi$`S_French-1`)
## table(qiaomei$`S_French-1`)

## table(steffi$`S_French-1`)/sum(table(steffi$`S_French-1`))
## table(qiaomei$`S_French-1`)/sum(table(qiaomei$`S_French-1`))

