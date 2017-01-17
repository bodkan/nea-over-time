#!/usr/bin/env Rscript

#
# Process data in the EIGENSTRAT format, stored in files specified in
# the command line arguments:
#
# - genotype file
# - snp file
# - indiv file.
#
# Convert these into a single table in a TSV format (filename given as
# a fourth command line argument) that's easier to process with something
# else than ADMIXTOOLS.
#
# Detailed description of the EIGENSTRAT format can be found here:
# https://github.com/DReichLab/AdmixTools/tree/master/convertf


library(tidyverse)
library(magrittr)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg19)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
    stop("Genotype file, snp file, indiv. file and output file have to be specified.")
}

geno_file <- "archaic.geno" #args[1]
snp_file <- "archaic.snp" #args[2]
ind_file <- "archaic.ind" #args[3]
output_file <- "t.tsv" #args[4]

snp <- read_fwf(snp_file, fwf_widths(c(20, 6, 16, 16, 2, 2),
                                     col_names=c("id", "chrom", "gen", "pos",
                                                 "alt", "ref")))
ind <- read_table(ind_file, col_names=c("id", "sex", "label"))
geno <- readLines(geno_file) %>% str_split("", simplify=TRUE)

# replace dashes and dots in sample names with underscores
colnames(geno) <- ind$id %>% str_replace("-|\\.", "_")

# perform random calling on het calls - for each het call (1), flip a
# coin and replace it with 0 or 2
hets <- geno == 1
geno <- replace(geno, hets, 2 * rbinom(sum(hets), 1, 0.5))

# merge the geno/snp file data into a single dataframe
snp_table <-
    snp %>%
    select(chrom, pos, ref, alt) %>% # select only columns of interest
    bind_cols(as_tibble(geno))       # merge with the snp calls

# fetch hg19 reference alleles at positions from the array
hg19_refs <-
    mutate(snp_table, start=pos, end=start, chrom=paste0("chr", chrom)) %>%
    makeGRangesFromDataFrame %>% # convert dataframe to GRanges
    getSeq(Hsapiens, .) %>%      # fetch hg19 alleles at sites in the GRanges object
    as.vector

# select only those SNPs from the EIGENSTRAT data where the true
# hg19 allele matches the reference in the data (some 1% of sites
# seem to have ref/alt columns flipped)
snp_table %<>% filter(ref == hg19_refs)

write_tsv(snp_table, output_file)
