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
library(stringr)

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
                                                 "ref", "alt")))
ind <- read_table(ind_file, col_names=c("id", "sex", "label"))
geno <- readLines(geno_file) %>% str_split("", simplify=TRUE)

# replace dashes and dots in sample names with underscores
colnames(geno) <- ind$id %>% str_replace("-|\\.", "_")

snp %>%
    select(chrom, pos, ref, alt) %>% # select only columns of interest
    bind_cols(as_tibble(geno)) %>%   # merge with the snp calls
    write_tsv(output_file)


x <- read_tsv("../tmp/t.tsv")
y <- read_tsv("../tmp/i.tsv")
colnames(y) <- str_replace(colnames(y), "Dolni", "Vestonice")

ids <- colnames(x)[5:length(colnames(x))]

xi <- sapply(ids, function(i) table(x[i]))
yi <- sapply(ids, function(i) table(y[i]))
