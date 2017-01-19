#!/usr/bin/env Rscript

#
# Merge the processed TSV data from the Ice Age paper with the table
# of genotype calls (0/1/2) of the archaics (Altai, Vindija and
# Denisovan).
#

library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
    stop("Table with EMH data, table of archaic GTs and output file must be specified")
}

iceage_file   <- args[1]
archaics_file <- args[2]
output_file   <- args[3]

iceage <- read_tsv(iceage_file)
archaics <- read_tsv(archaics_file)

inner_join(iceage, archaics, by=c("chrom", "pos", "ref")) %>%
    filter(alt.x == alt.y) %>% # filter out triallelic sites
    rename(alt=alt.x) %>%
    select(-alt.y) %>%
    write_tsv(output_file)
