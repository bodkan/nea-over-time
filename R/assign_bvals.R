library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(stringr)

assign_bvals <- function(sites, bval_path, chain_path) {
    # load the dataframes with McVicker scores and convert into a BED-like format
    bval_files <- list.files(bval_path, full.names=TRUE, ".*.bkgd")
    bval_df_list <- lapply(bval_files, function(filename) {
        read.table(filename, col.names=c("bval", "length")) %>%
            mutate(chr=str_replace(basename(filename), ".bkgd", ""),
                   end=cumsum(length),
                   start=c(1, (end + 1)[-n()])) %>%
            select(chr, start, end, bval)
    })

    # convert the list of dataframes into a GRanges hg18 object
    bval_regions_hg18 <- bind_rows(bval_df_list) %>%
        makeGRangesFromDataFrame(keep.extra.columns=TRUE)
    seqinfo(bval_regions_hg18) <- Seqinfo(genome="hg18")

    # perform the liftOver from hg18 to hg19
    chain <- import.chain(chain_path)
    bval_regions <- liftOver(bval_regions_hg18, chain) %>% unlist

    # take only sites which fall within some B value annotated region
    # (because not the whole genome is annotated)
    bval_sites <- subsetByOverlaps(sites, bval_regions)

    # get the pairings of SNP sites and corresponding B value regions
    hits <- findOverlaps(sites, bval_regions)

    # take the B values from the regions overlapping sites
    bvals <- IntegerList(split(bval_regions$bval[subjectHits(hits)], queryHits(hits)))

    # assign B values to each SNP
    mcols(bval_sites)[["bval"]] <- bvals

    # there are some weird cases where a SNP falls within two different regions,
    # which could be a bug in the liftOver function - remove these values
    mcols(bval_sites[which(elementNROWS(bval_sites$bval) > 1)]) <- NA

    # convert an IntegerList to a normal integer vector
    bval_sites$bval <- unlist(bval_sites$bval)

    bval_sites
}
