library(VariantAnnotation)
library(tidyverse)

#' Read a VCF file simulated by SLiM. Return as a GRanges object.
read_vcf <- function(path) {
  readVcf(path, param=ScanVcfParam(info=c("S", "DOM", "PO", "GO", "MT"), geno="GT"))
}

#' Load mutations of a given type from SLiM.
#' Mutation type 0 are deleterious mutations, mutation type 1 are neutral
#' markers.
get_muts <- function(vcf, mut_type) {
  mut_pos <- info(vcf)$MT == mut_type
  gr <- granges(vcf)[mut_pos]
  
  gt_mat <- geno(vcf)$GT[mut_pos, ]
  
  ind_gts <- apply(gt_mat, 2, function(gt) { str_count(gt, "1") })
  freq <- apply(ind_gts, 1, sum) / (2 * ncol(ind_gts))
  
  info_cols <- as.data.frame(info(vcf)[mut_pos, c("S", "DOM", "PO", "GO")])
  
  mcols(gr) <- bind_cols(info_cols, as.data.frame(freq), as.data.frame(ind_gts))
  names(gr) <- NULL
  
  # shift VCF coordinates back to the SLiM 0-based system
  gr <- shift(gr, shift=-1)
  
  sort(gr)
}

#' Read BED file with coordinates of a given type of genomic elements.
#' Allowed regions are: exon, promoter, protein_coding, tf_binding_site and utr3.
read_regions <- function(regions_bed) {
  gr <- read_tsv(regions_bed, col_types="ciicdiii") %>%
    makeGRangesFromDataFrame(starts.in.df.are.0based=TRUE)
  gr
}

#' Read the real coordinates of Neanderthal fixed sites.
read_sites <- function(sites_bed) {
  gr <- read_tsv(sites_bed, col_types="ciiii") %>%
    select(real_chrom=chrom, real_start=start, real_end=end, start=slim_start, end=slim_end) %>%
    mutate(chrom=1) %>% 
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)
  gr
}

#' Transpose the simulated sites in the SLiM 0-based coordinate system into
#' the realistic coordinates.
transpose_sites <- function(sim_sites, real_sites) {
  hits <- findOverlaps(real_sites, sim_sites)
  transposed <- as.data.frame(mcols(real_sites)) %>%
    setNames(c("chrom", "start", "end")) %>%
    makeGRangesFromDataFrame(starts.in.df.are.0based=TRUE) %>%
    .[queryHits(hits)]
  mcols(transposed) <- mcols(sim_sites)["freq"]
  
  transposed
}

#' Load the simulated Neanderthal fixed markers and their frequencies.
#' Return as a GRanges object.
get_markers <- function(vcf, regions_bed) {
  real_sites <- read_sites(regions_bed)
  sim_sites <- get_muts(vcf, mut_type=1)
  
  trans_sites <- transpose_sites(sim_sites, real_sites) %>%
    as.data.frame %>% select(chrom=seqnames, pos=start, freq) %>%
    mutate(chrom=as.character(chrom))

  all_sites <- mcols(real_sites) %>%
    as.data.frame %>% select(chrom=real_chrom, pos=real_end)

  full_join(trans_sites, all_sites, by=c("chrom", "pos")) %>%
    mutate(freq=ifelse(is.na(freq), 0, freq),
           chrom=factor(chrom, levels=paste0("chr", 1:22))) %>%
    arrange(chrom, pos) %>%
    inner_join(as.data.frame(import.bed("/mnt/scratch/mp/nea-over-time/data/bed/regions/gap_sites.bed")) %>% select(chrom=seqnames, pos=end))
}

#' Calculate the number of Nea. mutations (given in a GRanges object) in all
#' simulated individuals. Return result as an integer vector.
#' 
#' Nea. ancestry proportion can be calculated as: result / (2 * mut_count)
nea_per_ind <- function(gr) {
  ind_counts <- as.data.frame(mcols(gr)) %>%
    select(starts_with("i")) %>%
    summarise_all(sum) %>%
    t %>%
    as.vector
  
  ind_counts
}
