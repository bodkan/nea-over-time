library(VariantAnnotation)
library(rtracklayer)
library(tidyverse)

#' Split the GT matrix into haplotype GT matrix.
split_chromosomes <- function(gt_mat) {
    hap_mat <- matrix(nrow=nrow(gt_mat), ncol=ncol(gt_mat) * 2)
    colnames(hap_mat) <- paste0("chr", 1:ncol(hap_mat))
    for (i in seq_len(ncol(gt_mat))) {
        k <- (2 * i) - 1
        hap_mat[, c(k, k + 1)] <- str_split_fixed(gt_mat[, i], "\\|", 2) %>% as.integer
    }
    hap_mat
}

#' Load mutation info about a given mutation type.
mut_info <- function(vcf, mut_type=NULL, pop_origin=NULL, t_min=-Inf, t_max=Inf) {
  mut_pos <- filter_muts(vcf, mut_type, pop_origin, t_min, t_max)

  gr <- granges(vcf)[mut_pos]
  names(gr) <- NULL

  mcols(gr) <- as.data.frame(info(vcf)) %>%
    filter(mut_pos) %>%
    mutate(freq = AC / (2 * length(samples(header(vcf)))))

  gr
}

#' Load mutations of a given type from SLiM.
#' Mutation type 0 are deleterious mutations, mutation type 1 are neutral
#' markers.
mut_genotypes <- function(vcf, mut_type=NULL, pop_origin=NULL, t_min=-Inf, t_max=Inf) {
  mut_pos <- filter_muts(vcf, mut_type, pop_origin, t_min, t_max)

  gr <- granges(vcf)[mut_pos]
  
  gt_mat <- geno(vcf)$GT[mut_pos, ]
  hap_mat <- split_chromosomes(gt_mat)
  info_cols <- as.data.frame(mcols(mut_info(vcf, mut_type, pop_origin, t_min, t_max)))
  
  mcols(gr) <- bind_cols(info_cols, as.data.frame(hap_mat))
  names(gr) <- NULL
  
  # shift VCF coordinates back to the SLiM 0-based system
  gr <- shift(gr, shift=-1)
  
  sort(gr)
}

#' Load the simulated Neanderthal fixed markers and their frequencies
#' (either all, or the ones within regions or gaps). Return as a GRanges
#' object.
get_markers <- function(vcf, sites_coords, within_region=NULL) {
  if (!(is.null(within_region) || within_region %in% c("region", "gap")))
      stop("Invalid region specified")

  # coordinates of all sites (real & SLiM coords)
  real_sites <- read_sites(sites_coords)
  # simulated data
  sim_sites <- mut_info(vcf, mut_type=1) %>% shift(shift=-1)
  
  trans_sites <- transpose_sites(sim_sites, real_sites) %>%
    as.data.frame %>% select(chrom=seqnames, pos=start, freq) %>%
    mutate(chrom=as.character(chrom))

  all_sites <- mcols(real_sites) %>%
    as.data.frame %>% select(chrom=real_chrom, pos=real_end, within)

  markers <- full_join(trans_sites, all_sites, by=c("chrom", "pos")) %>%
    mutate(freq=ifelse(is.na(freq), 0, freq),
           chrom=factor(chrom, levels=paste0("chr", 1:22))) %>%
    arrange(chrom, pos)

  if (!is.null(within_region))
    markers <- markers[markers$within == within_region, ]

  markers
}

#' Calculate the number of Nea. mutations (given in a GRanges object) in all
#' simulated individuals. Return result as an integer vector.
#' 
#' Nea. ancestry proportion can be calculated as: result / (2 * mut_count)
nea_per_ind <- function(gr) {
  ind_counts <- as.data.frame(mcols(gr)) %>%
    select(-c(S, DOM, PO, GO)) %>%
    summarise_all(sum) %>%
    t %>%
    as.vector
  
  ind_counts
}

#' Extract coordinates of deserts from a given VCF file.
get_deserts <- function(markers, cutoff=0) {
    all_chrom <- list()
    for (chrom in paste0("chr", 1:22)) {
        sites <- markers[markers$chrom == chrom, ]
        desert_runs <- rle(as.integer(sites$freq > cutoff))
        if (length(desert_runs$values) < 2) {
            all_chrom[[chrom]] <- NULL
            next
        }

        block_idx <- c(0, desert_runs$lengths %>% cumsum)
        block_start <- block_idx[-length(block_idx)]
        block_end <- block_idx[2:length(block_idx)]

        desert_start <- block_start[desert_runs$values == 0]
        desert_end <- block_end[desert_runs$values == 0]

        all_chrom[[chrom]] <- GRanges(chrom, IRanges(start=sites[desert_start + 1, ]$pos,
                                                     end=sites[desert_end, ]$pos))
    }
    Reduce(c, GRangesList(all_chrom))
}

#' Download the coordinates of centromeres.
get_centromeres <- function() {
    session <- browserSession("UCSC")
    genome(session) <- "hg19"
    query <- ucscTableQuery(session, "gap", GRangesForUCSCGenome("hg19", chrom=paste0("chr", 1:22)))

    tbl <- getTable(query)

    gaps <- filter(tbl, type %in% c("centromere"), chrom %in% paste0("chr", 1:22)) %>%
        makeGRangesFromDataFrame(keep.extra.columns=TRUE)

    gaps
}

#' Read BED file with coordinates of a given type of genomic elements.
read_regions <- function(regions_bed, slim=FALSE) {
  df <- read_tsv(regions_bed, col_types="ciicdiii")
  if (slim) {
    df %>% mutate(chrom="1") %>% select(chrom, slim_start, slim_end) %>%
      makeGRangesFromDataFrame(starts.in.df.are.0based=TRUE)
  } else {
      makeGRangesFromDataFrame(df, starts.in.df.are.0based=TRUE)
  }
}

#' Read the real coordinates of Neanderthal fixed sites.
read_sites <- function(sites_bed) {
  gr <- read_tsv(sites_bed, col_types="ciiiic") %>%
    select(real_chrom=chrom, real_start=start, real_end=end, start=slim_start, end=slim_end, within) %>%
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
  mcols(transposed) <- mcols(sim_sites)
  
  transposed
}

#' Filter mutations' positions based on given criteria.
filter_muts <- function(vcf, mut_type=NULL, pop_origin=NULL, t_min=-Inf, t_max=Inf) {
  mut_pos <- info(vcf)$GO >= t_min & info(vcf)$GO <= t_max

  if (!is.null(mut_type)) {
    mut_pos <- mut_pos & info(vcf)$MT == mut_type
  }

  if (!is.null(pop_origin)) {
    mut_pos <- mut_pos & (info(vcf)$PO == pop_origin)
  }

  mut_pos
}

#' Get names of haplotypes in a given GT GRanges object.
chrom_ids <- function(gr) {
  keep(colnames(mcols(gr)), str_detect, "^chr")
}
