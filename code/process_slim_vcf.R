library(VariantAnnotation)
library(rtracklayer)
library(tidyverse)

#' Load mutation info about a given mutation type.
mut_info <- function(vcf, mut_type, pop_origin=NULL, t_min=-Inf, t_max=Inf) {
  mut_pos <- info(vcf)$MT == mut_type &
    info(vcf)$GO >= t_min &
    info(vcf)$GO <= t_max

  if (!is.null(pop_origin)) {
    mut_pos <- mut_pos & (info(vcf)$PO == pop_origin)
  }

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
mut_gt <- function(vcf, mut_type, pop_origin=NULL, t_min=-Inf, t_max=Inf) {
  mut_pos <- info(vcf)$MT == mut_type &
    info(vcf)$GO >= t_min &
    info(vcf)$GO <= t_max

  if (!is.null(pop_origin)) {
    mut_pos <- mut_pos & (info(vcf)$PO == pop_origin)
  }

  gr <- granges(vcf)[mut_pos]
  
  gt_mat <- geno(vcf)$GT[mut_pos, ]
  
  ind_gts <- apply(gt_mat, 2, function(gt) { str_count(gt, "1") })
  
  info_cols <- as.data.frame(info(vcf)[mut_pos, c("S", "DOM", "PO", "GO")])
  
  mcols(gr) <- bind_cols(info_cols, as.data.frame(ind_gts))
  names(gr) <- NULL
  
  # shift VCF coordinates back to the SLiM 0-based system
  gr <- shift(gr, shift=-1)
  
  sort(gr)
}

#' Read BED file with coordinates of a given type of genomic elements.
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
#' The argument names are kind of confusing :(
#'    - real_coords = table of coordinates of both region & gap sites
#'                    (in realistic coordinate system & SLiM system)
#'    - gap_coords = table of coordinates of gap sites in realistic system
#' Both are needed because the real_coord files don't specify which are gaps
#' and which are not :/
get_markers <- function(vcf, real_coords, gap_coords) {
  # coordinates of all sites (region and gap sites) - real & SLiM coords
  real_sites <- read_sites(real_sites)
  # gap sites
  gap_sites <- import.bed(gap_sites) %>% as.data.frame %>% select(chrom=seqnames, pos=end)
  # simulated data
  sim_sites <- mut_info(vcf, mut_typ=1) %>% shift(shift=-1)
  
  trans_sites <- transpose_sites(sim_sites, real_sites) %>%
    as.data.frame %>% select(chrom=seqnames, pos=start, freq) %>%
    mutate(chrom=as.character(chrom))

  all_sites <- mcols(real_sites) %>%
    as.data.frame %>% select(chrom=real_chrom, pos=real_end)

  full_join(trans_sites, all_sites, by=c("chrom", "pos")) %>%
    mutate(freq=ifelse(is.na(freq), 0, freq),
           chrom=factor(chrom, levels=paste0("chr", 1:22))) %>%
    arrange(chrom, pos) %>%
    inner_join(gap_sites, by=c("chrom", "pos"))
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

#' Extract coordinates of deserts from a given VCF file.
get_deserts <- function(markers) {
    all_chrom <- list()
    for (chrom in paste0("chr", 1:22)) {
        sites <- markers[markers$chrom == chrom, ]
        desert_runs <- rle(as.integer(sites$freq > 0))

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
