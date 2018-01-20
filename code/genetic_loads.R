#' Find all Neanderthal haplotypes on a given chromosome using the fixed
#' neutral Neanderthal markers.
find_nea_haps <- function(chrom_id, markers, nea_muts) {
  nea_runs <- rle(as.integer(mcols(markers)[[chrom_id]] == 1))
  
  # if there are no N markers, then this haplotype is 100% human
  if (!any(nea_runs$values)) return(NULL)
  
  block_idx <- c(0, nea_runs$lengths %>% cumsum)
  block_start <- block_idx[-length(block_idx)]
  block_end <- block_idx[2:length(block_idx)]
  
  # if there is only one run of N markers, this haplotype is 100% Neanderthal
  nea_pos <- if (length(nea_runs$values) > 1) nea_runs$values == 1 else 1
  hap_start <- block_start[nea_pos]
  hap_end <- block_end[nea_pos]

  haps <- GRanges(unique(seqnames(markers)),
                  IRanges(start=start(markers[hap_start + 1, ]),
                          end=start(markers[hap_end, ])))
  
  hits <- findOverlaps(nea_muts[mcols(nea_muts)[[chrom_id]] == 1], haps)
  func_haps <- haps[unique(subjectHits(hits))]
  func_haps$S <- split(nea_muts[mcols(nea_muts)[[chrom_id]] == 1][queryHits(hits)]$S, subjectHits(hits))
  func_haps
}

#' Given a Neanderthal haplotype, find a set of overlapping modern human
#' haplotypes and calculate their genetic loads.
find_mh_haps <- function(chrom_id, nea_haps, markers, mh_muts) {
  # find already detected N haplotypes overlapping this chromosome
  chrom_markers <- markers[, chrom_id]
  hits <- findOverlaps(chrom_markers, nea_haps)

  nea_states_per_hap <- !sapply(split(mcols(chrom_markers[queryHits(hits)])[[chrom_id]], subjectHits(hits)), all)
  if (!any(nea_states_per_hap)) return(NULL)
  nea_haps_overlap <- nea_haps[nea_states_per_hap]

  # find deleterious mutations falling within pure MH haplotypes
  mh_muts_hits <- findOverlaps(mh_muts[mcols(mh_muts)[[chrom_id]] == 1, chrom_id], nea_haps_overlap)
  nea_haps_overlap$mh_genload <- exp(sum(split(mh_muts[mcols(mh_muts)[[chrom_id]] == 1][queryHits(mh_muts_hits)]$S, subjectHits(mh_muts_hits))[[1]]))
  nea_haps_overlap$mh_chrom_id <- chrom_id

  as_tibble(as.data.frame(nea_haps_overlap)[c("nea_hap_id", "nea_genload", "mh_genload", "mh_chrom_id")])
}


loads_in_gen <- function(gen, vcf_path) {
  vcf <- readVcf(vcf_path)
  
  # load the GTs of neutral markers on each simulated chromosome
  markers <- mut_genotypes(vcf=vcf, mut_type=1)
  # load the GTs of MH and N deleterious mutations
  mh_muts <- mut_genotypes(vcf=vcf, mut_type=0, pop_origin=1, t_min=70000)
  nea_muts <- mut_genotypes(vcf=vcf, mut_type=0, pop_origin=2)
  
  # find N haplotypes on each simulated chromosome
  nea_haps <- chrom_ids(markers) %>%
    lapply(function(chrom_id) find_nea_haps(chrom_id, markers, nea_muts)) %>%
    setNames(chrom_ids(markers)) %>% compact
  # calculate the genetic load of each N haplotype in the population
  nea_hap_loads <- Reduce(c, nea_haps) %>% as.data.frame %>% as_tibble %>%
    mutate(nea_genload=map_dbl(S, ~ exp(sum(.x)))) %>%
    arrange(seqnames, start, end) %>% select(-S) %>% distinct %>%
    rownames_to_column("nea_hap_id") %>% mutate(nea_hap_id=as.integer(nea_hap_id)) %>% 
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)
  
  chrom_ids(markers) %>%
    lapply(function(chrom_id) find_mh_haps(chrom_id, nea_hap_loads, markers, mh_muts)) %>%
    setNames(chrom_ids(markers)) %>% bind_rows
}



# mutation frequencies before introgression
# ggplot() +
#   geom_point(aes(log10(-mh_muts$S), mh_muts$freq / max(mh_muts$freq), color="blue"), alpha=1/10) +
#   geom_point(aes(log10(-nea_muts$S), nea_muts$freq / max(nea_muts$freq), color="red"), alpha=1/10) +
#   xlim(-15, 0)
