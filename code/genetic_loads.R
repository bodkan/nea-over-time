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
  func_haps$nea_fit <- sapply(split(nea_muts[mcols(nea_muts)[[chrom_id]] == 1][queryHits(hits)]$S, subjectHits(hits)),
                                  function(muts) exp(sum(muts)))
  func_haps
}

#' Given a Neanderthal haplotype, find a set of overlapping modern human
#' haplotypes and calculate their genetic loads.
matching_mh_hap <- function(nea_hap, markers, mh_muts) {
  for (chrom_id in sample(chrom_ids(markers))) {
    # get the sites of fixed markers overlapping this N haplotype
    marker_sites <- subsetByOverlaps(markers[, chrom_id], nea_hap)
    # this MH chromosome doesn't overlap any fixed N markers - it's pure MH
    if (sum(mcols(marker_sites)[[chrom_id]]) == 0) {
      # => find deleterious MH mutations of the MH haplotype overlapping a N hap
      mh_del_muts <- mh_muts[mcols(mh_muts)[[chrom_id]] == 1, c("S", chrom_id)]
      overlap_del_muts <- subsetByOverlaps(mh_del_muts, nea_hap)
      # calculate the fitness of that haplotype caused by these mutations
      nea_hap$mh_fit <- exp(sum(overlap_del_muts$S))
      return(nea_hap)
    }
  }
  stop("This should never happen!")
}



nea_mh_fitness <- function(vcf_path) {
  # vcf_path <- "data/simulations/exon_h_0.5_rep_1_gen_1.vcf.gz"
  vcf <- readVcf(vcf_path)

  # load the GTs of neutral markers on each simulated chromosome
  markers <- mut_genotypes(vcf=vcf, mut_type=1)
  # load the GTs of new non-African mutations
  new_muts <- mut_genotypes(vcf, mut_type=0, pop_origin=3)
  # load the GTs of MH and N deleterious mutations and merge them with
  # the new non-African mutations
  mh_muts <- c(mut_genotypes(vcf=vcf, mut_type=0, pop_origin=1, t_min=70000), new_muts) %>% sort
  nea_muts <- c(mut_genotypes(vcf=vcf, mut_type=0, pop_origin=2), new_muts) %>% sort

  # find N haplotypes on each simulated chromosome and calculate their fitness
  nea_haps <- chrom_ids(markers) %>%
    lapply(find_nea_haps, markers, nea_muts) %>%
    compact %>%
    Reduce(c, .)

  # find "competing" MH haplotypes for each identified N haplotype
  nea_vs_mh_hap <-
      lapply(nea_haps, matching_mh_hap, markers, mh_muts) %>%
      Reduce(c, .) %>%
      as.data.frame %>%
      select(nea_fit, mh_fit) %>%
      as_tibble

  nea_vs_mh_hap
}

#
# ggplot(fits, aes(nea_fit, mh_mean)) +
#     geom_point() +
#     geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci)) +
#     geom_abline(slope=1) +
#     xlim(0, 1) + ylim(0, 1)

# mutation frequencies before introgression
# ggplot() +
#   geom_point(aes(log10(-mh_muts$S), mh_muts$freq / max(mh_muts$freq), color="blue"), alpha=1/10) +
#   geom_point(aes(log10(-nea_muts$S), nea_muts$freq / max(nea_muts$freq), color="red"), alpha=1/10) +
#   xlim(-15, 0)

