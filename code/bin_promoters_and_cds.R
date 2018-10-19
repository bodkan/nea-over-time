suppressPackageStartupMessages({
  library(rtracklayer)
  library(tidyverse)
})

gr_to_bed <- function(gr, path) {
  gr %>% as.data.frame %>%
    rename(chrom = seqnames) %>%
    select(-width, strand) %>%
    mutate(start = as.integer(start - 1), end = as.integer(end)) %>%
    write_tsv(path)
}

promoters <- suppressMessages(read_tsv("data/slim_coords/promoter_regions.bed") %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE, starts.in.df.are.0based=TRUE))
merged <- suppressMessages(read_tsv("data/slim_coords/merged_regions.bed") %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE, starts.in.df.are.0based=TRUE))

hits <- findOverlaps(merged, promoters)

h_scale <- c("0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0")
s_scale <- c("0.0", "0.25", "0.5", "0.75", "1.0", "1.25", "1.5", "1.75", "2.0")

for (h in h_scale) {
  # generate regions for dominance differences between promoters and CDS
  h_regions <- merged
  h_regions$bin_h <- 0.5
  h_regions[queryHits(hits)]$bin_h <- h
  gr_to_bed(h_regions, file.path("data/slim_coords/", glue::glue("different_h_{h}_merged_regions.bed")))
}

for (s in s_scale) {
  # generate regions for selective differences between promoters and CDS
  s_regions <- merged
  s_regions$bin_s <- 1.0
  s_regions[queryHits(hits)]$bin_s <- s
  gr_to_bed(s_regions, file.path("data/slim_coords/", glue::glue("different_s_{s}_merged_regions.bed")))
}
