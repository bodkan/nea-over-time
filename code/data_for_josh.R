source("code/process_slim_vcf.R")

out_dir <- "tmp/josh_sfs"
dir.create(out_dir)

regions <- c("exon", "promoter", "protein_coding", "tf_binding_site", "utr3", "neutral_Ne_10000")

for (rep_i in 1:3) {
  for (region in regions) {
    vcf <- readVcf(paste0("data/simulations/deserts_", region, "_h_0.5_rep_", rep_i, "_gen_2200.vcf.gz"))
    r <- ifelse(str_detect(region, "neutral"), "tf_binding_site", region)

    in_markers <- get_markers(vcf, paste0("data/slim_coords/", r, "_all_sites.bed"), within_region="region")
    out_markers <- get_markers(vcf, paste0("data/slim_coords/", r, "_all_sites.bed"), within_region="gap")

    markers <- bind_rows(in_markers, out_markers) %>%
      arrange(chrom, pos) %>%
      select(chrom, pos, freq, within) %>%
      write_tsv(file.path(out_dir, paste0(ifelse(str_detect(region, "neutral"), "neutral", region), "_rep", rep_i, ".tsv")))
  }
}

system("cd tmp; tar -zcvf josh_sfs.tar.gz josh_sfs/")