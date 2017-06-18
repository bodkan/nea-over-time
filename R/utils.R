#
# Replace het calls in a set of samples with randomly called REF or
# ALT hom calls.
#
# 'df' is simply a dataframe with columns: chrom, pos, ref, alt and
# then one column per individual, containing values:
#    - 0 - no ALT alleles observed (hom reference)
#    - 1 - 1 ALT allele observed (het)
#    - 2 - 2 ALT alleles observed (hom alternative)
#    - 9 - missing data
#
random_call <- function(df, sample_names=NULL) {
    # if no samples were specified, perform random calls on all
    # samples
    if (is.null(sample_names)) {
        sample_names <- colnames(df)[5 : ncol(df)]
    }
    
    for (s in sample_names) {
        # get the positions of het sites in a given individual
        het_pos <- df[[s]] == 1

        # randomly sample hom REF/ALT alleles at het positions
        df[[s]] <- replace(df[[s]],
                           het_pos,
                           2 * rbinom(sum(het_pos, na.rm=TRUE), 1, 0.5)) %>% as.integer
    }

    df
}

#
# Calculate the proportion of alleles shared between sample_a and
# sample_b.
#
calc_sharing_prop <- function(snps, sample_a, sample_b) {
    # take positions of SNPs where both samples have a valid allele
    # (9 encodes a missing position)
    available_positions <- !is.na(snps[[sample_a]]) & !is.na(snps[[sample_b]])
    
    mean(snps[available_positions, sample_a] == snps[available_positions, sample_b])
}

#
# Remove SNPs that could be enriched for aDNA errors.
#
remove_transitions <- function(snps) {
    filter(snps, !(
        (ref == "C" & alt == "T") |
        (ref == "T" & alt == "C") |
        (ref == "G" & alt == "A") |
        (ref == "A" & alt == "G")
    ))
}

#
# Estimate Neanderthal ancestry in a given set of samples.
#
estimate_nea <- function(snps, sample_names) {
    nea_estimates <- sapply(sample_names,
                            function(s) {
                                calc_sharing_prop(snps, "archaic_Altai", s)
                            })
    nea_estimates
}

#
# Load the SGDP metadata info table.
#
load_sgdp_info <- function(path) {
    read_tsv(path) %>%
        select(Panel, name=SGDP_ID, Region, Country, Latitude, Longitude) %>%
        filter(complete.cases(.)) %>%
        filter(Panel == "C") %>%
        mutate(name=str_replace(name, "-", "_")) %>%
        select(-Panel)
}

load_annotations <- function(annotations_path) {
    read_tsv(annotations_path, progress=FALSE) %>%
        rename(chrom=Chrom, pos=Pos, ref=Ref, alt=Alt)
}

#
# Load the whole data set of EMH, SGDP and archaic human SNPs
# and combine it with the CADD annotations.
#
load_dataset <- function(ice_age_path,
                         sgdp_path,
                         archaics_path,
                         filter_damage,
                         metadata_path,
                         random_sample=T) {
    ## ice_age_path <- "clean_data/ice_age.tsv"
    ## sgdp_path <- "clean_data/sgdp.tsv"
    ## archaics_path <- "clean_data/archaics.tsv"
    ## annotations_path <- "clean_data/annotations.tsv"
    
    # load the genotypes from all individuals and call random alleles for
    # humans with diploid calls
    ice_age <-
        read_tsv(ice_age_path, progress=FALSE)
    if (random_sample) ice_age %<>% random_call(c("Loschbour", "Stuttgart", "UstIshim"))

    if (filter_damage) {
        ice_age <- remove_transitions(ice_age)
    }

    # read the list of samples with available metadata
    sgdp_info <- load_sgdp_info(metadata_path)

    sgdp <-
        read_tsv(sgdp_path, progress=FALSE) %>%
        select(c(chrom, pos, ref, alt,
                 one_of(sgdp_info$name)))
    if (random_sample) sgdp %<>% random_call
    names(sgdp)[-(1 : 4)] %<>% str_replace("^S_", "")

    archaics <-
        read_tsv(archaics_path, progress=FALSE) %>%
        filter(Altai == 2, Vindija == 2)
    names(archaics)[-(1 : 4)] %<>% str_replace("^", "archaic_")    

    all_samples <-
        right_join(ice_age, archaics) %>% # join archaic data with EMH
        { replace(., is.na(.), 9L) } %>%  # replace missing EMH data with 9 "alleles"
        left_join(sgdp) %>%               # add variable SGDP sites
        { replace(., is.na(.), 0L) }      # fill in missing SGDP alleles as REF

    # replace all 9 values with NA
    pos_bases <- all_samples[, 1:4] # ignore the first 4 columns
    snps <- all_samples[, -(1:4)]   # this is the part to modify
    snps[snps == 9] <- -1   # first replace all 9s in the SNP matrix with -1
    snps[snps == -1] <- NA  # then substitute for NA
    
    # ice_age %>% {sapply(colnames(.)[5:ncol(.)], function(s) {table(.[[s]])})}
    # sgdp %>% {sapply(colnames(.)[5:ncol(.)], function(s) {table(.[[s]])})}
    # archaics %>% {sapply(colnames(.)[5:ncol(.)], function(s) {table(.[[s]])})}

    bind_cols(pos_bases, snps)
}
