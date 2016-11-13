library(magrittr)
library(dplyr)
library(stringr)
library(ggplot2)


## Load information about samples.
load_sample_info <- function(path)
{
    exclude_samples <- c("Altai", "Chimp", "Denisova", "B_Dinka_3",
                         "aurig", "Href", "Kostenki14_SG",
                         "B_Mandenka_3", "B_Mbuti_4", "MezE",
                         "I0908_published", "Oase1_pu_own",
                         "Vi_merge", "B_Yoruba_3", "B_Han_3",
                         "B_Papuan_15", "B_Australian_4",
                         "B_Karitiana_3", "B_Dai_4", "Oase1_d")
    
    read.delim(path,
               stringsAsFactors=FALSE,
               col.names=c("x", "sample_id", "sample_name",
                           "date", "x", "x")) %>%
        .[order(.$date, decreasing=TRUE), ] %>%                    # order by age
        mutate(sample_id=str_replace_all(sample_id, "-", "_")) %>% # replace `-`s with `_`s
        mutate(time_after_admix=55000 - date) %>%                  # time since introgression
        mutate(estimate_nea=! sample_id %in% exclude_samples) %>% 
        select(-matches("^x")) %>%                                 # ignore useless columns
        return()
}


## Calculate the proportion of sharing with between given
## two individuals
calc_proportion <- function(snps, sample_a, sample_b) {
    # take positions of SNPs where both samples have a valid allele
    # (9 encodes a missing position)
    available_positions <- (snps[, sample_a] != 9) & (snps[, sample_b] != 9)
    
    return(mean(snps[available_positions, sample_a] == snps[available_positions, sample_b]))
}


## Split the given list of SNPs into a defined number of blocks.
split_into_blocks <- function(snps, n_blocks) {
    block_size <- nrow(snps) %/% n_blocks
    block_breaks <- seq(1, nrow(snps), block_size)

    # split the list of SNPs into blocks
    snps_blocks <- split(snps, findInterval(seq(1, nrow(snps)), block_breaks))[1 : n_blocks]
       
    return(snps_blocks)
}


## Calculate Neanderthal ancestry estimates in leave-one-block-out
## subsets of the data.
calc_nea_partials <- function(snps, n_blocks, sample_ids) {
    # initialize the matrix of blocks (rows) for each sample (columns)
    results <- matrix(, nr=n_blocks, nc=length(sample_ids))
    colnames(results) <- sample_ids
    
    # split the given set of SNPs into blocks and calculate proportions on 
    # leave-one-out concatenated blocks
    snp_blocks <- split_into_blocks(snps, n_blocks)

    for (i in seq(n_blocks)) {
        # concatenate blocks without the i-th block
        leave_one_block_out_snps <- do.call(rbind, snp_blocks[-i])
        
        # calculate the proportion of sharing with Altai in Dinka (assumed 0% Nea.)
        # and Vindija (assumed 100% Nea.)
        E_dinka <- calc_proportion(leave_one_block_out_snps, "Altai", "B_Dinka_3")
        E_vindija <- calc_proportion(leave_one_block_out_snps, "Altai", "Vi_merge")
        
        for (sample_id in sample_ids) {
            E_sample <- calc_proportion(leave_one_block_out_snps, "Altai", sample_id)
            
            # convert the raw proportions of matching Altai into Nea. ancestry estimates
            results[i, sample_id] <- (E_sample - E_dinka) / (E_vindija - E_dinka)
        }
    }

    return(results)
}


## Calculate whole-genome Nea% estimate using rates of matching
## of Altai with Dinka and Vindija.
calc_mean_nea <- function(snps, sample_ids) {
    E_dinka <- calc_proportion(snps, "Altai", "B_Dinka_3")
    E_vindija <- calc_proportion(snps, "Altai", "Vi_merge")

    mean_nea <- sapply(include_samples, function(s) {
        (calc_proportion(snps, "Altai", s) - E_dinka) / (E_vindija - E_dinka)
    })

    return(mean_nea)
}


# Calculate genome-wide statistics of Nea% using jackknife on
# the precalculated values of the leave-one-block-out subsamples.
calc_jackknife_nea <- function(nea_in_partials, mean_nea) {
    data.frame(
        sample_id=colnames(nea_in_partials),

        # mean Nea. ancestry over all SNPs
        mean_nea=mean_nea,

        # mean of Nea. ancestries across all leave-one-out subsamples
        mean_partials=apply(nea_in_partials, 2, mean),
        
        stringsAsFactors=FALSE
    ) %>%
    mutate(
        # bias of the estimate and the bias-corrected jackknife estimate
        jack_bias=(n_blocks - 1) * (mean_partials - mean_nea),
        jack_mean=mean_nea - jack_bias,

        # standard error of the Nea. ancestry
        jack_se=apply(nea_in_partials, 2, function(sample_partials) {
             sqrt((n_blocks - 1) / n_blocks * sum((sample_partials - mean(sample_partials))^2))
        }),

        # confidence intervals
        jack_ci_low=jack_mean - qt(0.975, n_blocks - 1) * jack_se,
        jack_ci_up=jack_mean + qt(0.975, n_blocks - 1) * jack_se,
        
        # confidence intervals
        paper_ci_low=mean_partials - qt(0.975, n_blocks - 1) * jack_se,
        paper_ci_up=mean_partials + qt(0.975, n_blocks - 1) * jack_se
    ) %>%
    return()
}
