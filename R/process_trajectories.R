library(magrittr)
library(dplyr)
library(reshape2)
library(ggplot2)


##
## Create a data frame with Nea. trajectories from all replicates.
##
load_trajectories <- function(sim_dir, scenario, h, init_nea, sites) {
    # compose the pattern describing all input tables
    pattern <- paste0(scenario,
                      "__h_", h,
                      "__init_nea_", init_nea,
                      "__rep_.*",
                      "__", sites, "_sites.txt")

    # get paths of all input tables
    files <- list.files(sim_dir, pattern, full.names=TRUE)

    # load each table individually
    tables <- lapply(seq_along(files),
           function(i) {
                   read.table(files[i], header=TRUE) %>%
                       mutate(rep=i,
                              model=scenario,
                              sites=sites)
           })

    # merge all dataframes into one
    do.call(rbind, tables)
}


##
## Plot one statistic of Nea. trajectories over time
## separately for each replicate simulation.
##
plot_one_stat <- function(df, stat="mean", title="") {
    if (! stat %in% c("mean", "sd", "median", "min", "max")) {
        error(paste("Invalid admixture statistic provided:", stat, "!"))
    }
    
    df <- melt(df, id=c("gen", "rep"), measure=c(stat))
    p1 <-
        ggplot(df, aes(x=gen, y=value, group=rep, color=rep)) +
        geom_line(alpha=0.6) +
        labs(title="Average stats of Nea. ancestry over time",
             x="generation since the admixture",
             y="proportion of Nea. ancestry") +
        ylim(0, 0.15) +
        theme(legend.position="none")
    
    p2 <- p1 + scale_x_log10()
    gridExtra::grid.arrange(p1, p2, ncol=2)
}


##
## Plot means of all statistics of Nea. trajectories over time
## from all replicates.
##
plot_all_stats <- function(df, log_scale=FALSE, title="")
{
    df %<>%
        select(-model, -sites) %>%
        group_by(gen) %>%
        summarize_each(funs(mean)) %>%
        melt(id=c("gen"), measure=c("mean", "sd", "min", "max"))

    p <-
        ggplot(df, aes(x=gen, y=value, color=variable)) +
        geom_line() +
        labs(title=paste("Average stats of Nea. ancestry over time", title),
         x="generation since the admixture",
         y="proportion of Nea. ancestry") + ylim(0, 0.15)

    if (log_scale) {
        p <- p + scale_x_log10()
    }
 
    p
}
