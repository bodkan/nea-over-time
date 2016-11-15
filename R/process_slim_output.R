library(magrittr)
library(dplyr)
library(stringr)
library(ggplot2)

MUTATION_COLS <- c("file_id", "mut_id", "mut_type", "pos", "s", "h", "pop_origin",
                    "gen_origin", "freq")

get_section_delims <- function(slim_file) {
    section_delim <- c(which(!is.na(str_extract(slim_file, "^[a-zA-X]+:$"))),
                       length(slim_file))

    data.frame(section=c("populations", "mutations", "individuals", "genomes"),
                         start=section_delim[-length(section_delim)] + 1,
                         end=section_delim[-1] - 1)
}
    
read_section_data <- function(section_name, delims, slim_file) {
    pos_in_file <- filter(delims, section == section_name)
    return(slim_file[pos_in_file$start:pos_in_file$end])
}

read_mutations <- function(slim_file, m, p, t=0) {
    delims <- get_section_delims(slim_file)
    read_section_data("mutations", delims, slim_file) %>%
        read.table(text=., col.names=MUTATION_COLS, stringsAsFactors=FALSE) %>%
        filter(mut_type == m & pop_origin == p & gen_origin > t)
}

slim_file <- readLines("exome_and_sites__h_0.5__seed_6977220333793.txt")
nea_mut <- read_mutations(slim_file, "m0", "p2")


ggplot(mut, aes(pop_origin, log_s, fill=pop_origin)) +
    geom_violin() +
    geom_hline(yintercept=log10(1/(2 * c(1000, 10000))), color=c("green", "red")) +
    coord_flip() +
    ylim(-1, -12)

group_by(mut, pop_origin) %>% summarise(mean(s), median(s), max(s), min(s))
