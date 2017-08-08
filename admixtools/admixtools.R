# eurasian lineage was there
library(tidyverse)
library(stringr)

# read a table of F4 ratio alpha estimates
read_f4_ratios <- function(log_filename) {
    log_lines <- readLines(log_filename)

    # extract the number of analyzed test populations/individuals
    # (corresponding to the number of rows of the results table)
    n_pops <- log_lines[which(str_detect(log_lines, "^nplist:"))] %>%
        str_extract("[0-9]+$") %>%
        as.integer

    # parse the lines of the results section and extract the names of
    # tested populations/individuals, estimated admixture proportions
    # alpha, std. errors and Z-score
    res <- log_lines[(length(log_lines) - n_pops) : (length(log_lines) - 1)] %>%
        str_replace_all(":", "") %>%
        str_replace_all(" +", " ") %>%
        str_replace("result", "") %>%
        paste(collapse="\n") %>%
        textConnection %>%
        read.table %>%
        setNames(c("A", "O", "X", "C", "A", "O", "B", "C", "alpha", "stderr", "z")) %>%
        .[c("A", "B", "X", "C", "O", "alpha", "stderr", "z")]

    res
}

# generate a vector of poplist entries, optionally saving it to a file
create_f4_poplist_file <- function(X, A, B, C, O, filename) {
    lines <- sprintf("%s %s : %s %s :: %s %s : %s %s", A, O, X, C, A, O, B, C)
    writeLines(lines, filename)
}

create_f4_param_file <- function(param_file, poplist_file,
                                 eigenstrat_prefix=NULL,
                                 geno_file=NULL, snp_file=NULL, ind_file= NULL) {
    if (is.null(eigenstrat_prefix) & all(is.null(c(geno_file, snp_file, ind_file)))) {
        stop("Either the eigenstrat_prefix or paths to geno/snp/ind files must be specified")
    }

    if (!is.null(eigenstrat_prefix)) {
        geno_file <- paste0(eigenstrat_prefix, ".geno")
        snp_file <- paste0(eigenstrat_prefix, ".snp")
        ind_file <- paste0(eigenstrat_prefix, ".ind")
    }

    writeLines(sprintf("genotypename: %s\nsnpname: %s\nindivname: %s\npopfilename: %s",
                       geno_file, snp_file, ind_file, poplist_file),
               con=param_file)
}

run_f4 <- function(param_file, log_file, admixtools_path="/Users/martin_petr/local/Admixtools/bin/") {
    system(command=paste(file.path(admixtools_path, "qpF4ratio"),
                         "-p",
                         param_file,
                         ">",
                         log_file))
}

