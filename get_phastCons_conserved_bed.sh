#!/usr/bin/env bash

# This script converts a given gzipped wig file with phastCons track
# values into a BED file, converts the values of the conserved
# posterior probability into "+" (conserved) and "-" (not conserved)
# based on a specified cutoff and outputs the coordinates of the
# conserved elements in a BED format.
#
# phastConst from primate alignments can be downloaded from here:
# http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/

if [ "$#" -ne 3 ]; then
    echo "Usage: ./get_phastCons_conserved_bed.sh prob_cutoff wiggz_file output_bed"
    exit
fi

prob_cutoff=$1
wiggz_file=$2
output_bed=$3

gunzip -c $wiggz_file \
    | wig2bed -       \
    | awk -v OFS="\t" -v prob_cutoff=$prob_cutoff '
        {
            if ($5 >= prob_cutoff) { $5="+" } else {$5="-"};
            print $1, $2, $3, $4, $4, $5;
        }' \
    | bedtools merge -i - -S + \
    > $output_bed
