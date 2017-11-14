#!/usr/bin/env python3

import sys

in_bed = open(sys.argv[1], 'r')

prev_chrom = ''
prev_end = ''

with open(sys.argv[2], 'w') as out_bed:
    for line in in_bed:
        chrom, start, end = line.split()

        if chrom == prev_chrom and prev_end != start:
            print(chrom, prev_end, start, 'spacer', sep='\t', file=out_bed)
        print(chrom, start, end, 'original', sep='\t', file=out_bed)

        prev_chrom, prev_end = chrom, end
