import pandas as pd
import numpy as np
import pybedtools
from pybedtools import BedTool
import argparse

parser = argparse.ArgumentParser(description='Generate annotations of exonic distances/densities')

parser.add_argument('--input-file', required=True,
                    help='File with SNP coordinates (needs to have chromosome '
                    'and position in its first two columns')
parser.add_argument('--gtf-file', metavar='FILE', required=True,
                    help='Genome annotation file')
parser.add_argument('--window-sizes', nargs='+', required=True, type=int,
                    help='Window sizes to calculate exonic densities for')
parser.add_argument('--output-file', required=True,
                    help='Output file name')

args = parser.parse_args()

snp_pos = pd.read_table(args.input_file)[['chrom', 'pos']]

snp_pos['start'] = snp_pos.pos - 1
snp_pos['end'] = snp_pos.pos
snp_bed = BedTool.from_dataframe(snp_pos[['chrom', 'start', 'end']]).sort()

gtf = pd.read_table(args.gtf_file,
                    header=None, sep='\t', skipinitialspace=True, skiprows=5, compression='gzip',
                    names=['chrom', 'source', 'feature', 'start', 'end',
                           'score', 'strand', 'frame', 'attribute'], low_memory=False)

# subset to exon annotations only and create a BED object
exons = BedTool.from_dataframe(
    gtf[(gtf.source == "protein_coding") &
        (gtf.feature == "exon")].query('end - start > 10')
).sort().merge().sort()

# calculate the distance of each SNP to the nearest exon
closest = snp_bed.closest(exons, t='first', d=True).to_dataframe()
snp_pos = snp_pos.merge(closest, on=["chrom", "start", "end"])[["chrom", "pos", "start", "end", "thickStart"]] \
                 .rename(columns={"thickStart": "exon_distance"})


# calculate the densities of exon sequences in windows surrounding each SNP
# flank specifies how far upstream/downstream to extend the window
# Window size will be (2 * flank + 1) bp
for flank in args.window_sizes:
    # generate the BED object of windows flanking both sides of all SNPs
    snp_windows = snp_bed.slop(b=flank, genome='hg19')
    
    # for each exon overlapping a given window, count the number of bases
    # of the overlap -- one row for each potential exon
    # If there are no exons overlapping a window, report 0 bp overlap
    exon_overlaps = snp_windows.intersect(exons, wao=True)   \
                               .to_dataframe(low_memory=False)[['chrom', 'start', 'end', 'thickStart']]
        
    # thickStart contains the number of bases overlapping each exon in
    # a window around each SNP -- Summing up this column for each unique
    # region using groupby (i.e. window around each SNP) gives the total
    # number of coding sequence surrounding each SNP
    exon_total = exon_overlaps.groupby(['chrom', 'start', 'end'])['thickStart'].sum().reset_index()['thickStart']
    snp_pos['exon_density_' + str(flank)] = exon_total / (2 * flank + 1)

snp_pos \
    .drop(["pos"], axis=1) \
    .to_csv(args.output_file, sep='\t', index=False, na_rep='-')
