import argparse

import pandas as pd
from pybedtools import BedTool
 

def create_bed(snps):
    """Convert a Pandas table with chrom and pos columns into
    a BedTool object.
    """
    snps["start"] = snps.pos - 1
    snps["end"] = snps.pos

    return BedTool.from_dataframe(snps[["chrom", "start", "end"]])


def get_genic_coords():
    """Download the GTF annotation file and extract the coordinates of genic regions."""
    gtf = pd.read_table(
        'ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz',
        header=None, sep='\t', skipinitialspace=True, skiprows=5, compression='gzip',
        names=['chrom', 'source', 'feature', 'start', 'end',
               'score', 'strand', 'frame', 'attribute'], low_memory=False
    )
    genic = BedTool.from_dataframe(gtf[gtf.source == "protein_coding"]).sort().merge()

    return genic


def main(argv=None):
    parser = argparse.ArgumentParser(description="Add window-averaged"
                                     "PhyloP as a column of a given table")
    parser.add_argument('--input-table', metavar='FILE', required=True,
                        help='Tab delimited text file with SNP positions')
    Parser.add_argument('--genic-window', required=True, type=int,
                        help='Size of the window around each SNP for genic density calculation')

    # if there were no arguments supplied to the main function, use sys.argv
    # (skipping the first element, i.e. the name of this script)
    args = parser.parse_args(argv if argv else sys.argv[1:])

    pos = pd.read_table(args.input_table, usecols=["chrom", "pos"])
    
    return pos


#if __name__ == "__main__":
snp_table = main("--input-table clean_data/ice_age.tsv --genic-window 10000".split())
snps = create_bed(snp_table).sort()
genes = get_genic_coords()

# distance to genes
snp_table["gene_distance"] = snps.closest(genes, t='first', d=True).to_dataframe()[[-1]].values

# genic density
# generate the BED object of windows flanking both sides of all SNPs
snp_windows = snps.slop(b=int(args.genic_window / 2), genome='hg19')

# for each exon overlapping a given window, count the number of bases
# of the overlap, one row for each potential exon - if there are no
# exons overlapping a window, report 0 bp overlap
exon_overlaps = snp_windows.intersect(genes, wao=True)   \
                           .to_dataframe()[['chrom', 'start', 'end', 'thickStart']]

# thickStart contains the number of genic bases in a window around
# each SNP - summing up these values for each window (i.e. unique
# region by calling group_by) gives the total number of coding
# sequence surrounding each SNP
total_overlaps = exon_overlaps.groupby(['chrom', 'start', 'end'])['thickStart'] \
                              .sum() \
                              .reset_index()['thickStart']
snp_table['exon_overlap'] = total_overlaps / (2 * args.genic_window + 1)
