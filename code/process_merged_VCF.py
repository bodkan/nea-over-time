import sys
import gzip

import pysam

f = "/mnt/scratch/steffi/D/Vcfs/mergedArchaics/merged_archaics_manifesto_chr22.vcf.gz"
vcf = gzip.open(f, "rt")#sys.argv[1])
bed = open("/mnt/scratch/mp/nea-over-time/data/bed/admixture_array.bed", "r") #sys.argv[2], "r")

# read the coordinates of the BED file
bed_coords = set()
for bed_line in bed:
    bed_fields = bed.readline().rstrip().split()
    bed_coords.add((bed_fields[0], bed_fields[1], bed_fields[2]))

header = vcf.readline().rstrip().split()

vcf_line = vcf.readline().rstrip()
vcf_fields = vcf_line.split()
