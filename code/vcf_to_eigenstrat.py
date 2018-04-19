import re
import sys
from collections import Counter
from itertools import chain
import gzip

trans_tbl = {
    "./." : 9,
    "0/0" : 2,
    "1/1" : 0,
    "1/0" : 1,
    "0/1" : 1
}

if len(sys.argv) != 3:
    print("Insufficient number of arguments")
    sys.exit(1)

ind, snp, geno = sys.argv[2] + ".ind", sys.argv[2] + ".snp", sys.argv[2] + ".geno"

with gzip.open(sys.argv[1], "rt") as vcf_file, \
           open(ind, "w") as ind_file, \
           open(snp, "w") as snp_file, \
           open(geno, "w") as geno_file:
    for line in vcf_file:
        if line.startswith("##"): continue
        
        fields = line.strip().split("\t")
        if line.startswith("#CHROM"):
            samples = fields[9:]
            for s in samples:
                if s.startswith("Altai"): pop = "new_Altai"
                elif s.startswith("Vindija"): pop = "new_Vindija"
                elif s.startswith("Denisova"): pop = "new_Denisova"
                else:
                    pop = re.sub("\d+|-|S_", "", s)
                print(s, "F", pop, file=ind_file)

        else:
            gts = fields[9:]
            chimp_gt = gts[0]
            hum_gts = Counter(gts[1:])

            if chimp_gt != "./.":
                chimp_allele = chimp_gt[0]
                hum_alleles = set(chain.from_iterable(gt.split("/") for gt in (hum_gts.keys())))
                hum_alleles.discard('9')

                if chimp_allele not in hum_alleles:
                    continue

            # write a record to a geno file
            print("".join(str(trans_tbl[gt]) for gt in gts), file=geno_file)

            # write a record to a snp file
            print("{}_{}\t{}\t0\t{}\t{}\t{}" \
                    .format(fields[0], fields[1], fields[0], fields[1], fields[3], fields[4]),
                  file=snp_file)

