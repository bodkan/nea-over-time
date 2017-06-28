#!/bin/bash

# Email from Felix Kay:
# # Turn the ompressed eigenstrat files from the Reich Webpage (genotype info for all 230ancient inds) into vcf
# Code from MI's github: https://github.com/mathii/gdc
# To get the code working I also downloaded gdc.py (same location) and pyEigenstrat.py (https://github.com/mathii/pyEigenstrat/blob/master/pyEigenstrat.py)
# NOTE: File misses shebang!

# python /home/felix_schulze/projects/trpm8_proj/mathiseson_data/eigenstrat2vcf.py -r /home/felix_schulze/projects/trpm8_proj/mathiseson_data/MathiesonEtAl_genotypes/full230 > /home/felix_schulze/projects/trpm8_proj/mathiseson_data/MathiesonEtAl_genotypes/full230.vcf
# 
# bgzip /home/felix_schulze/projects/trpm8_proj/mathiseson_data/MathiesonEtAl_genotypes/full230.vcf 
# tabix -p vcf /home/felix_schulze/projects/trpm8_proj/mathiseson_data/MathiesonEtAl_genotypes/full230.vcf.gz

if [ "$#" -ne 4 ]; then
    echo "This script accepts 4 arguments!"
    exit
fi

array_bed=$1
vcf=$2
mathieson_dir=$3
output_tsv=$4

# create header of the SNP table
echo "chrom pos ref alt" `cat ${mathieson_dir}/full230.ind | tr -s ' ' | cut -f2 -d ' ' | tr '\n' ' '` | tr ' ' '\t' > $output_tsv

# process the full230 VCF file
bcftools query -R $array_bed -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' $vcf \
    | sed 's/0\/0/0/g; s/0\/1/1/g; s/1\/1/2/g; s#\./\.#9#g' \
    >> $output_tsv
