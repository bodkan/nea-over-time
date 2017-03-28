#!/usr/bin/env bash

# download Fu et al. data
mkdir -p raw_data
cd raw_data
wget https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/FuQ.zip
unzip FuQ.zip
cd ../

mkdir -p clean_data

# process Fu et al. data
./../../R/process_eigenstrat.R raw_data/archaic.geno raw_data/archaic.snp raw_data/archaic.ind clean_data/ice_age.tsv

# format the Altai/Vindija/Denisovan VCFs into a simple 0/1/2 table
echo -e "chrom\tpos\tref\talt\tAltai\tVindija\tDenisovan" > clean_data/archaics.tsv
for chrom in `seq 1 22`; do
    bcftools query -R <(tail -n+2 clean_data/ice_age.tsv) -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' /mnt/scratch/steffi/D/Vcfs/Altai_Vindija_Denis/altai_vindija_denis_chr${chrom}.vcf.gz \
	| sed 's/0\/0/0/g; s/0\/1/1/g; s/1\/1/2/g'
done >> clean_data/archaics.tsv

# get the list of SGDP (C-team)
sgdp_samples=`bcftools query -l /mnt/sequencedb/SGDP_May2016/combined_vcf/c_team_chr1.vcf.gz | grep "^S_" | tr '\n' ' '`
# process the SGDP VCF files to get the numbers of archaic-like
# alleles for each samples
echo "chrom pos ref alt "`echo -e $sgdp_samples` | tr ' ' '\t' | tr '-' '_' > clean_data/sgdp.tsv
time for chrom in `seq 1 22`; do
    bcftools query -R <(tail -n+2 clean_data/ice_age.tsv) -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' /mnt/sequencedb/SGDP_May2016/combined_vcf/c_team_chr${chrom}.vcf.gz -s `echo -e $sgdp_samples | tr ' ' ','` \
    | awk '$4 !~ "-"' \
    | sed 's/0\/0/0/g; s/0\/1/1/g; s/1\/1/2/g; s#\./\.#9#g'
done \
    | grep -v "/" >> clean_data/sgdp.tsv # remove remaining triallelic sites

# format the SGDP individuals into a simple 0/1/2 table
echo -e "chrom\tpos\tref\talt\tS_Dinka_1\tS_Dinka_2\tS_Mbuti_1\tS_Mbuti_2\tS_Mbuti_3\tS_Yoruba_1\tS_Yoruba_2\tS_French_1\tS_French_2\tS_Sardinian_1\tS_Sardinian_2\tS_Papuan_1\tS_Papuan_2\tS_Papuan_3" > clean_data/sgdp.tsv
for chrom in `seq 1 22`; do
    bcftools query -R <(tail -n+2 clean_data/ice_age.tsv) -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' /mnt/sequencedb/SGDP_May2016/combined_vcf/c_team_chr${chrom}.vcf.gz -s S_Dinka-1,S_Dinka-2,S_Mbuti-1,S_Mbuti-2,S_Mbuti-3,S_Yoruba-1,S_Yoruba-2,S_French-1,S_French-2,S_Sardinian-1,S_Sardinian-2,S_Papuan-1,S_Papuan-2,S_Papuan-3 \
        | awk '$4 !~ "-"' \
        | sed 's/0\/0/0/g; s/0\/1/1/g; s/1\/1/2/g; s#\./\.#9#g'
done \
    | grep -v "/" >> clean_data/sgdp.tsv # remove remaining triallelic sites

# annotate the SNPs using CADD
gunzip v1.3_anno*
sed 's/#//;2q;d' v1.3_anno*.tsv > annotations.tsv_unsorted
for file in v1.3_anno*.tsv; do
    tail -n+3 $file >> annotations.tsv_unsorted
done
sort -k1,1n -k2,2n annotations.tsv_unsorted > annotations.tsv
rm annotations.tsv_unsorted
