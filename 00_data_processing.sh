#!/usr/bin/env bash

eigenstrat_dir=data/eigenstrat
mkdir -p $eigenstrat_dir

# ---------------------------------------------------------------------- 
# process merged SGDP VCF files
vcf_dir=data/vcf
mkdir -p $vcf_dir

sample_ids="panTro4,AltaiNeandertal,Vindija33.19,Denisova,Ust_Ishim,Loschbour,LBK,`cut -f1-4 data/10_24_2014_SGDP_metainformation_update.txt | egrep "Africa|WestEurasia|EastAsia" | grep "^C" | cut -f2 | tr '\n' ',' | sed 's/,$//'`"
bcftools view -h -s $sample_ids /mnt/scratch/steffi/D/Vcfs/mergedArchModernApes_new/merged_ancient_apes_sgdp1_chr1.vcf.gz | sed 's/panTro4/Chimp/; s/Ust_Ishim/new_UstIshim/; s/Loschbour/new_Loschbour/; s/LBK/new_Stuttgart/' \
    > ${vcf_dir}/whole_genome.vcf

for chrom in `seq 1 22`; do
    bcftools view -a -s $sample_ids /mnt/scratch/steffi/D/Vcfs/mergedArchModernApes_new/merged_ancient_apes_sgdp1_chr${chrom}.vcf.gz \
        | bcftools view -H -m2 -M2 -v snps \
        | awk '$5 != "-"' >> ${vcf_dir}/whole_genome.vcf
done

bgzip ${vcf_dir}/whole_genome.vcf
tabix ${vcf_dir}/whole_genome.vcf.gz

mkdir -p ${eigenstrat_dir}/whole_genome
python3 code/vcf_to_eigenstrat.py ${vcf_dir}/whole_genome.vcf.gz ${eigenstrat_dir}/whole_genome/whole_genome

# use the whole-genome VCF to create EIGENSTRAT format data of 2.2M sites
python3 code/vcf_to_eigenstrat.py <(bcftools view -R data/bed/2.2M.bed ${vcf_dir}/whole_genome.vcf.gz -Oz) \
    ${eigenstrat_dir}/whole_genome/2.2M




# ---------------------------------------------------------------------- 
# convert the Ice Age genotypes from EIGENSTRAT into a normal genotype table
# and do the same with the SGDP VCF files and with archaic VCF files

cd $eigenstrat_dir

# download Fu et al. data
mkdir iceage
cd iceage
wget https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/FuQ.zip
unzip FuQ.zip
cd ../../../

# extract the coordinates of archaic admixture array sites (for later)
mkdir -p data/bed
awk -v OFS="\t" '{print $2, $4 - 1, $4}' ${eigenstrat_dir}/iceage/archaic.snp > data/bed/admixture_array.bed

gt_dir=data/genotypes
mkdir -p $gt_dir

# process Fu et al. data into a simple 0/1/2 table
Rscript code/process_eigenstrat.R ${eigenstrat_dir}/iceage/archaic.geno ${eigenstrat_dir}/iceage/archaic.snp ${eigenstrat_dir}/iceage/archaic.ind ${gt_dir}/ice_age.tsv
chmod -w ${gt_dir}/ice_age.tsv





########
# fix this to use the newest merge of archaics & SGDP data from Steffi
#
# TODO - now she has proper VCF files so I can use bcftools to subset
# to the positions from the admixture array (without writing a custom
# Python script to filter for non-Chimp-unique sites)

# get the list of SGDP (C-team)
sgdp_samples=`bcftools query -l /mnt/sequencedb/SGDP_May2016/combined_vcf/c_team_chr1.vcf.gz | grep "^S_" | tr '\n' ' '`
# process the SGDP VCF files to get the numbers of archaic-like
# alleles for each samples
echo "chrom pos ref alt "`echo -e $sgdp_samples` | tr ' ' '\t' | tr '-' '_' > ${gt_dir}/sgdp.tsv
for chrom in `seq 1 22`; do
    bcftools query -R <(tail -n+2 ${gt_dir}/ice_age.tsv) -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' /mnt/sequencedb/SGDP_May2016/combined_vcf/c_team_chr${chrom}.vcf.gz -s `echo -e $sgdp_samples | tr ' ' ','` \
    | awk '$4 !~ "-"' \
    | sed 's/0\/0/0/g; s/0\/1/1/g; s/1\/1/2/g; s#\./\.#9#g'
done \
    | grep -v "," >> ${gt_dir}/sgdp.tsv # remove remaining triallelic sites
chmod -w ${gt_dir}/sgdp.tsv

# format the Altai/Vindija/Denisovan VCFs into a simple 0/1/2 table
echo -e "chrom\tpos\tref\talt\tAltai\tVindija\tDenisovan" > ${gt_dir}/archaics.tsv
for chrom in `seq 1 22`; do
    bcftools query -R <(tail -n+2 ${gt_dir}/ice_age.tsv) -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' /mnt/scratch/steffi/D/Vcfs/Altai_Vindija_Denis/altai_vindija_denis_chr${chrom}_filtered_derived.gz \
	| sed 's/0\/0/0/g; s/0\/1/1/g; s/1\/1/2/g'
done >> ${gt_dir}/archaics.tsv
chmod -w ${gt_dir}/archaics.tsv





# ---------------------------------------------------------------------- 
#  get the coordinates of the Big Yoruba array sites
cp /r1/people/public/AncientDNA/probe_designs/AA77-81_bigYoruba/big_yoruba_and_altai_filtN_printed_probes_with_annotation.txt data/

cat data/big_yoruba_and_altai_filtN_printed_probes_with_annotation.txt | awk '$12 == "0/1"' | awk '{print $1, $2 - 1, $2}' > data/bed/bigyri_Altai.bed
cat data/big_yoruba_and_altai_filtN_printed_probes_with_annotation.txt | awk '$6 == "0/1" || $9 == "0/1"' | awk '{print $1, $2 - 1, $2}' > data/bed/bigyri_YRI.bed
cat data/big_yoruba_and_altai_filtN_printed_probes_with_annotation.txt | awk -v OFS="\t" '{print $1, $2 - 1, $2}' data/big_yoruba_and_altai_filtN_printed_probes_with_annotation.txt > data/bed/bigyri_both.bed

# ---------------------------------------------------------------------- 
# generate updated EIGENSTRAT 2.2 M sites data that include the new
# processing of Altai and Vindija (genotyped using snpAD)

mkdir data/eigenstrat/bigyri_ho
cd data/eigenstrat/bigyri_ho

# make a copy of Qiaomei's combined Eigenstrat dataset
cp /mnt/454/Carbon_beast_QM/TY/snp/UPA_all.{snp,ind,geno} .

# save the coordinates of all 2.2M sites
less UPA_all.snp | awk -v OFS="\t" '{print $2, $4-1, $4}' > ../../bed/2.2M.bed

# ------------------------------
# VCF processing

# Vindija
seq 1 22 | xargs -P 22 -I {} bash -c "bcftools view -R ../../bed/2.2M.bed -M 2 /mnt/454/Vindija/high_cov/genotypes/Vindija33.19/chr{}_mq25_mapab100.vcf.gz | bcftools query -f '%CHROM\t%POS\t[%GT]\n' | sed 's/0\/0/2/g; s/0\/1/1/g; s/1\/1/0/g' > vindija_chr{}.tmp"
cat vindija_chr{1..22}.tmp > vindija.tmp
rm vindija_chr*.tmp

# Altai
seq 1 22 | xargs -P 22 -I {} bash -c "bcftools view -R ../../bed/2.2M.bed -M 2 /mnt/454/Vindija/high_cov/genotypes/Altai/chr{}_mq25_mapab100.vcf.gz | bcftools query -f '%CHROM\t%POS\t[%GT]\n' | sed 's/0\/0/2/g; s/0\/1/1/g; s/1\/1/0/g' > altai_chr{}.tmp"
cat altai_chr{1..22}.tmp > altai.tmp
rm altai_chr*.tmp

# Denisova
seq 1 22 | xargs -P 22 -I {} bash -c "bcftools view -R ../../bed/2.2M.bed -M 2 /mnt/454/Vindija/high_cov/genotypes/Denisova/chr{}_mq25_mapab100.vcf.gz | bcftools query -f '%CHROM\t%POS\t[%GT]\n' | sed 's/0\/0/2/g; s/0\/1/1/g; s/1\/1/0/g' > denisova_chr{}.tmp"
cat denisova_chr{1..22}.tmp > denisova.tmp
rm denisova_chr*.tmp

# UstIshim
seq 1 22 | xargs -P 22 -I {} bash -c "bcftools view -R ../../bed/2.2M.bed -M 2 /mnt/454/Vindija/high_cov/genotypes/Ust_Ishim/chr{}_mq25_mapab100.vcf.gz | bcftools query -f '%CHROM\t%POS\t[%GT]\n' | sed 's/0\/0/2/g; s/0\/1/1/g; s/1\/1/0/g' > ustishim_chr{}.tmp"
cat ustishim_chr{1..22}.tmp > ustishim.tmp
rm ustishim_chr*.tmp

# Loschbour
seq 1 22 | xargs -P 22 -I {} bash -c "bcftools view -R ../../bed/2.2M.bed -M 2 /mnt/454/Vindija/high_cov/genotypes/Loschbour/chr{}_mq25_mapab100.vcf.gz | bcftools query -f '%CHROM\t%POS\t[%GT]\n' | sed 's/0\/0/2/g; s/0\/1/1/g; s/1\/1/0/g' > loschbour_chr{}.tmp"
cat loschbour_chr{1..22}.tmp > loschbour.tmp
rm loschbour_chr*.tmp

# ------------------------------
# conversion to EIGENSTRAT

# Vindija
R --no-save < <(echo '
library(tidyverse)
all <- read_table2("UPA_all.snp", col_names=c("id", "chrom", "gen", "pos", "alt", "ref"))
vindija <- read_tsv("vindija.tmp", col_names=c("chrom", "pos", "geno"))
merged <- left_join(all, vindija)
merged$geno[is.na(merged$geno)] <- 9
write_tsv(select(merged, -geno), "vindija.snp", col_names=FALSE)
write_tsv(select(merged, geno), "vindija.geno", col_names=FALSE)
')

# Altai
R --no-save < <(echo '
library(tidyverse)
all <- read_table2("UPA_all.snp", col_names=c("id", "chrom", "gen", "pos", "alt", "ref"))
altai <- read_tsv("altai.tmp", col_names=c("chrom", "pos", "geno"))
merged <- left_join(all, altai)
merged$geno[is.na(merged$geno)] <- 9
write_tsv(select(merged, -geno), "altai.snp", col_names=FALSE)
write_tsv(select(merged, geno), "altai.geno", col_names=FALSE)
')

# Denisova
R --no-save < <(echo '
library(tidyverse)
all <- read_table2("UPA_all.snp", col_names=c("id", "chrom", "gen", "pos", "alt", "ref"))
denisova <- read_tsv("denisova.tmp", col_names=c("chrom", "pos", "geno"))
merged <- left_join(all, denisova)
merged$geno[is.na(merged$geno)] <- 9
write_tsv(select(merged, -geno), "denisova.snp", col_names=FALSE)
write_tsv(select(merged, geno), "denisova.geno", col_names=FALSE)
')

# UstIshim
R --no-save < <(echo '
library(tidyverse)
all <- read_table2("UPA_all.snp", col_names=c("id", "chrom", "gen", "pos", "alt", "ref"))
tmp <- read_tsv("ustishim.tmp", col_names=c("chrom", "pos", "geno"))
merged <- left_join(all, tmp)
merged$geno[is.na(merged$geno)] <- 9
write_tsv(select(merged, -geno), "ustishim.snp", col_names=FALSE)
write_tsv(select(merged, geno), "ustishim.geno", col_names=FALSE)
')

# Loschbour
R --no-save < <(echo '
library(tidyverse)
all <- read_table2("UPA_all.snp", col_names=c("id", "chrom", "gen", "pos", "alt", "ref"))
tmp <- read_tsv("loschbour.tmp", col_names=c("chrom", "pos", "geno"))
merged <- left_join(all, tmp)
merged$geno[is.na(merged$geno)] <- 9
write_tsv(select(merged, -geno), "loschbour.snp", col_names=FALSE)
write_tsv(select(merged, geno), "loschbour.geno", col_names=FALSE)
')

# ------------------------------
# mergeit runs

# merging Vindija

# create 'ind' EIGENSTRAT file
echo "new_Vindija F new_Vindija" > vindija.ind
# generate a mergit parameter file
echo "outputformat: EIGENSTRAT
geno1: UPA_all.geno
snp1: UPA_all.snp
ind1: UPA_all.ind
geno2: vindija.geno
snp2: vindija.snp
ind2: vindija.ind
genooutfilename: UPA_Vindija.geno
snpoutfilename: UPA_Vindija.snp
indoutfilename: UPA_Vindija.ind" > mergeit_Vindija.par
mergeit -p mergeit_Vindija.par

# merging Altai

# create 'ind' EIGENSTRAT file
echo "new_Altai F new_Altai" > altai.ind
# generate a mergit parameter file
echo "outputformat: EIGENSTRAT
geno1: UPA_Vindija.geno
snp1: UPA_Vindija.snp
ind1: UPA_Vindija.ind
geno2: altai.geno
snp2: altai.snp
ind2: altai.ind
genooutfilename: UPA_Altai.geno
snpoutfilename: UPA_Altai.snp
indoutfilename: UPA_Altai.ind" > mergeit_Altai.par
mergeit -p mergeit_Altai.par

# merging Denisova

# create 'ind' EIGENSTRAT file
echo "new_Denisova F new_Denisova" > denisova.ind
# generate a mergit parameter file
echo "outputformat: EIGENSTRAT
geno1: UPA_Altai.geno
snp1: UPA_Altai.snp
ind1: UPA_Altai.ind
geno2: denisova.geno
snp2: denisova.snp
ind2: denisova.ind
genooutfilename: UPA_Denisova.geno
snpoutfilename: UPA_Denisova.snp
indoutfilename: UPA_Denisova.ind" > mergeit_Denisova.par
mergeit -p mergeit_Denisova.par

# merging UstIshim

# create 'ind' EIGENSTRAT file
echo "new_UstIshim F new_UstIshim" > ustishim.ind
# generate a mergit parameter file
echo "outputformat: EIGENSTRAT
geno1: UPA_Denisova.geno
snp1: UPA_Denisova.snp
ind1: UPA_Denisova.ind
geno2: ustishim.geno
snp2: ustishim.snp
ind2: ustishim.ind
genooutfilename: UPA_UstIshim.geno
snpoutfilename: UPA_UstIshim.snp
indoutfilename: UPA_UstIshim.ind" > mergeit_UstIshim.par
mergeit -p mergeit_UstIshim.par

# merging Loschbour

# create 'ind' EIGENSTRAT file
echo "new_Loschbour F new_Loschbour" > loschbour.ind
# generate a mergit parameter file
echo "outputformat: EIGENSTRAT
geno1: UPA_UstIshim.geno
snp1: UPA_UstIshim.snp
ind1: UPA_UstIshim.ind
geno2: loschbour.geno
snp2: loschbour.snp
ind2: loschbour.ind
genooutfilename: all.geno
snpoutfilename: all.snp
indoutfilename: all.ind" > mergeit_Loschbour.par
mergeit -p mergeit_Loschbour.par




# ---------------------------------------------------------------------- 
# download the McVicker B values
cd data
wget http://www.phrap.org/software_dir/mcvicker_dir/bkgd.tar.gz
tar xvf bkgd.tar.gz
rm bkgd.tar.gz





# ---------------------------------------------------------------------- 
# download the hg18-to-hg19 liftover chain
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
gunzip hg18ToHg19.over.chain.gz





# ---------------------------------------------------------------------- 
# YRI admixture array frequencies from Ben's table
less /mnt/expressions/benjamin_vernot/martin_neand_over_time_emh/1kg_p3_allele_freqs/YRI.freq \
    | tr -s ' ' | tr '\t' ' ' | awk -v OFS="\t" '{print $4, $1, $8}' | sort -k1,1n -k2,2n | bgzip \
    > data/YRI.freq.gz
