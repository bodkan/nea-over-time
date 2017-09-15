#!/usr/bin/env bash

# This script processes the published data from the Ice Age paper,
# SGDP VCF files and also VCF files with Altai, Vindija and Denisova genotypes.
# Finally, it processes the CADD functional annotation data of SNPs
# from the archaic admixture array.





mkdir -p raw_data
cd raw_data




# download Fu et al. data
mkdir Fu
cd Fu
wget https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/FuQ.zip
unzip FuQ.zip
cd ../../




mkdir -p clean_data


# process Fu et al. data
./R/process_eigenstrat.R raw_data/Fu/archaic.geno raw_data/Fu/archaic.snp raw_data/Fu/archaic.ind clean_data/ice_age.tsv
chmod -w clean_data/ice_age.tsv

# format the Altai/Vindija/Denisovan VCFs into a simple 0/1/2 table
echo -e "chrom\tpos\tref\talt\tAltai\tVindija\tDenisovan" > clean_data/archaics.tsv
for chrom in `seq 1 22`; do
    bcftools query -R <(tail -n+2 clean_data/ice_age.tsv) -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' /mnt/scratch/steffi/D/Vcfs/Altai_Vindija_Denis/altai_vindija_denis_chr${chrom}.vcf.gz \
	| sed 's/0\/0/0/g; s/0\/1/1/g; s/1\/1/2/g'
done >> clean_data/archaics.tsv
chmod -w clean_data/archaics.tsv

# get the list of SGDP (C-team)
sgdp_samples=`bcftools query -l /mnt/sequencedb/SGDP_May2016/combined_vcf/c_team_chr1.vcf.gz | grep "^S_" | tr '\n' ' '`
# process the SGDP VCF files to get the numbers of archaic-like
# alleles for each samples
echo "chrom pos ref alt "`echo -e $sgdp_samples` | tr ' ' '\t' | tr '-' '_' > clean_data/sgdp.tsv
for chrom in `seq 1 22`; do
    bcftools query -R <(tail -n+2 clean_data/ice_age.tsv) -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' /mnt/sequencedb/SGDP_May2016/combined_vcf/c_team_chr${chrom}.vcf.gz -s `echo -e $sgdp_samples | tr ' ' ','` \
    | awk '$4 !~ "-"' \
    | sed 's/0\/0/0/g; s/0\/1/1/g; s/1\/1/2/g; s#\./\.#9#g'
done \
    | grep -v "," >> clean_data/sgdp.tsv # remove remaining triallelic sites
chmod -w clean_data/sgdp.tsv






# make a copy of Qiaomei's combined Eigenstrat dataset
mkdir raw_data/eigenstrat_all; cd raw_data/eigenstrat_all
cp /mnt/454/Carbon_beast_QM/TY/snp/UPA_all.{snp,ind,geno} .

less UPA_all.snp | tr -s ' ' | cut -d ' ' -f 3,5 | tr ' ' '\t' > 2.2M.pos





# generate the new high coverage Vindija "geno" and "snp" files
seq 1 22 | xargs -P 22 -I {} bash -c "bcftools view -R 2.2M.pos -M 2 /mnt/454/Vindija/high_cov/genotypes/Vindija33.19/chr{}_mq25_mapab100.vcf.gz | bcftools query -f '%CHROM\t%POS\t[%GT]\n' | sed 's/0\/0/0/g; s/0\/1/1/g; s/1\/1/2/g' > chr{}.tmp"
cat chr{1..22}.tmp > vindija.tmp
rm chr*.tmp

library(tidyverse)
source("../../R/admixr.R")
all <- read_table2("UPA_all.snp", col_names=c("id", "chrom", "gen", "pos", "alt", "ref"))
vindija <- read_tsv("vindija.tmp", col_names=c("chrom", "pos", "geno"))
merged <- left_join(all, vindija)
merged$geno[is.na(merged$geno)] <- 9
write_tsv(select(merged, -geno), "vindija.snp", col_names=FALSE)
write_tsv(select(merged, geno), "vindija.geno", col_names=FALSE)


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
genooutfilename: UPA_merged.geno
snpoutfilename: UPA_merged.snp
indoutfilename: UPA_merged.ind" > mergeit.par

mergeit -p mergeit.par





# generate the new high coverage Altai "geno" and "snp" files
seq 1 22 | xargs -P 22 -I {} bash -c "bcftools view -R 2.2M.pos -M 2 /mnt/454/Vindija/high_cov/genotypes/Altai/chr{}_mq25_mapab100.vcf.gz | bcftools query -f '%CHROM\t%POS\t[%GT]\n' | sed 's/0\/0/0/g; s/0\/1/1/g; s/1\/1/2/g' > chr{}.tmp"
cat chr{1..22}.tmp > altai.tmp
rm chr*.tmp

library(tidyverse)
all <- read_table2("UPA_all.snp", col_names=c("id", "chrom", "gen", "pos", "alt", "ref"))
altai <- read_tsv("altai.tmp", col_names=c("chrom", "pos", "geno"))
merged <- left_join(all, altai)
merged$geno[is.na(merged$geno)] <- 9
write_tsv(select(merged, -geno), "altai.snp", col_names=FALSE)
write_tsv(select(merged, geno), "altai.geno", col_names=FALSE)


# create 'ind' EIGENSTRAT file
echo "new_Altai F new_Altai" > altai.ind

# generate a mergit parameter file
echo "outputformat: EIGENSTRAT
geno1: UPA_merged.geno
snp1: UPA_merged.snp
ind1: UPA_merged.ind
geno2: altai.geno
snp2: altai.snp
ind2: altai.ind
genooutfilename: UPA_merged_all.geno
snpoutfilename: UPA_merged_all.snp
indoutfilename: UPA_merged_all.ind" > mergeit_all.par

mergeit -p mergeit_all.par

















# annotate the SNPs using CADD
# the SNPs were annotated "manually" by uploading the ice_age.tsv pseudo-VCF file
# to the CADD online server and results were downloaded to raw_data/annotations
cd raw_data
mkdir annotations
cd annotations
gunzip v1.3_anno*
sed 's/#//;2q;d' v1.3_anno*.tsv > annotations.tsv_unsorted
for file in v1.3_anno*.tsv; do
    tail -n+3 $file >> annotations.tsv_unsorted
done
sort -k1,1n -k2,2n annotations.tsv_unsorted > ../../clean_data/annotations.tsv
cd ../../
chmod -w clean_data/annotations.tsv


# download the SGDP metainformation table
cd raw_data/
curl -O http://simonsfoundation.s3.amazonaws.com/share/SCDA/datasets/10_24_2014_SGDP_metainformation_update.txt
cd ../





######################################################################
## conservation annotations
######################################################################

annotation_dir=clean_data/annotations
mkdir -p $annotation_dir

# download the GTF file and annotate SNPs with distances to the nearest
# exon and densities of exons in different windows
curl ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz -o raw_data/annotations/gtf.gz
python3 exon_annotations.py --input-file clean_data/ice_age.tsv --gtf-file raw_data/annotations/gtf.gz --window-sizes 10000 25000 50000 100000 200000 --output-dir $annotation_dir

cadd_file=${annotation_dir}/genome_wide_cadd.bed.gz
akey_file=${annotation_dir}/genome_wide_phylop_akey.bed.gz

# generate whole-genome BED files with different annotations for
# later window-based analysis
# extract chrom, pos, priPhCons, mamPhCons, priPhyloP, mamPhyloP, bStatistic
# column numbers from here: http://cadd.gs.washington.edu/static/ReleaseNotes_CADD_v1.3.pdf
for chr in 1 {10..19} 2 20 21 22 {3..9}; do
    tabix /mnt/expressions/cadd/whole_genome_SNVs_inclAnno.tsv.gz ${chr} \
        | awk -vOFS="\t" '{ $1=$1"\t"$2-1; print $1, $2, $19, $20, $22, $23, $29 }' \
        | uniq
done | gzip > $cadd_file

# concatenate the phyloP data generated in the Akey lab
for chr in 1 {10..19} 2 20 21 22 {3..9}; do
    sed 's/^chr//' /mnt/expressions/benjamin_vernot/phyloP_no_human/Compara.36_eutherian_mammals_EPO_LOW_COVERAGE.chr${chr}_*
done | gzip > $akey_file

# calculate conservation statistics for different tracks and different window sizes
window_average() {
    track=$1      # track name
    column=$2     # track's column in the "annot_file"
    annot_file=$4 # BED file with the conservation track
    snp_file=$5   # table with SNP coordinates (needs to have 2 columns: chrom and pos)
    flank_size=`echo $3 / 2 | bc`

    zcat $annot_file \
	| cut -f1-3,${column} \
	| awk -vOFS='\t' '{ print $1, $2, $3, $4, $4 }' \
	| grep -v "NA" \
	| bedmap --ec --delim '\t' --range $flank_size --echo --mean --kth .05 --kth .95 --median --count \
		 <(tail -n+2 $snp_file | awk -vOFS='\t' '{print $1, $2-1, $2}' | sort-bed -) - \
        > $annotation_dir/${track}__${3}bp.bed
}

snp_file=clean_data/ice_age.tsv

for window_size in 10000 25000 50000 100000; do
    window_average priPhCons      4 $window_size $cadd_file $snp_file &
    window_average priPhyloP      6 $window_size $cadd_file $snp_file &
    window_average bval           8 $window_size $cadd_file $snp_file &
    window_average phyloP_nohuman 4 $window_size $akey_file $snp_file &
done







# download the McVicker B values
cd raw_data
wget http://www.phrap.org/software_dir/mcvicker_dir/bkgd.tar.gz
tar xvf bkgd.tar.gz
rm bkgd.tar.gz
# download the hg18-to-hg19 liftover chain
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
gunzip hg18ToHg19.over.chain.gz

