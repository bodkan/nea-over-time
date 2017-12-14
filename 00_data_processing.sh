# ---------------------------------------------------------------------- 
# convert the Ice Age genotypes from EIGENSTRAT into a normal genotype table
# and do the same with the SGDP VCF files and with archaic VCF files

eigenstrat_dir=data/eigenstrat
mkdir -p $eigenstrat_dir
cd $eigenstrat_dir

# download Fu et al. data
mkdir Fu
cd Fu
wget https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/FuQ.zip
unzip FuQ.zip
cd ../../../

# extract the coordinates of archaic admixture array sites (for later)
mkdir -p data/bed
awk -v OFS="\t" '{print $2, $4 - 1, $4}' ${eigenstrat_dir}/Fu/archaic.snp > data/bed/admixture_array.bed

gt_dir=data/genotypes
mkdir -p $gt_dir

# process Fu et al. data into a simple 0/1/2 table
Rscript code/process_eigenstrat.R ${eigenstrat_dir}/Fu/archaic.geno ${eigenstrat_dir}/Fu/archaic.snp ${eigenstrat_dir}/Fu/archaic.ind ${gt_dir}/ice_age.tsv
chmod -w ${gt_dir}/ice_age.tsv

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
awk -v OFS="\t" '{print $1, $2 - 1, $2}' data/big_yoruba_and_altai_filtN_printed_probes_with_annotation.txt > data/bed/big_yoruba_array.bed



# ---------------------------------------------------------------------- 
# generate updated EIGENSTRAT 2.2 M sites data that include the new
# processing of Altai, Vindija and Denisova (genotyped using snpAD)

# make a copy of Qiaomei's combined Eigenstrat dataset
mkdir data/eigenstrat/big_yoruba_ho; cd data/eigenstrat/big_yoruba_ho;
cp /mnt/454/Carbon_beast_QM/TY/snp/UPA_all.{snp,ind,geno} .

# make a directory for different coordinate files
less UPA_all.snp | awk -v OFS="\t" '{print $2, $4-1, $4}' > ../coordinates/2.2M.bed

# ---------------------------------------------------------------------- 
# generate the new high coverage Vindija "geno" and "snp" files

seq 1 22 | xargs -P 22 -I {} bash -c "bcftools view -R 2.2M.pos -M 2 /mnt/454/Vindija/high_cov/genotypes/Vindija33.19/chr{}_mq25_mapab100.vcf.gz | bcftools query -f '%CHROM\t%POS\t[%GT]\n' | sed 's/0\/0/2/g; s/0\/1/1/g; s/1\/1/0/g' > chr{}.tmp"
cat chr{1..22}.tmp > vindija.tmp
rm chr*.tmp

library(tidyverse)
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
genooutfilename: UPA_vi.geno
snpoutfilename: UPA_vi.snp
indoutfilename: UPA_vi.ind" > mergeit_vi.par

mergeit -p mergeit_vi.par


# ---------------------------------------------------------------------- 
# generate the new high coverage Altai "geno" and "snp" files
seq 1 22 | xargs -P 22 -I {} bash -c "bcftools view -R 2.2M.pos -M 2 /mnt/454/Vindija/high_cov/genotypes/Altai/chr{}_mq25_mapab100.vcf.gz | bcftools query -f '%CHROM\t%POS\t[%GT]\n' | sed 's/0\/0/2/g; s/0\/1/1/g; s/1\/1/0/g' > chr{}.tmp"
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
genooutfilename: UPA_vialtai.geno
snpoutfilename: UPA_vialtai.snp
indoutfilename: UPA_vialtai.ind" > mergeit_vialtai.par

mergeit -p mergeit_vialtai.par



# ---------------------------------------------------------------------- 
# generate the new high coverage Denisova "geno" and "snp" files
seq 1 22 | xargs -P 22 -I {} bash -c "bcftools view -R 2.2M.pos -M 2 /mnt/454/Vindija/high_cov/genotypes/Denisova/chr{}_mq25_mapab100.vcf.gz | bcftools query -f '%CHROM\t%POS\t[%GT]\n' | sed 's/0\/0/2/g; s/0\/1/1/g; s/1\/1/0/g' > chr{}.tmp"
cat chr{1..22}.tmp > denisova.tmp
rm chr*.tmp

library(tidyverse)
all <- read_table2("UPA_all.snp", col_names=c("id", "chrom", "gen", "pos", "alt", "ref"))
denisova <- read_tsv("denisova.tmp", col_names=c("chrom", "pos", "geno"))
merged <- left_join(all, denisova)
merged$geno[is.na(merged$geno)] <- 9
write_tsv(select(merged, -geno), "denisova.snp", col_names=FALSE)
write_tsv(select(merged, geno), "denisova.geno", col_names=FALSE)


# create 'ind' EIGENSTRAT file
echo "new_Denisova F new_Denisova" > denisova.ind

# generate a mergit parameter file
echo "outputformat: EIGENSTRAT
geno1: UPA_merged_all.geno
snp1: UPA_merged_all.snp
ind1: UPA_merged_all.ind
geno2: denisova.geno
snp2: denisova.snp
ind2: denisova.ind
genooutfilename: UPA_merged.geno
snpoutfilename: UPA_merged.snp
indoutfilename: UPA_merged.ind" > mergeit_merged.par

mergeit -p mergeit_merged.par




# # get the merged VCFs of SGDP + archaics
# mkdir raw_data/merged_vcfs

# for f in /mnt/scratch/steffi/D/Vcfs/mergedArchaics_SGDP0/*chr{1..22}.vcf.gz; do zcat $f | ./filter_vcf_with_bed.py <(awk -vOFS="\t" '{print $2, $4-1, $4}' raw_data/eigenstrat_all/UPA_merged_all.snp); done | awk '($5 != "." && $5 != "-" && length($5) == 1) { print $0 }' | bgzip > raw_data/merged_vcfs/mergedArchaics_qual0.vcf.gz

# for f in /mnt/scratch/steffi/D/Vcfs/mergedArchaics/*chr{1..22}.vcf.gz; do zcat $f | ./filter_vcf_with_bed.py <(awk -vOFS="\t" '{print $2, $4-1, $4}' raw_data/eigenstrat_all/UPA_merged_all.snp); done | awk '($5 != "." && $5 != "-" && length($5) == 1) { print $0 }' | bgzip > raw_data/merged_vcfs/mergedArchaics_qual1.vcf.gz






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
# GT archaics and africans table
echo -e "chrom\tstart\tend\tref\talt\tAltai\tVindija\tDenisova\tChimp\tMbuti\tYoruba" > data/genotypes/genload.bed
for chr in `seq 1 22`; do
    tabix /mnt/scratch/steffi/D/Vcfs/mergedArchaics/merged_archaics_manifesto_chr${chr}.vcf.gz \
        -R <( sed 's/chr//' data/bed/admixture_array.bed | cut -f 1,3) \
        | cut -f 1,2,4,5,10,11,12,13,196,302 \
        | grep -v "," \
        | sed 's/0\/0/0/g; s/0\/1/1/g; s/1\/1/2/g; s#\./\.# #g'
done \
    | awk -v OFS="\t" '{ $2=$2-1"\t"$2; print $0 }' \
    >> data/genotypes/genload.bed



less /mnt/expressions/benjamin_vernot/martin_neand_over_time_emh/1kg_p3_allele_freqs/YRI.freq \
    | tr -s ' ' | tr '\t' ' ' | awk -v OFS="\t" '{print $4, $1, $8}' | sort -k1,1n -k2,2n | bgzip \
    > data/YRI.freq.gz
