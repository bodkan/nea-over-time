# ---------------------------------------------------------------------- 
# convert the Ice Age genotypes from EIGENSTRAT into a normal genotype table
# and do the same with the SGDP VCF files and with archaic VCF files

eigenstrat_dir=data/eigenstrat
mkdir -p $eigenstrat_dir
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
awk -v OFS="\t" '{print $1, $2 - 1, $2}' data/big_yoruba_and_altai_filtN_printed_probes_with_annotation.txt > data/bed/bigyri_array.bed



# ---------------------------------------------------------------------- 
# generate updated EIGENSTRAT 2.2 M sites data that include the new
# processing of Altai, Vindija and Denisova (genotyped using snpAD)

mkdir data/eigenstrat/bigyri_ho
cd data/eigenstrat/bigyri_ho

# make a copy of Qiaomei's combined Eigenstrat dataset
cp /mnt/454/Carbon_beast_QM/TY/snp/UPA_all.{snp,ind,geno} .

# make a directory for different coordinate files
less UPA_all.snp | awk -v OFS="\t" '{print $2, $4-1, $4}' > ../../bed/2.2M.bed

seq 1 22 | xargs -P 22 -I {} bash -c "bcftools view -a /mnt/scratch/steffi/D/Vcfs/mergedArchModernApes/merged_high_apes_low_sgdp1_chr{}.vcf.gz -R ../../bed/2.2M.bed -s AltaiNeandertal,Vindija33.19,Mezmais1Deam,Denisova | bcftools query -f '%CHROM\t%POS[\t%GT]\n' | sed 's/\.\/\./9/g; s/0\/0/0/g; s/0\/1/1/g; s/1\/1/2/g' | grep -v '/' > chr{}.tmp"
cat chr{1..22}.tmp > all_chr.tmp
rm chr*.tmp

R --no-save < <(echo '
library(tidyverse)
library(admixr)
all <- read_snp("UPA_all.snp")
archaics <- read_tsv("asd.tmp", col_names=c("chrom", "pos", "Altai", "Vindija", "Mez1", "Denisovan"), col_types="ciiiii")
joined <- left_join(all, archaics)
joined$Altai[is.na(joined$Altai)] <- 9
joined$Vindija[is.na(joined$Vindija)] <- 9
joined$Mez1[is.na(joined$Mez1)] <- 9
joined$Denisovan[is.na(joined$Denisovan)] <- 9

write_geno("archaics.geno", select(joined, Altai:Denisovan))
write_snp("archaics.snp", select(joined, -(Altai:Denisovan)))
write_ind("archaics.ind", tibble(name=c("new_Altai", "new_Vindija", "new_Mez1", "new_Densiovan"), sex=rep("F", 4), pop=c("new_Altai", "new_Vindija", "new_Mez1", "new_Denisovan")))
')

# generate a mergit parameter file
echo "outputformat: EIGENSTRAT
geno1: UPA_all.geno
snp1: UPA_all.snp
ind1: UPA_all.ind
geno2: archaics.geno
snp2: archaics.snp
ind2: archaics.ind
genooutfilename: UPA_merged.geno
snpoutfilename: UPA_merged.snp
indoutfilename: UPA_merged.ind" > mergeit_archaics.par

mergeit -p mergeit_archaics.par


R --no-save < <(echo '
library(tidyverse)
library(admixr)

all_data <- read_eigenstrat("UPA_all")
archaics <- read_tsv("archaics.tmp", col_names=c("chrom", "pos", "Altai", "Vindija", "Mez1", "Denisovan"), col_types="ciiiii")

joined <- left_join(all_data$snp, archaics)
joined$Altai[is.na(joined$Altai)] <- 9
joined$Vindija[is.na(joined$Vindija)] <- 9
joined$Mez1[is.na(joined$Mez1)] <- 9
joined$Denisovan[is.na(joined$Denisovan)] <- 9

write_geno("archaics.geno", select(joined, Altai:Denisovan))
write_snp("all.snp", all_data$snp)
write_ind("all.ind", bind_rows(all_data$ind, tibble(id=c("new_Altai", "new_Vindija", "new_Mez1", "new_Densiovan"), sex=rep("F", 4), label=c("new_Altai", "new_Vindija", "new_Mez1", "new_Denisovan"))))
')

paste -d '' UPA_all.geno archaics.geno > all.geno





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
