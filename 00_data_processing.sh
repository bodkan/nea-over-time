#!/usr/bin/env bash

# This script processes the published data from the Ice Age paper,
# SGDP VCF files and also VCF files with Altai, Vindija and Denisova genotypes.
# Finally, it processes the CADD functional annotation data of SNPs
# from the archaic admixture array.

# download Fu et al. data
mkdir -p raw_data
cd raw_data
wget https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/FuQ.zip
unzip FuQ.zip
cd ../

mkdir -p clean_data

# process Fu et al. data
./R/process_eigenstrat.R raw_data/archaic.geno raw_data/archaic.snp raw_data/archaic.ind clean_data/ice_age.tsv
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


# annotate the SNPs using CADD
# the SNPs were annotated "manually" by uploading the ice_age.tsv pseudo-VCF file
# to the CADD online server and results were downloaded to raw_data/
cd raw_data
gunzip v1.3_anno*
sed 's/#//;2q;d' v1.3_anno*.tsv > annotations.tsv_unsorted
for file in v1.3_anno*.tsv; do
    tail -n+3 $file >> annotations.tsv_unsorted
done
sort -k1,1n -k2,2n annotations.tsv_unsorted > ../clean_data/annotations.tsv
cd ../
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
curl ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz -o raw_data/gtf.gz
python3 exon_annotations.py --input-file clean_data/ice_age.tsv --gtf-file raw_data/gtf.gz --window-sizes 10000 25000 50000 100000 --output-file $annotation_dir/exon_distance_and_density.bed

# generate whole-genome BED files with different annotations for
# later window-based analysis
# extract chrom, pos, priPhCons, mamPhCons, priPhyloP, mamPhyloP, bStatistic
# column numbers from here: http://cadd.gs.washington.edu/static/ReleaseNotes_CADD_v1.3.pdf
for chr in 1 {10..19} 2 20 21 22 {3..9}; do
    tabix /mnt/expressions/cadd/whole_genome_SNVs_inclAnno.tsv.gz ${chr} \
        | awk -vOFS="\t" '{ $1=$1"\t"$2-1; print $1, $2, $19, $20, $22, $23, $29 }' \
        | uniq
done | gzip > ${annotation_dir}/genome_wide_annotations.sorted.bed.gz


# calculate conservation statistics for different tracks and different window sizes
window_average() {
    track=$1      # track name
    column=$2     # track's column in the "annot_file"
    annot_file=$4 # BED file with the conservation track
    snp_file=$5   # table with SNP coordinates (needs to have 2 columns: chrom and pos)
    flank_size=`echo $3 / 2 | bc`

    zcat $annot_file \
	| tr '\t' ' ' \
	| cut -d' ' -f1-3,${column} \
	| awk -vOFS='\t' '{ print $1, $2, $3, $4, $4 }' \
	| grep -v "NA" \
	| bedmap --ec --delim '\t' --range $flank_size --echo --count --mean --kth .05 --kth .95 --median \
		 <(tail -n+2 $snp_file | awk -vOFS='\t' '{print $1, $2-1, $2}' | sort-bed -) - \
        > $annotation_dir/${track}__window_${3}bp.bed
}

snp_file=clean_data/ice_age.tsv
annot_file=raw_data/genome_wide_annotations.sorted.bed.gz

for window_size in 10000 25000 50000 100000; do
    window_average priPhCons      4 $window_size $annot_file $snp_file &
    window_average priPhyloP      6 $window_size $annot_file $snp_file &
    window_average bval           8 $window_size $annot_file $snp_file &
    window_average phyloP_nohuman 4 $window_size <(for chr in 1 {10..19} 2 20 21 22 {3..9}; do cat /mnt/expressions/benjamin_vernot/phyloP_no_human/Compara.36_eutherian_mammals_EPO_LOW_COVERAGE.chr${chr}_*; done | sed 's/^chr//' | gzip) $snp_file &
done
