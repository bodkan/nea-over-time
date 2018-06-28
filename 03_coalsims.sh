#!/usr/bin/env bash

mkdir data/coalsims

t=20_000
m=0.0001

sample_times=`seq 0 1000 50000`

# no EUR <-> AFR migration
for rep in `seq 1 100`; do
    N="no_migration_${rep}"
    qsub -V -cwd -j y -l virtual_free=20G,h_vmem=20G -N $N -o tmp/${N}.txt \
        ./code/coalsim.py \
            --time 0 --eur-to-afr 0 --afr-to-eur 0 --nea-rate 0.03 \
            --hap-length 500_000_000 --eur-ages ${sample_times} \
            --output-file data/coalsims/${N}.tsv --calc-stats
done

# EUR -> AFR migration
for rep in `seq 1 100`; do
    N="eur_to_afr_${rep}"
    qsub -V -cwd -j y -l virtual_free=20G,h_vmem=20G -N $N -o tmp/${N}.txt \
        ./code/coalsim.py \
            --time ${t} --eur-to-afr ${m} --afr-to-eur 0 --nea-rate 0.03 \
            --hap-length 500_000_000 --eur-ages ${sample_times} \
            --output-file data/coalsims/${N}.tsv --calc-stats
done


# EUR <- AFR migration
for rep in `seq 1 100`; do
    N="afr_to_eur_${rep}"
    qsub -V -cwd -j y -l virtual_free=20G,h_vmem=20G -N $N -o tmp/${N}.txt \
        ./code/coalsim.py \
            --time ${t} --eur-to-afr 0 --afr-to-eur ${m} --nea-rate 0.03 \
            --hap-length 500_000_000 --eur-ages ${sample_times} \
            --output-file data/coalsims/${N}.tsv --calc-stats
done

# EUR <-> AFR migration
for rep in `seq 1 100`; do
    N="both_directions_${rep}"
    qsub -V -cwd -j y -l virtual_free=20G,h_vmem=20G -N $N -o tmp/${N}.txt \
        ./code/coalsim.py \
            --time ${t} --eur-to-afr ${m} --afr-to-eur ${m} --nea-rate 0.03 \
            --hap-length 500_000_000 --eur-ages ${sample_times} \
            --output-file data/coalsims/${N}.tsv --calc-stats
done


#
# SNP table dumps for investigating BABA/ABBA changes
#
t=20_000
m=0.0001
sample_times=`seq 0 1000 50000`

./code/coalsim.py \
    --time ${t} --eur-to-afr 0 --afr-to-eur 0 --nea-rate 0.03 \
    --hap-length 500_000_000 --eur-ages ${sample_times} \
    --output-file data/coalsims/snps_no_migration.tsv --dump-snps &
./code/coalsim.py \
    --time ${t} --eur-to-afr ${m} --afr-to-eur 0 --nea-rate 0.03 \
    --hap-length 500_000_000 --eur-ages ${sample_times} \
    --output-file data/coalsims/snps_eur_to_afr.tsv --dump-snps &
./code/coalsim.py \
    --time ${t} --eur-to-afr 0 --afr-to-eur ${m} --nea-rate 0.03 \
    --hap-length 500_000_000 --eur-ages ${sample_times} \
    --output-file data/coalsims/snps_afr_to_eur.tsv --dump-snps &
./code/coalsim.py \
    --time ${t} --eur-to-afr ${m} --afr-to-eur ${m} --nea-rate 0.03 \
    --hap-length 500_000_000 --eur-ages ${sample_times} \
    --output-file data/coalsims/snps_both_directions.tsv --dump-snps &
./code/coalsim.py \
    --time ${t} --eur-to-afr ${m} --afr-to-eur 0 --nea-rate 0 \
    --hap-length 500_000_000 --eur-ages ${sample_times} \
    --output-file data/coalsims/snps_eur_to_afr_no_nea.tsv --dump-snps &

