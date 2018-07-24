#!/usr/bin/env bash

mkdir data/msprime_sims




t=20_000
m=0.0001

sample_times=`seq 50000 -1000 0`

# no EUR <-> AFR migration
for rep in `seq 1 10`; do
    N="snps_no_migration_${rep}"
    qsub -V -cwd -j y -l virtual_free=30G,h_vmem=30G -N $N -o tmp/${N}.txt \
        ./code/coalsim.py \
            --time 0 --eur-to-afr 0 --afr-to-eur 0 --nea-rate 0.03 \
            --hap-length 500_000_000 --eur-ages ${sample_times} \
            --output-file data/msprime_sims/${N}.tsv --dump-snp
done

# EUR -> AFR migration
for rep in `seq 1 10`; do
    N="snps_eur_to_afr_${rep}"
    qsub -V -cwd -j y -l virtual_free=30G,h_vmem=30G -N $N -o tmp/${N}.txt \
        ./code/coalsim.py \
            --time ${t} --eur-to-afr ${m} --afr-to-eur 0 --nea-rate 0.03 \
            --hap-length 500_000_000 --eur-ages ${sample_times} \
            --output-file data/msprime_sims/${N}.tsv --dump-snp
done


# EUR <- AFR migration
for rep in `seq 1 10`; do
    N="snps_afr_to_eur_${rep}"
    qsub -V -cwd -j y -l virtual_free=30G,h_vmem=30G -N $N -o tmp/${N}.txt \
        ./code/coalsim.py \
            --time ${t} --eur-to-afr 0 --afr-to-eur ${m} --nea-rate 0.03 \
            --hap-length 500_000_000 --eur-ages ${sample_times} \
            --output-file data/msprime_sims/${N}.tsv --dump-snp
done

# EUR <-> AFR migration
for rep in `seq 1 10`; do
    N="snps_both_directions_${rep}"
    qsub -V -cwd -j y -l virtual_free=30G,h_vmem=30G -N $N -o tmp/${N}.txt \
        ./code/coalsim.py \
            --time ${t} --eur-to-afr ${m} --afr-to-eur ${m} --nea-rate 0.03 \
            --hap-length 500_000_000 --eur-ages ${sample_times} \
            --output-file data/msprime_sims/${N}.tsv --dump-snp
done









# no EUR <-> AFR migration
for rep in `seq 1 20`; do
    N="symmetry_no_migration_${rep}"
    qsub -V -cwd -j y -l virtual_free=30G,h_vmem=30G -N $N -o tmp/${N}.txt \
        ./code/coalsim.py \
            --time 0 --eur-to-afr 0 --afr-to-eur 0 --nea-rate 0.03 \
            --hap-length 500_000_000 --eur-ages ${sample_times} \
            --output-file data/msprime_sims/${N}.tsv --calc-stats
done

# EUR -> AFR migration
for rep in `seq 1 20`; do
    N="symmetry_eur_to_afr_${rep}"
    qsub -V -cwd -j y -l virtual_free=30G,h_vmem=30G -N $N -o tmp/${N}.txt \
        ./code/coalsim.py \
            --time ${t} --eur-to-afr ${m} --afr-to-eur 0 --nea-rate 0.03 \
            --hap-length 500_000_000 --eur-ages ${sample_times} \
            --output-file data/msprime_sims/${N}.tsv --calc-stats
done


# EUR <- AFR migration
for rep in `seq 1 20`; do
    N="symmetry_afr_to_eur_${rep}"
    qsub -V -cwd -j y -l virtual_free=30G,h_vmem=30G -N $N -o tmp/${N}.txt \
        ./code/coalsim.py \
            --time ${t} --eur-to-afr 0 --afr-to-eur ${m} --nea-rate 0.03 \
            --hap-length 500_000_000 --eur-ages ${sample_times} \
            --output-file data/msprime_sims/${N}.tsv --calc-stats
done

# EUR <-> AFR migration
for rep in `seq 1 20`; do
    N="symmetry_both_directions_${rep}"
    qsub -V -cwd -j y -l virtual_free=30G,h_vmem=30G -N $N -o tmp/${N}.txt \
        ./code/coalsim.py \
            --time ${t} --eur-to-afr ${m} --afr-to-eur ${m} --nea-rate 0.03 \
            --hap-length 500_000_000 --eur-ages ${sample_times} \
            --output-file data/msprime_sims/${N}.tsv --calc-stats
done
