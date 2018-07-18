#!/usr/bin/env bash

mkdir data/msprime_sims




t=20_000
m=0.0001

sample_times=`seq 50000 -1000 0`

# no EUR <-> AFR migration
for rep in `seq 1 100`; do
    N="no_migration_${rep}"
    qsub -V -cwd -j y -l virtual_free=30G,h_vmem=30G -N $N -o tmp/${N}.txt \
        ./code/coalsim.py \
            --time 0 --eur-to-afr1 0 --eur-to-afr2 0 --afr-to-eur 0 --nea-rate 0.03 \
            --hap-length 500_000_000 --eur-ages ${sample_times} \
            --output-file data/msprime_sims/symmetry_${N}.tsv --calc-stats
done

# EUR -> AFR migration
for rep in `seq 1 100`; do
    N="eur_to_afr_${rep}"
    qsub -V -cwd -j y -l virtual_free=30G,h_vmem=30G -N $N -o tmp/${N}.txt \
        ./code/coalsim.py \
            --time ${t} --eur-to-afr1 ${m} --eur-to-afr2 ${m} --afr-to-eur 0 --nea-rate 0.03 \
            --hap-length 500_000_000 --eur-ages ${sample_times} \
            --output-file data/msprime_sims/symmetry_${N}.tsv --calc-stats
done


# EUR <- AFR migration
for rep in `seq 1 100`; do
    N="afr_to_eur_${rep}"
    qsub -V -cwd -j y -l virtual_free=30G,h_vmem=30G -N $N -o tmp/${N}.txt \
        ./code/coalsim.py \
            --time ${t} --eur-to-afr1 0 --eur-to-afr2 0 --afr-to-eur ${m} --nea-rate 0.03 \
            --hap-length 500_000_000 --eur-ages ${sample_times} \
            --output-file data/msprime_sims/symmetry_${N}.tsv --calc-stats
done

# EUR <-> AFR migration
for rep in `seq 1 100`; do
    N="both_directions_${rep}"
    qsub -V -cwd -j y -l virtual_free=30G,h_vmem=30G -N $N -o tmp/${N}.txt \
        ./code/coalsim.py \
            --time ${t} --eur-to-afr1 ${m} --eur-to-afr2 ${m} --afr-to-eur ${m} --nea-rate 0.03 \
            --hap-length 500_000_000 --eur-ages ${sample_times} \
            --output-file data/msprime_sims/symmetry_${N}.tsv --calc-stats
done








t=10_000
m1=0.0005
m2=0.0001
m=0.0001

sample_times=`seq 50000 -1000 0`

# no EUR <-> AFR migration
for rep in `seq 1 50`; do
    N="asymmetry_no_migration_${rep}"
    qsub -V -cwd -j y -l virtual_free=30G,h_vmem=30G -N $N -o tmp/${N}.txt \
        ./code/coalsim.py \
            --time 0 --eur-to-afr1 0 --eur-to-afr2 0 --afr-to-eur 0 --nea-rate 0.03 \
            --hap-length 500_000_000 --eur-ages ${sample_times} \
            --output-file data/msprime_sims/symmetry_${N}.tsv --calc-stats
done

# EUR -> AFR migration
for rep in `seq 1 50`; do
    N="asymmetry_eur_to_afr_${rep}"
    qsub -V -cwd -j y -l virtual_free=30G,h_vmem=30G -N $N -o tmp/${N}.txt \
        ./code/coalsim.py \
            --time ${t} --eur-to-afr1 ${m1} --eur-to-afr2 ${m2} --afr-to-eur 0 --nea-rate 0.03 \
            --hap-length 500_000_000 --eur-ages ${sample_times} \
            --output-file data/msprime_sims/symmetry_${N}.tsv --calc-stats
done


# EUR <- AFR migration
for rep in `seq 1 50`; do
    N="asymmetry_afr_to_eur_${rep}"
    qsub -V -cwd -j y -l virtual_free=30G,h_vmem=30G -N $N -o tmp/${N}.txt \
        ./code/coalsim.py \
            --time ${t} --eur-to-afr1 0 --eur-to-afr2 0 --afr-to-eur ${m} --nea-rate 0.03 \
            --hap-length 500_000_000 --eur-ages ${sample_times} \
            --output-file data/msprime_sims/symmetry_${N}.tsv --calc-stats
done

# EUR <-> AFR migration
for rep in `seq 1 50`; do
    N="asymmetry_both_directions_${rep}"
    qsub -V -cwd -j y -l virtual_free=30G,h_vmem=30G -N $N -o tmp/${N}.txt \
        ./code/coalsim.py \
            --time ${t} --eur-to-afr1 ${m1} --eur-to-afr2 ${m2} --afr-to-eur ${m} --nea-rate 0.03 \
            --hap-length 500_000_000 --eur-ages ${sample_times} \
            --output-file data/msprime_sims/symmetry_${N}.tsv --calc-stats
done
