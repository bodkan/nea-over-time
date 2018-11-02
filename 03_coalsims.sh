#!/usr/bin/env bash

mkdir -p data/msprime_sims

sample_times=`seq 50000 -2000 0`

r=0.1

for t in 2000 5000 10000 15000 20000; do
    scenario="geneflow_${t}"

    m=`echo "$r / ($t / 25)" | bc -l`

    # EUR -> AFR migration
    for rep in `seq 1 500`; do
        N="${scenario}_eur_to_afr_${rep}"
        qsub -V -cwd -j y -l virtual_free=10G,h_vmem=10G -N $N -o tmp/${N}.txt \
            ./code/coalsim.py \
                --time ${t} --duration ${t} --eur-to-dinka ${m} --eur-to-yoruba ${m} --dinka-to-eur 0 --nea-rate 0.03 \
                --hap-length 100_000_000 --eur-ages ${sample_times} \
                --output-file data/msprime_sims/${N}.tsv --calc-stats f4
    done


    # EUR <- AFR migration
    for rep in `seq 1 500`; do
        N="${scenario}_afr_to_eur_${rep}"
        qsub -V -cwd -j y -l virtual_free=10G,h_vmem=10G -N $N -o tmp/${N}.txt \
            ./code/coalsim.py \
                --time ${t} --duration ${t} --eur-to-dinka 0 --eur-to-yoruba 0 --dinka-to-eur ${m} --nea-rate 0.03 \
                --hap-length 100_000_000 --eur-ages ${sample_times} \
                --output-file data/msprime_sims/${N}.tsv --calc-stats f4
    done

    # EUR <-> AFR migration
    for rep in `seq 1 500`; do
        N="${scenario}_both_directions_${rep}"
        qsub -V -cwd -j y -l virtual_free=10G,h_vmem=10G -N $N -o tmp/${N}.txt \
            ./code/coalsim.py \
                --time ${t} --duration ${t} --eur-to-dinka ${m} --eur-to-yoruba ${m} --dinka-to-eur ${m} --nea-rate 0.03 \
                --hap-length 100_000_000 --eur-ages ${sample_times} \
                --output-file data/msprime_sims/${N}.tsv --calc-stats f4
    done

    # no EUR <-> AFR migration
    for rep in `seq 1 500`; do
        N="${scenario}_no_migration_${rep}"
        qsub -V -cwd -j y -l virtual_free=10G,h_vmem=10G -N $N -o tmp/${N}.txt \
            ./code/coalsim.py \
                --time 0 --duration 0 --eur-to-dinka 0 --eur-to-yoruba 0 --dinka-to-eur 0 --nea-rate 0.03 \
                --hap-length 100_000_000 --eur-ages ${sample_times} \
                --output-file data/msprime_sims/${N}.tsv --calc-stats f4
    done

done
