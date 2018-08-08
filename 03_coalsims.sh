#!/usr/bin/env bash

mkdir data/msprime_sims

sample_times=`seq 50000 -2000 0`

r=0.05

for asc in all dinka yoruba eur eas; do 

for t in 5000 10000 15000; do

    m=`echo "$r / ($t / 25)" | bc -l`

    scenario="ascertainment_${asc}_geneflow_both_afr_${t}"

    # no EUR <-> AFR migration
    for rep in `seq 1 50`; do
        N="${scenario}_no_migration_${rep}"
        qsub -V -cwd -j y -l virtual_free=30G,h_vmem=30G -N $N -o tmp/${N}.txt \
            ./code/coalsim.py \
                --time 0 --eur-to-dinka 0 --eur-to-yoruba 0 --dinka-to-eur 0 --nea-rate 0.03 \
                --hap-length 500_000_000 --eur-ages ${sample_times} \
                --output-file data/msprime_sims/${N}.tsv --calc-stats --ascertainment $asc
    done
    
    # EUR -> AFR migration
    for rep in `seq 1 50`; do
        N="${scenario}_eur_to_afr_${rep}"
        qsub -V -cwd -j y -l virtual_free=30G,h_vmem=30G -N $N -o tmp/${N}.txt \
            ./code/coalsim.py \
                --time ${t} --eur-to-dinka ${m} --eur-to-yoruba ${m} --dinka-to-eur 0 --nea-rate 0.03 \
                --hap-length 500_000_000 --eur-ages ${sample_times} \
                --output-file data/msprime_sims/${N}.tsv --calc-stats --ascertainment $asc
    done
    
    
    # EUR <- AFR migration
    for rep in `seq 1 50`; do
        N="${scenario}_afr_to_eur_${rep}"
        qsub -V -cwd -j y -l virtual_free=30G,h_vmem=30G -N $N -o tmp/${N}.txt \
            ./code/coalsim.py \
                --time ${t} --eur-to-dinka 0 --eur-to-yoruba 0 --dinka-to-eur ${m} --nea-rate 0.03 \
                --hap-length 500_000_000 --eur-ages ${sample_times} \
                --output-file data/msprime_sims/${N}.tsv --calc-stats --ascertainment $asc
    done
    
    # EUR <-> AFR migration
    for rep in `seq 1 50`; do
        N="${scenario}_both_directions_${rep}"
        qsub -V -cwd -j y -l virtual_free=30G,h_vmem=30G -N $N -o tmp/${N}.txt \
            ./code/coalsim.py \
                --time ${t} --eur-to-dinka ${m} --eur-to-yoruba ${m} --dinka-to-eur ${m} --nea-rate 0.03 \
                --hap-length 500_000_000 --eur-ages ${sample_times} \
                --output-file data/msprime_sims/${N}.tsv --calc-stats --ascertainment $asc
    done

done
done
