#!/bin/bash

# ----------------------------------------------------------------------

mkdir data/burnins

for region in protein_coding promoter utr3 tf_binding_site; do
for h in 0.0 0.5 1.0; do
    python3 code/run_mutation_accumulation.py \
        --regions data/slim_coords/${region}_regions.bed \
        --sites data/slim_coords/${region}_sites.bed \
        --recomb-map data/slim_coords/${region}_recomb_map.bed \
        --mut-rate 1e-8 \
        --dominance-coef $h \
        --output data/burnins/${region}_h_${h}.txt &
done
done

# ----------------------------------------------------------------------

mkdir data/slim_coords/subset

# generate chr1-5 subsets of the coordinate files 
for f in data/slim_coords/*.bed; do
    head -n 1 $f > data/slim_coords/subset/`basename ${f}`
    for chr in `seq 1 5`; do
        grep -w "^chr${chr}" $f
    done >> data/slim_coords/subset/`basename ${f}`
done

for region in protein_coding promoter utr3 tf_binding_site; do
for h in `seq 0.0 0.1 1.0`; do
    python3 code/run_mutation_accumulation.py \
        --regions data/slim_coords/subset/${region}_regions.bed \
        --sites data/slim_coords/subset/${region}_sites.bed \
        --recomb-map data/slim_coords/subset/${region}_recomb_map.bed \
        --mut-rate 1e-8 \
        --dominance-coef $h \
        --output data/burnins/subset_${region}_h_${h}.txt &
done
done

