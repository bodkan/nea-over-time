#!/bin/bash

mkdir data/burnins

for region in exon protein_coding promoter utr3 tf_binding_site; do
for h in 0.0 0.5 1.0; do
    python3 code/run_mutation_accumulation.py \
        --regions data/slim_coords/${region}_regions.bed \
        --sites data/slim_coords/${region}_all_sites.bed \
        --recomb-map data/slim_coords/${region}_recomb_map.bed \
        --mut-rate 1e-8 \
        --dominance-coef $h \
        --output data/burnins/${region}_h_${h}.txt &
done
done

for region in exon protein_coding promoter utr3 tf_binding_site; do
for h in 0.5; do
    python3 code/run_mutation_accumulation.py \
        --regions data/slim_coords/${region}_regions.bed \
        --sites data/slim_coords/${region}_all_sites.bed \
        --recomb-map data/slim_coords/uniform_${region}_recomb_map.bed \
        --mut-rate 1e-8 \
        --dominance-coef $h \
        --output data/burnins/uniform_${region}_h_${h}.txt &
done
done
