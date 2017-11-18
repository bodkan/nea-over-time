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

mkdir data/burnins/toy_genome

for h in `seq 0.0 0.1 1.0`; do
    python3 code/run_mutation_accumulation.py \
        --regions data/slim_coords/toy_genome/independent_haps_regions.bed \
        --sites data/slim_coords/toy_genome/independent_haps_sites.bed \
        --recomb-map data/slim_coords/toy_genome/independent_haps_recomb_map.bed \
        --mut-rate 1e-8 \
        --dominance-coef $h \
        --output data/burnins/toy_genome/independent_haps_h_${h}.txt &
done
python3 code/run_mutation_accumulation.py \
    --regions data/slim_coords/toy_genome/independent_haps_regions.bed \
    --sites data/slim_coords/toy_genome/independent_haps_sites.bed \
    --recomb-map data/slim_coords/toy_genome/independent_haps_recomb_map.bed \
    --mut-rate 0 \
    --dominance-coef 0.0 \
    --output data/burnins/toy_genome/independent_haps_neutral.txt &


# ----------------------------------------------------------------------

for h in `seq 0.0 0.1 1.0`; do
    python3 code/run_mutation_accumulation.py \
        --regions data/slim_coords/toy_genome/linked_haps_regions.bed \
        --sites data/slim_coords/toy_genome/linked_haps_sites.bed \
        --recomb-map data/slim_coords/toy_genome/linked_haps_recomb_map.bed \
        --mut-rate 1e-8 \
        --dominance-coef $h \
        --output data/burnins/toy_genome/linked_haps_h_${h}.txt &
done
python3 code/run_mutation_accumulation.py \
    --regions data/slim_coords/toy_genome/linked_haps_regions.bed \
    --sites data/slim_coords/toy_genome/linked_haps_sites.bed \
    --recomb-map data/slim_coords/toy_genome/linked_haps_recomb_map.bed \
    --mut-rate 0 \
    --dominance-coef 0.0 \
    --output data/burnins/toy_genome/linked_haps_neutral.txt &
