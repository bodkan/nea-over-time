#!/bin/bash

mkdir data/burnins

# normal burnins for different regions
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

# burnin of a merge of regions showing functional significance
region="merged"; h="0.5"
python3 code/run_mutation_accumulation.py \
    --regions data/slim_coords/${region}_regions.bed \
    --sites data/slim_coords/${region}_all_sites.bed \
    --recomb-map data/slim_coords/${region}_recomb_map.bed \
    --mut-rate 1e-8 \
    --dominance-coef $h \
    --output data/burnins/${region}_h_${h}.txt &

# testing the effect of different Neanderthal population sizes
region="exon"; h=0.5
for Ne in 100 500 10000; do
    python3 code/run_mutation_accumulation.py \
        --regions data/slim_coords/${region}_regions.bed \
        --sites data/slim_coords/${region}_all_sites.bed \
        --recomb-map data/slim_coords/${region}_recomb_map.bed \
        --nea-size $Ne \
        --mut-rate 1e-8 \
        --dominance-coef $h \
        --output data/burnins/nea_Ne_${Ne}_${region}_h_${h}.txt &
done
cp data/burnins/exon_h_0.5.txt data/burnins/nea_Ne_1000_exon_h_0.5.txt

# uniform recombination rate test
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

# "regression test" - uniform recombination rate and lower mutation rate (7e-9)
region="exon"; h=0.5
python3 code/run_mutation_accumulation.py \
    --regions data/slim_coords/${region}_regions.bed \
    --sites data/slim_coords/${region}_all_sites.bed \
    --recomb-map data/slim_coords/uniform_${region}_recomb_map.bed \
    --mut-rate 7e-9 \
    --dominance-coef $h \
    --output data/burnins/nonsyn_uniform_${region}_h_${h}.txt &
