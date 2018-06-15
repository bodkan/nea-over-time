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

# Neanderthal and Denisovan deserts
for region in exon protein_coding; do
for chrom in chr1 chr7; do
for rep in `seq 1 100`; do
  N="${chrom}_${region}_${rep}"
  qsub -V -cwd -j y -l virtual_free=20G,h_vmem=20G -N $N -o tmp/${N}.txt \
  ./code/run_mutation_accumulation.py \
      --chrom $chrom \
      --regions data/slim_coords/${region}_unif_regions.bed \
      --sites data/slim_coords/${region}_unif_all_sites.bed \
      --recomb-map data/slim_coords/${region}_unif_recomb_map.bed \
      --mut-rate 7e-9 \
      --dominance-coef 0.5 \
      --nea-den-split 400000 \
      --output data/burnins/${chrom}_${region}_${rep}.txt
done
done
done

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
