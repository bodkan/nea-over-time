#!/bin/bash

# Neanderthal vs Denisovan deserts
seq 1 50 | xargs -I {} -c "
  python3 code/run_mutation_accumulation.py \
      --regions data/slim_coords/${region}_regions.bed \
      --sites data/slim_coords/${region}_all_sites.bed \
      --recomb-map data/slim_coords/${region}_recomb_map.bed \
      --mut-rate 1e-8 \
      --dominance-coef 0.5 \
      --nea-den-split 400000 \
      --output data/burnins/${region}_h_${h}.txt &
"