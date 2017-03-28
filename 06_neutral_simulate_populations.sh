#!/bin/bash

input_dir=clean_data
sims_dir=simulations
tmp_dir=tmp


#
# Simulations of AMH and Neanderthal demography until the introgression event.
# Exome and sites from the archaic admixture array.
#

exome_and_sites_exon_coords=${input_dir}/exome_and_sites_exon_coordinates.txt
exome_and_sites_recomb_map=${input_dir}/exome_and_sites_recombination_map.txt
exonic_array_sites=${input_dir}/admixture_array_coordinates_exonic.txt
nonexonic_array_sites=${input_dir}/admixture_array_coordinates_nonexonic.txt


h=1.0
python3 run_mutation_accumulation.py --exon-coordinates $exome_and_sites_exon_coords \
                                     --recomb-map $exome_and_sites_recomb_map \
				     --mut-rate 0 \
                                     --exonic-sites $exonic_array_sites \
                                     --nonexonic-sites $nonexonic_array_sites \
                                     --dominance-coef $h \
                                     --output-prefix ${sims_dir}/exome_and_sites__neutral__ \
                                     --anc-size 1 --hum-nea-split 50 --out-of-africa 25 &
