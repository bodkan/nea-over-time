#!/bin/bash

input_dir=input_data
clean_dir=clean_data
sims_dir=simulations

exome_and_sites_exon_coords=${clean_dir}/exome_and_sites_exon_coordinates.txt
exome_and_sites_recomb_map=${clean_dir}/exome_and_sites_recombination_map.txt
exonic_array_sites=${clean_dir}/admixture_array_coordinates_exonic.txt
nonexonic_array_sites=${clean_dir}/admixture_array_coordinates_nonexonic.txt

./run_introgression.py \
    	    --population-file ${sims_dir}/exome_and_sites__h_0.5__seed_*.txt \
    	    --exon-coordinates $exome_and_sites_exon_coords \
    	    --recomb-map $exome_and_sites_recomb_map \
    	    --admixture-rate $1 \
    	    --dominance-coef 0.5 \
    	    --exonic-sites $exonic_array_sites \
    	    --nonexonic-sites $nonexonic_array_sites \
	    --founder-size $2 \
    	    --output-prefix $3 \
    	    --model constant \
    	    --save-trajectories \
