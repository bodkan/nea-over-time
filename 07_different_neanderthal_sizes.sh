#!/bin/bash

input_dir=input_data
clean_dir=clean_data
sims_dir=simulations
tmp_dir=tmp

exome_and_sites_exon_coords=${clean_dir}/exome_and_sites_exon_coordinates.txt
exome_and_sites_recomb_map=${clean_dir}/exome_and_sites_recombination_map.txt
exonic_array_sites=${clean_dir}/admixture_array_coordinates_exonic.txt
nonexonic_array_sites=${clean_dir}/admixture_array_coordinates_nonexonic.txt


#
# Simulations of AMH and Neanderthal demography until the introgression event.
# Exome and sites from the archaic admixture array.


# h=0.5
# python3 run_mutation_accumulation.py --exon-coordinates $exome_and_sites_exon_coords \
# 	                             --recomb-map $exome_and_sites_recomb_map \
#                                      --exonic-sites $exonic_array_sites \
#                                      --nonexonic-sites $nonexonic_array_sites \
# 				     --dominance-coef $h \
# 				     --nea-size 100 \
# 				     --output-prefix ${sims_dir}/exome_and_sites__h_${h}__Ne_100__ &

# h=0.5
# python3 run_mutation_accumulation.py --exon-coordinates $exome_and_sites_exon_coords \
# 	                             --recomb-map $exome_and_sites_recomb_map \
#                                      --exonic-sites $exonic_array_sites \
#                                      --nonexonic-sites $nonexonic_array_sites \
# 				     --dominance-coef $h \
# 				     --nea-size 500 \
# 				     --output-prefix ${sims_dir}/exome_and_sites__h_${h}__Ne_500__ &

# h=0.5
# python3 run_mutation_accumulation.py --exon-coordinates $exome_and_sites_exon_coords \
# 	                             --recomb-map $exome_and_sites_recomb_map \
#                                      --exonic-sites $exonic_array_sites \
#                                      --nonexonic-sites $nonexonic_array_sites \
# 				     --dominance-coef $h \
# 				     --nea-size 10000 \
# 				     --output-prefix ${sims_dir}/exome_and_sites__h_${h}__Ne_10000__ &



# initial Neanderthal proportions
init_nea=0.1

# number of replicates of each model
num_replicates=10

traject_dir=${sims_dir}/nea_pop_sizes

# wrapper for submitting introgression simulations to SGE
run_introgression() {
    run_id=Ne_$8__h_${2}__init_nea_${5}__rep_${1}

    python3 run_introgression.py \
	    --population-file ${sims_dir}/exome_and_sites__h_${2}__Ne_${8}__seed_*.txt \
	    --exon-coordinates $3 \
	    --recomb-map $4 \
	    --admixture-rate $5 \
	    --dominance-coef $2 \
	    --exonic-sites $6 \
	    --nonexonic-sites $7 \
	    --output-prefix ${traject_dir}/${run_id} \
	    --model constant \
	    --save-trajectories \
	    ${9}
}


#
# Simulations of Neanderthal ancestry trajectories in the European
# population using previously simulated populations with accumulated
# deleterious variants.
#
# Neutral Neanderthal markers are exonic and non-exonic sites from
# the archaic admixture array.
#

dominance_coef=0.5
for nea_Ne in 100 500 1000 10000; do
    run_introgression 1 $dominance_coef $exome_and_sites_exon_coords $exome_and_sites_recomb_map $init_nea $exonic_array_sites $nonexonic_array_sites $model $nea_Ne --save-mutations &
    for rep_i in `seq 2 $num_replicates`; do
        run_introgression $rep_i $dominance_coef $exome_and_sites_exon_coords $exome_and_sites_recomb_map $init_nea $exonic_array_sites $nonexonic_array_sites $model $nea_Ne &
    done
done
