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
# introgression simulations
#
traject_dir=${sims_dir}/fixed_nea_s
mkdir -p $traject_dir


# initial Neanderthal proportions
init_nea=0.1

# number of replicates of each model
num_replicates=10

# wrapper for submitting introgression simulations to SGE
run_introgression() {
    run_id=${5}__s_${4}__h__${2}__init_nea_${3}__rep_${1}

    ./run_introgression.py \
    	    --population-file ${sims_dir}/exome_and_sites__h_${2}__seed_*.txt \
    	    --exon-coordinates $exome_and_sites_exon_coords \
    	    --recomb-map $exome_and_sites_recomb_map \
    	    --admixture-rate $3 \
    	    --dominance-coef $2 \
    	    --fixed-nea-s=-$4 \
    	    --exonic-sites $exonic_array_sites \
    	    --nonexonic-sites $nonexonic_array_sites \
    	    --output-prefix ${traject_dir}/${run_id} \
    	    --model $5 \
    	    --save-trajectories
}


# 
# Simulations of Neanderthal ancestry trajectories in the European
# population using previously simulated populations with accumulated
# deleterious variants.
#
# Neutral Neanderthal markers are exonic and non-exonic sites from
# the archaic admixture array.
#
init_nea=0.1


h=0.5

for s in 4.7e-4 4.7e-3 1e-3; do
for model in constant gravel; do
    for rep_i in `seq 1 $num_replicates`; do
	run_introgression $rep_i $h $init_nea $s $model &
    done
done
done

cp ${sims_dir}/different_models/{constant,gravel}__h_0.5_*sites.txt $traject_dir
