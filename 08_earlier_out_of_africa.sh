#!/bin/bash

input_dir=input_data
clean_dir=clean_data
sims_dir=simulations
tmp_dir=tmp

exome_and_sites_exon_coords=${clean_dir}/exome_and_sites_exon_coordinates.txt
exome_and_sites_recomb_map=${clean_dir}/exome_and_sites_recombination_map.txt
exonic_array_sites=${clean_dir}/admixture_array_coordinates_exonic.txt
nonexonic_array_sites=${clean_dir}/admixture_array_coordinates_nonexonic.txt


# initial Neanderthal proportions
init_nea=0.1

# number of replicates of each model
num_replicates=10

traject_dir=${sims_dir}/earlier_ooa
mkdir -p $traject_dir

# wrapper for submitting introgression simulations to SGE
run_introgression() {
    run_id=ooa_$8__h_${2}__init_nea_${5}__rep_${1}

    python3 run_introgression.py \
	    --population-file ${sims_dir}/exome_and_sites__h_${2}__seed_*.txt \
	    --exon-coordinates $3 \
	    --recomb-map $4 \
	    --admixture-rate $5 \
	    --dominance-coef $2 \
	    --exonic-sites $6 \
	    --nonexonic-sites $7 \
	    --output-prefix ${traject_dir}/${run_id} \
	    --model gravel \
	    --out-of-africa $8 \
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

for dominance_coef in 0.5 0.1 0.0; do
    for ooa in 80000 100000; do
	for rep_i in `seq 1 $num_replicates`; do
            run_introgression $rep_i $dominance_coef $exome_and_sites_exon_coords $exome_and_sites_recomb_map $init_nea $exonic_array_sites $nonexonic_array_sites $ooa &
	done
    done
done

cp ${sims_dir}/different_models/gravel__h_*_sites.txt $traject_dir
rename 's/gravel__/ooa_55000__/' ${traject_dir}/gravel__*
