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

traject_dir=${sims_dir}/bump_up_coding_s
mkdir -p $traject_dir

mkdir ${tmp_dir}/sge

# wrapper for submitting introgression simulations to SGE
run_introgression() {
    run_id=modify_by_${4}__fraction_${5}__h_${2}__init_nea_${3}__rep_${1}

    qsub -V -cwd -j y -S /bin/bash \
	 -l virtual_free=10G,h_vmem=10G \
	 -o ${tmp_dir}/sge/sge__${run_id}.out \
	 -N ${run_id} \
	 run_introgression.py \
	    --population-file ${sims_dir}/exome_and_sites__h_${2}__seed_*.txt \
	    --exon-coordinates $exome_and_sites_exon_coords \
	    --recomb-map $exome_and_sites_recomb_map \
	    --admixture-rate $3 \
	    --dominance-coef $2 \
	    --exonic-sites $exonic_array_sites \
	    --nonexonic-sites $nonexonic_array_sites \
	    --output-prefix ${traject_dir}/${run_id} \
	    --model constant \
	    --coding-modifier $4 \
	    --coding-fraction $5 \
	    --terminate-after 500 \
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

dominance_coef=0.5
for modifier in 1.5 2.5 5.0 10.0 50.0; do
    for fraction in 1.0 0.75 0.5 0.25 0.1; do
	for rep_i in `seq 6 $num_replicates`; do
            run_introgression $rep_i $dominance_coef $init_nea $modifier $fraction &
	done
    done
done

cp ${sims_dir}/different_models/gravel__h_*_sites.txt $traject_dir
rename 's/gravel__/modify_by_unmodified__fraction_unmodified__/' ${traject_dir}/gravel__*
