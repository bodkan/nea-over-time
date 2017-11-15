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

traject_dir=${sims_dir}/fixed_s
mkdir -p $traject_dir

mkdir -p ${tmp_dir}/sge

# wrapper for submitting introgression simulations to SGE
run_introgression() {
    run_id=fixed_s_${3}__h_0.5__init_nea_${2}__rep_${1}

    # echo qsub -V -cwd -j y -S /bin/bash \
    # 	 -l virtual_free=10G,h_vmem=10G \
    # 	 -o ${tmp_dir}/sge/sge__${run_id}.out \
    # 	 -N ${run_id} \
	 ./run_introgression.py \
	    --population-file ${sims_dir}/exome_and_sites__h_0.5__seed_*.txt \
	    --exon-coordinates $exome_and_sites_exon_coords \
	    --recomb-map $exome_and_sites_recomb_map \
	    --admixture-rate $2 \
	    --dominance-coef 0.5 \
	    --exonic-sites $exonic_array_sites \
	    --nonexonic-sites $nonexonic_array_sites \
	    --output-prefix ${traject_dir}/${run_id} \
	    --model constant \
	    --fix-s=-$3 \
	    --modify-fraction 1.0 \
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

# dominance_coef=0.5
# fix_s=4.7e-4
# for rep_i in `seq 1 $num_replicates`; do
#         run_introgression $rep_i $init_nea $fix_s &
# done

cp ${sims_dir}/different_models/constant__h_0.5__*_sites.txt $traject_dir
rename 's/constant__/unmodified__/' ${traject_dir}/constant__*
