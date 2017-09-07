#!/bin/bash

input_dir=input_data
clean_dir=clean_data
sims_dir=simulations
tmp_dir=tmp

exome_and_sites_exon_coords=${clean_dir}/exome_and_sites_exon_coordinates.txt
exome_and_sites_recomb_map=${clean_dir}/exome_and_sites_recombination_map.txt
exonic_array_sites=${clean_dir}/admixture_array_coordinates_exonic.txt
nonexonic_array_sites=${clean_dir}/admixture_array_coordinates_nonexonic.txt

traject_dir=${sims_dir}/delayed_fixed_s
mkdir -p $traject_dir

# wrapper for submitting introgression simulations to SGE
run_introgression() {
    run_id=${5}__mut_count_${4}__fixed_s_${3}__init_nea_${2}__rep_${1}

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
            --founder-size 2500 \
	    --fix-s=-$3 \
            --modify-at 300 \
	    --modify-count $4 \
            --modify-what $5 \
	    --save-trajectories
}

# mut_type=m0
# init_nea=0.1
# for mut_count in {100,500,1000,2000}; do
#     for fix_s in {0.0001,0.001,0.005,0.01}; do
#         for rep_i in `seq 1 $num_replicates`; do
#             run_introgression $rep_i $init_nea $fix_s $mut_count $mut_type &
#         done
#     done
# done


mut_type=m2
init_nea=0.1
rep_i=$1
for mut_count in {100,500,1000,2000,3000,10000}; do
    for fix_s in {0.0001,0.001,0.005,0.01}; do
        run_introgression $rep_i $init_nea $fix_s $mut_count $mut_type &
    done
done

# cp ${sims_dir}/different_models/constant__h_0.5__*_sites.txt $traject_dir
# rename 's/constant__/unmodified__/' ${traject_dir}/constant__*
