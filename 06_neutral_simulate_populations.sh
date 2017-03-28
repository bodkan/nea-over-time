#!/bin/bash

input_dir=input_data
clean_dir=clean_data
sims_dir=simulations
tmp_dir=tmp


#
# Simulations of AMH and Neanderthal demography until the introgression event.
# Exome and sites from the archaic admixture array.
#

exome_and_sites_exon_coords=${clean_dir}/exome_and_sites_exon_coordinates.txt
exome_and_sites_recomb_map=${clean_dir}/exome_and_sites_recombination_map.txt
exonic_array_sites=${clean_dir}/admixture_array_coordinates_exonic.txt
nonexonic_array_sites=${clean_dir}/admixture_array_coordinates_nonexonic.txt


h=1.0
python3 run_mutation_accumulation.py --exon-coordinates $exome_and_sites_exon_coords \
                                     --recomb-map $exome_and_sites_recomb_map \
				     --mut-rate 0 \
                                     --exonic-sites $exonic_array_sites \
                                     --nonexonic-sites $nonexonic_array_sites \
                                     --dominance-coef $h \
                                     --output-prefix ${sims_dir}/exome_and_sites__neutral__ \
                                     --anc-size 1 --hum-nea-split 50 --out-of-africa 25 &






#
# Introgression simulations
#
traject_dir=${sims_dir}/neutral_models

mkdir -p $traject_dir

# initial Neanderthal proportions
init_nea=0.1

# number of replicates of each model
num_replicates=10

# dates of samples
dates=`tail -n+2 ${input_dir}/nea_ancestry_direct.tsv | cut -f2 | tr '\n' ' '`

# wrapper for submitting introgression simulations to SGE
run_introgression() {
    run_id=$8__h_${2}__init_nea_${5}__rep_${1}
    
    ./run_introgression.py \
	    --population-file ${sims_dir}/exome_and_sites__neutral__seed_*.txt \
	    --exon-coordinates $3 \
	    --recomb-map $4 \
	    --mut-rate 0 \
	    --admixture-rate $5 \
	    --dominance-coef $2 \
	    --sampling-times $dates \
	    --exonic-sites $6 \
	    --nonexonic-sites $7 \
	    --output-prefix ${traject_dir}/${run_id} \
	    --model $8 \
	    --save-trajectories \
	    $9
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
for model in 'constant'; do #'constant' 'linear' 'gravel'; do
    run_introgression 1 $dominance_coef $exome_only_exon_coords $exome_only_recomb_map $init_nea $exonic_sites $nonexonic_sites $model --save-mutations &
    # for rep_i in `seq 2 $num_replicates`; do
    #     run_introgression $rep_i $dominance_coef $exome_only_exon_coords $exome_only_recomb_map $init_nea $exonic_sites $nonexonic_sites $model &
    # done
done
