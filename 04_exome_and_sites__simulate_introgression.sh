#!/bin/bash

input_dir=input_data
clean_dir=clean_data
sims_dir=simulations
traject_dir=${sims_dir}/different_models
tmp_dir=tmp

mkdir -p $traject_dir


# initial Neanderthal proportions
init_nea=0.1

# number of replicates of each model
num_replicates=10

# dates of samples
dates=`tail -n+2 ${input_dir}/nea_ancestry_direct.tsv | cut -f2 | tr '\n' ' '`

# wrapper for submitting introgression simulations to SGE
run_introgression() {
    case $2 in
	0.5)
	    mem=8G
	    ;;
	0.1)
	    mem=10G
	    ;;
	0.0)
	    mem=12G
	    ;;
	*)
	    echo "Suspicious dominance coefficient: $2"
	    exit
    esac

    run_id=$8__h_${2}__init_nea_${5}__rep_${1}
    
    qsub -V -cwd -j y -S /bin/bash -l virtual_free=$mem,h_vmem=$mem \
	 -o ${tmp_dir}/sge__${run_id}.out -N ${run_id} \
	 run_introgression.py \
	    --population-file ${sims_dir}/exome_and_sites__h_${2}__seed_*.txt \
	    --exon-coordinates $3 \
	    --recomb-map $4 \
	    --admixture-rate $5 \
	    --dominance-coef $2 \
	    --sampling-times $dates \
	    --exonic-sites $6 \
	    --nonexonic-sites $7 \
	    --output-prefix ${traject_dir}/${run_id} \
	    --model $8 \
	    $9
}


# 
# Simulations of Neanderthal ancestry trajectories in the European
# population using previously simulated populations with accumulated
# deleterious variants.
#
# Neutral Neanderthal markers are exonic and non-exonic sites from
# the archaic admixture array. Dominance coefficient 0.5 only.
#

exome_only_exon_coords=${clean_dir}/exome_and_sites_exon_coordinates.txt
exome_only_recomb_map=${clean_dir}/exome_and_sites_recombination_map.txt
exonic_sites=${clean_dir}/admixture_array_coordinates_exonic.txt
nonexonic_sites=${clean_dir}/admixture_array_coordinates_nonexonic.txt

for model in 'constant' 'linear' 'gravel'; do
    run_introgression 1 0.5 $exome_only_exon_coords $exome_only_recomb_map $init_nea $exonic_sites $nonexonic_sites $model --save-mutations
    for rep_i in `seq 2 $num_replicates`; do
        run_introgression $rep_i 0.5 $exome_only_exon_coords $exome_only_recomb_map $init_nea $exonic_sites $nonexonic_sites $model
    done
done
