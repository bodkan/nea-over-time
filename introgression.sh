#!/bin/bash

input_dir=input_data
clean_dir=clean_data
sims_dir=simulations
traject_dir=${sims_dir}/trajectories
tmp_dir=tmp

mkdir -p $traject_dir


# initial Neanderthal proportions
init_props="0.01 0.025 0.05 0.075 0.1"

# number of replicates of each model
num_replicates=50

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

    run_id=exome_only__h_${2}__init_nea_${5}__rep_${1}
    
    qsub -V -cwd -j y -S /bin/bash -l virtual_free=$mem,h_vmem=$mem \
	 -o ${tmp_dir}/sge__${run_id}.out -N ${run_id} \
	 run_introgression.py \
	    --population-file ${sims_dir}/exome_only__h_${2}__seed_*.txt \
	    --exon-coordinates $3 \
	    --recomb-map $4 \
	    --neutral-spacing 10000 \
	    --admixture-rate $5 \
	    --dominance-coef $2 \
	    --sampling-times $dates \
	    --output-file ${traject_dir}/${run_id}.txt
}


# 
# Simulations of Neanderthal ancestry trajectories in the European
# population using previously simulated populations with accumulated
# deleterious variants.
#
# Exome only, not including sites from the archaic admixture array.
#

exome_only_exon_coords=${clean_dir}/exome_only_exon_coordinates.txt
exome_only_recomb_map=${clean_dir}/exome_only_recombination_map.txt

for rep_i in `seq 1 $num_replicates`; do
    for init_nea in $init_props; do
	run_introgression $rep_i 0.5 $exome_only_exon_coords $exome_only_recomb_map $init_nea
	run_introgression $rep_i 0.1 $exome_only_exon_coords $exome_only_recomb_map $init_nea
	run_introgression $rep_i 0.0 $exome_only_exon_coords $exome_only_recomb_map $init_nea
    done
done


# 
# Simulations of Neanderthal ancestry trajectories in the European
# population using previously simulated populations with accumulated
# deleterious variants.
#
# Exome only, not including sites from the archaic admixture array.
#
# Using SLiM recombination map by Harris and Nielsen, 2016.
#

# harris_recomb_map=${clean_dir}/harris_recombination_map.txt
# harris_exon_coords=${clean_dir}/harris_exon_coordinates.txt


# for rep_i in `seq 1 $num_replicates`; do
#     for init_nea in $init_props; do
# 	run_introgression $rep_i 0.5 $harris_exon_coords $harris_recomb_map $init_nea
# 	run_introgression $rep_i 0.1 $harris_exon_coords $harris_recomb_map $init_nea
# 	run_introgression $rep_i 0.0 $harris_exon_coords $harris_recomb_map $init_nea
#     done
# done
