#!/usr/bin/env bash

mkdir data/simulations



# ----------------------------------------------------------------------
# simulations of Neandertal trajectories over time
#

#
# human demographic models
#
region="exon"; h=0.5
for model in gravel linear constant; do
    for rep in `seq 1 30`; do
        N="${model}_${rep}"
        qsub -V -cwd -j y -l virtual_free=50G,h_vmem=50G -N $N -o tmp/${N}.txt \
        ./code/run_introgression.py \
            --regions data/slim_coords/${region}_regions.bed \
            --sites data/slim_coords/${region}_all_sites.bed \
            --recomb-map data/slim_coords/${region}_recomb_map.bed \
            --mut-rate 1e-8 \
            --dominance-coef 0.5 \
            --model $model \
            --gap-trajectories \
            --output-prefix data/simulations/traj_${model}_${region}_rep_${rep} \
            --population-file data/burnins/${region}_h_${h}.txt
    done
done

#
# deleterious regions of different sizes
#
h=0.5
for region in exon promoter tf_binding_site protein_coding utr3; do
    for rep in `seq 1 30`; do
        N="${region}_${rep}"
        qsub -V -cwd -j y -l virtual_free=50G,h_vmem=50G -N $N -o tmp/${N}.txt \
        ./code/run_introgression.py \
            --regions data/slim_coords/${region}_regions.bed \
            --sites data/slim_coords/${region}_all_sites.bed \
            --recomb-map data/slim_coords/${region}_recomb_map.bed \
            --mut-rate 1e-8 \
            --dominance-coef $h \
            --model constant \
            --gap-trajectories \
            --output-prefix data/simulations/traj_${region}_rep_${rep} \
            --population-file data/burnins/${region}_h_${h}.txt
    done
done

#
# different sizes of Neandertal population prior to the introgression
#
region="exon"; h=0.5
for Ne in 100 500 1000 10000; do
    for rep in `seq 1 30`; do # some of these jobs got killed
        N="Ne_${Ne}_${rep}"
        qsub -V -cwd -j y -l virtual_free=80G,h_vmem=80G -N $N -o tmp/${N}.txt \
        ./code/run_introgression.py \
            --regions data/slim_coords/${region}_regions.bed \
            --sites data/slim_coords/${region}_all_sites.bed \
            --recomb-map data/slim_coords/${region}_recomb_map.bed \
            --mut-rate 1e-8 \
            --dominance-coef 0.5 \
            --model constant \
            --gap-trajectories \
            --output-prefix data/simulations/traj_Ne_${Ne}_${region}_rep_${rep} \
            --population-file data/burnins/nea_Ne_${Ne}_${region}_h_${h}.txt
    done
done


#
# artificially making Neandertal mutations more deleterious
#
region="exon"; h=0.5
for modifier in 1.0 1.1 1.25 1.5 1.75 2.0 5.0 10.0; do
    for rep in `seq 1 10`; do
        N="modif_${modifier}_${rep}"
        qsub -V -cwd -j y -l virtual_free=50G,h_vmem=50G -N $N -o tmp/${N}.txt \
        ./code/run_introgression.py \
            --regions data/slim_coords/${region}_regions.bed \
            --sites data/slim_coords/${region}_all_sites.bed \
            --recomb-map data/slim_coords/${region}_recomb_map.bed \
            --mut-rate 1e-8 \
            --modify-fraction 1.0 \
            --multiply-s $modifier \
            --dominance-coef $h \
            --model constant \
            --gap-trajectories \
            --output-prefix data/simulations/traj_mult_${modifier}_rep_${rep} \
            --population-file data/burnins/${region}_h_${h}.txt
    done
done

# dominance mix trajectories
for prop_add in `seq 0 0.1 1`; do
    prop_rec=`echo "1 $prop_add" | awk '{print $1 - $2}'`
    for rep in `seq 1 20`; do
        N="mix_add_${prop_add}_rec_${prop_rec}_rep_${rep}"
        qsub -V -cwd -j y -l virtual_free=50G,h_vmem=50G -N $N -o tmp/${N}.txt \
            ./code/run_slim.sh "-d prop_add=$prop_add -d prop_rec=$prop_rec -d output_file='data/simulations/traj_${N}.txt' code/slim/dominance_mix__introgression.slim"
    done
done


# ----------------------------------------------------------------------
# simulations for analysis of frequency derivatives over time
#

# (ran 50 reps for exons, in fact)
for region in exon promoter tf_binding_site protein_coding utr3; do
for h in 0.5; do
for rep in `seq 1 20`; do
    N="${region}_${rep}"
    qsub -V -cwd -j y -l virtual_free=60G,h_vmem=60G -N $N -o tmp/${N}.txt \
    ./code/run_introgression.py \
        --regions data/slim_coords/${region}_regions.bed \
        --sites data/slim_coords/${region}_all_sites.bed \
        --recomb-map data/slim_coords/${region}_recomb_map.bed \
        --mut-rate 1e-8 \
        --dominance-coef $h \
        --model constant \
        --output-prefix data/simulations/delta_constant_${region}_h_${h}_rep_${rep} \
        --population-file data/burnins/${region}_h_${h}.txt \
        --vcf-times 1 2 3 4 5 6 7 8 9 10 20 50 100 `seq 200 200 2200` \
        --vcf-sample 500
done > /dev/null
done
done

# ----------------------------------------------------------------------
# simulations for comparison of real Neandertal ancestry in functional
# regions with distance-to-genes analysis

# multiplying selection coefficients of regions under selection
# to compare Nea ancestry levels to those in gap regions
region="exon"
h=0.5
for mult in 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0; do
for rep in `seq 1 20`; do
    N="${region}_${mult}_${rep}"
    qsub -V -cwd -j y -l virtual_free=60G,h_vmem=60G -N $N -o tmp/${N}.txt \
    ./code/run_introgression.py \
        --regions data/slim_coords/${region}_regions.bed \
        --sites data/slim_coords/${region}_all_sites.bed \
        --recomb-map data/slim_coords/${region}_recomb_map.bed \
        --mut-rate 1e-8 \
        --multiply-s $mult \
        --modify-fraction 1.0 \
        --dominance-coef $h \
        --model constant \
        --output-prefix data/simulations/${region}_mult_${mult}_h_${h}_rep_${rep} \
        --population-file data/burnins/${region}_h_${h}.txt \
        --vcf-times 2200 \
        --vcf-sample 500
done > /dev/null
done

# different levels of purifying selection in different classes of exons
# (each exon assigned into one of random bins)
region="exon"
for bin in bin_s bin_h; do
for rep in `seq 1 50`; do
    N="${region}_${bin}_${rep}"
    qsub -V -cwd -j y -l virtual_free=60G,h_vmem=60G -N $N -o tmp/${N}.txt \
    ./code/run_introgression.py \
        --regions data/slim_coords/${bin}_${region}_regions.bed \
        --sites data/slim_coords/${region}_all_sites.bed \
        --recomb-map data/slim_coords/${region}_recomb_map.bed \
        --mut-rate 1e-8 \
        --model constant \
        --output-prefix data/simulations/${region}_${bin}_rep_${rep} \
        --population-file data/burnins/${bin}_${region}.txt \
        --vcf-times 2200 \
        --vcf-sample 500
done > /dev/null
done

# different dominance of mutations in promoters and protein-coding regions
region="merged"
for h in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0; do
for rep in `seq 1 20`; do
    N="different_h_${h}_${rep}"
    qsub -V -cwd -j y -l virtual_free=40G,h_vmem=40G -N $N -o tmp/${N}.txt \
    ./code/run_introgression.py \
        --regions data/slim_coords/different_h_${h}_${region}_regions.bed \
        --sites data/slim_coords/${region}_all_sites.bed \
        --recomb-map data/slim_coords/${region}_recomb_map.bed \
        --mut-rate 1e-8 \
        --model constant \
        --output-prefix data/simulations/${region}_dom_${h}_rep_${rep} \
        --population-file data/burnins/different_h_${h}_${region}.txt \
        --vcf-times 2200 \
        --vcf-sample 500
done > /dev/null
done

# different selection coefficients of mutations in promoters and protein-coding regions
region="merged"
for s in 0.0 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0; do
for rep in `seq 1 20`; do
    N="different_s_${s}_${rep}"
    qsub -V -cwd -j y -l virtual_free=50G,h_vmem=50G -N $N -o tmp/${N}.txt \
    ./code/run_introgression.py \
        --regions data/slim_coords/different_s_${s}_${region}_regions.bed \
        --sites data/slim_coords/${region}_all_sites.bed \
        --recomb-map data/slim_coords/${region}_recomb_map.bed \
        --mut-rate 1e-8 \
        --model constant \
        --output-prefix data/simulations/${region}_sel_${s}_rep_${rep} \
        --population-file data/burnins/different_s_${s}_${region}.txt \
        --vcf-times 2200 \
        --vcf-sample 500
done > /dev/null
done
