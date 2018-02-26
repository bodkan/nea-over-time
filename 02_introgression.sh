mkdir data/simulations



# ----------------------------------------------------------------------
# simulations for trajectories for different regions over time

for model in gravel linear constant; do
for region in merged exon promoter tf_binding_site protein_coding utr3; do # not all h for merged simulated
for h in 0.0 0.5 1.0; do
for rep in `seq 1 5`; do
    python3 code/run_introgression.py \
        --regions data/slim_coords/${region}_regions.bed \
        --sites data/slim_coords/${region}_all_sites.bed \
        --recomb-map data/slim_coords/${region}_recomb_map.bed \
        --mut-rate 1e-8 \
        --dominance-coef $h \
        --model $model \
        --output-prefix data/simulations/traj_${model}_${region}_h_${h}_rep_${rep} \
        --population-file data/burnins/${region}_h_${h}.txt &
done
done
done
done



# ----------------------------------------------------------------------
# simulations for analysis of frequency derivatives over time

# constant model of exonic selection

region="exon"; h=0.5
for rep in `seq 1 3`; do
    python3 code/run_introgression.py \
        --regions data/slim_coords/${region}_regions.bed \
        --sites data/slim_coords/${region}_all_sites.bed \
        --recomb-map data/slim_coords/${region}_recomb_map.bed \
        --mut-rate 1e-8 \
        --dominance-coef $h \
        --model constant \
        --output-prefix data/simulations/delta_constant_${region}_h_${h}_rep_${rep} \
        --population-file data/burnins/${region}_h_${h}.txt \
        --vcf-times 1 2 3 4 5 6 7 8 9 10 20 50 100 `seq 200 200 2200` \
        --vcf-sample 500 &
done

# constant model of protein coding selection

region="protein_coding"; h=0.5
for rep in `seq 1 5`; do
    python3 code/run_introgression.py \
        --regions data/slim_coords/${region}_regions.bed \
        --sites data/slim_coords/${region}_all_sites.bed \
        --recomb-map data/slim_coords/${region}_recomb_map.bed \
        --mut-rate 7e-9 \
        --dominance-coef $h \
        --admixture-rate 0.05 \
        --model constant \
        --output-prefix data/simulations/delta_constant_${region}_h_${h}_rep_${rep} \
        --population-file data/burnins/nonsyn_${region}_h_${h}.txt \
        --vcf-times 1 2 3 4 5 6 7 8 9 10 20 50 `seq 100 100 1000` `seq 1200 200 2200` \
        --vcf-sample 500 &
done > /dev/null



# ----------------------------------------------------------------------
# analysis of desert sizes in the present based

# different amount of deleterious sequence

h=0.5
for region in merged exon promoter tf_binding_site protein_coding utr3; do
for rep in `seq 1 5`; do
    python3 code/run_introgression.py \
        --regions data/slim_coords/${region}_regions.bed \
        --sites data/slim_coords/${region}_all_sites.bed \
        --recomb-map data/slim_coords/${region}_recomb_map.bed \
        --mut-rate 1e-8 \
        --dominance-coef $h \
        --model constant \
        --output-prefix data/simulations/deserts_${region}_h_${h}_rep_${rep} \
        --population-file data/burnins/${region}_h_${h}.txt \
        --vcf-times 2200 \
        --vcf-sample 500 &
done
done



# under neutrality

region="tf_binding_site" # size of deleterious region doesn't matter here
h=0.5
for Ne in 10000; do
for rep in `seq 1 5`; do
    python3 code/run_introgression.py \
        --regions data/slim_coords/${region}_regions.bed \
        --sites data/slim_coords/${region}_all_sites.bed \
        --recomb-map data/slim_coords/${region}_recomb_map.bed \
        --force-neutral \
        --dominance-coef $h \
        --model constant \
        --founder-size $Ne \
        --output-prefix data/simulations/deserts_neutral_Ne_${Ne}_h_${h}_rep_${rep} \
        --population-file data/burnins/${region}_h_${h}.txt \
        --vcf-sample 500 \
        --vcf-times 2200 &
done
done



# ----------------------------------------------------------------------
# influence of continuous admixture (as opposed to a "single pulse")
# on the distribution of desert sizes

#
# analysis of desert sizes under weak negative selection model
#
h=0.5
for region in exon promoter tf_binding_site protein_coding utr3; do
for rep in `seq 1 3`; do
    python3 code/run_introgression.py \
        --regions data/slim_coords/${region}_regions.bed \
        --sites data/slim_coords/${region}_all_sites.bed \
        --recomb-map data/slim_coords/${region}_recomb_map.bed \
        --mut-rate 1e-8 \
        --dominance-coef $h \
        --model constant \
        --admixture-rate 0.0025 \
        --admixture-end 54000 \
        --output-prefix data/simulations/deserts_continuous_${region}_h_${h}_rep_${rep} \
        --population-file data/burnins/${region}_h_${h}.txt \
        --vcf-times 2200 \
        --vcf-sample 500 &
done
done

#
# analysis of desert sizes in the present under neutrality
#
region="tf_binding_site" # size of deleterious region doesn't matter here
h=0.5
for Ne in 10000; do
for rep in `seq 1 3`; do
    python3 code/run_introgression.py \
        --regions data/slim_coords/${region}_regions.bed \
        --sites data/slim_coords/${region}_all_sites.bed \
        --recomb-map data/slim_coords/${region}_recomb_map.bed \
        --force-neutral \
        --dominance-coef $h \
        --model constant \
        --founder-size $Ne \
        --admixture-rate 0.0025 \
        --admixture-end 54000 \
        --output-prefix data/simulations/deserts_continuous_neutral_Ne_${Ne}_h_${h}_rep_${rep} \
        --population-file data/burnins/${region}_h_${h}.txt \
        --vcf-sample 500 \
        --vcf-times 2200 &
done
done


# cd ../slim-neanderthal
# for h in 0.5; do
# python3 run_introgression.py \
#     --exon-coordinates /mnt/scratch/mp/slim-neanderthal/clean_data/exome_and_sites_exon_coordinates.txt \
#     --exonic-sites /mnt/scratch/mp/slim-neanderthal/clean_data/admixture_array_coordinates_exonic.txt \
#     --nonexonic-sites /mnt/scratch/mp/slim-neanderthal/clean_data/admixture_array_coordinates_nonexonic.txt \
#     --recomb-map /mnt/scratch/mp/slim-neanderthal/clean_data/exome_and_sites_recombination_map.txt \
#     --dominance-coef $h \
#     --model constant \
#     --terminate-after 210 \
#     --save-mutations \
#     --dump-slim asd.txt \
#     --output-prefix ../nea-over-time/data/simulations/nea_vs_mh_${h} \
#     --population-file /mnt/scratch/mp/slim-neanderthal/simulations/exome_and_sites__h_${h}__seed_* &
# done
# cd ../nea-over-time

# $ 2017-12-01 09:29:29 :: INFO :: Running simulation from SLiM input file "/tmp/tmpzvgjkb93"
# 2017-12-01 09:29:29 :: INFO :: Running simulation from SLiM input file "/tmp/tmpm77c43ut"
# 2017-12-01 09:29:29 :: INFO :: Running simulation from SLiM input file "/tmp/tmpzcza_60d"
