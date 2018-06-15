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

region="exon";
for h in 0.5; do
for rep in `seq 1 20`; do
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
done > /dev/null
done



# ----------------------------------------------------------------------
# simulate Neanderthal and Denisovan deserts

mkdir -p data/deserts
for region in exon; do
for chrom in chr1 chr7; do
for source in p2 p4; do
for rep in `seq 1 100`; do
    N="${chrom}_${source}_${region}_${rep}"
    qsub -V -cwd -j y -l virtual_free=30G,h_vmem=30G -N $N -o tmp/${N}.txt \
    ./code/run_introgression.py \
        --chrom $chrom \
        --regions data/slim_coords/${region}_unif_regions.bed \
        --sites data/slim_coords/${region}_unif_all_sites.bed \
        --recomb-map data/slim_coords/${region}_unif_recomb_map.bed \
        --mut-rate 7e-9 \
        --dominance-coef 0.5 \
        --model constant \
        --admixture-rate 0.025 \
        --admixture-source ${source} \
        --output-prefix data/deserts/selection_${chrom}_source_${source}_${region}_rep_${rep} \
        --population-file data/burnins/${chrom}_${region}_${rep}.txt \
        --vcf-times 2200 \
        --vcf-sample 500
done
done
done
done

for region in exon; do
for chrom in chr1 chr7; do
for source in p2 p4; do
for rep in `seq 1 100`; do
    N="neutral_${chrom}_${source}_${region}_${rep}"
    qsub -V -cwd -j y -l virtual_free=30G,h_vmem=30G -N $N -o tmp/${N}.txt \
    ./code/run_introgression.py \
        --chrom $chrom \
        --regions data/slim_coords/${region}_unif_regions.bed \
        --sites data/slim_coords/${region}_unif_all_sites.bed \
        --recomb-map data/slim_coords/${region}_unif_recomb_map.bed \
        --mut-rate 0 \
	--force-neutral \
        --dominance-coef 0.5 \
        --model constant \
        --admixture-rate 0.025 \
        --admixture-source ${source} \
        --output-prefix data/deserts/neutral_${chrom}_source_${source}_${region}_rep_${rep} \
        --population-file data/burnins/${chrom}_${region}_${rep}.txt \
        --vcf-times 2200 \
        --vcf-sample 500
done
done
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
