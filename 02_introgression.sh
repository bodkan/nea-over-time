burnin_dir=data/burnins
sim_dir=data/simulations

# TF was failing, probably a corrupted burnin file - I run the simulations
# from the whole_burnin TF files
for region in promoter tf_binding_site protein_coding; do
for h in 0.0 0.5 1.0; do
for rep in `seq 1 5`; do
    if [ $rep -eq 1 ]; then vcf_opt="--vcf-times 2200"; fi
    python3 code/run_introgression.py \
        --regions data/slim_coords/${region}_regions.bed \
        --sites data/slim_coords/${region}_all_sites.bed \
        --recomb-map data/slim_coords/${region}_recomb_map.bed \
        --mut-rate 1e-8 \
        --dominance-coef $h \
        --model constant \
        --output-prefix ${sim_dir}/${region}_h_${h}_rep_${rep} \
        --population-file ${burnin_dir}/${region}_h_${h}.txt \
        $vcf_opt &
done
done
done

cd ../slim-neanderthal
for h in 0.5; do
python3 run_introgression.py \
    --exon-coordinates /mnt/scratch/mp/slim-neanderthal/clean_data/exome_and_sites_exon_coordinates.txt \
    --exonic-sites /mnt/scratch/mp/slim-neanderthal/clean_data/admixture_array_coordinates_exonic.txt \
    --nonexonic-sites /mnt/scratch/mp/slim-neanderthal/clean_data/admixture_array_coordinates_nonexonic.txt \
    --recomb-map /mnt/scratch/mp/slim-neanderthal/clean_data/exome_and_sites_recombination_map.txt \
    --dominance-coef $h \
    --model constant \
    --terminate-after 210 \
    --save-mutations \
    --dump-slim asd.txt \
    --output-prefix ../nea-over-time/data/simulations/nea_vs_mh_${h} \
    --population-file /mnt/scratch/mp/slim-neanderthal/simulations/exome_and_sites__h_${h}__seed_* &
done
cd ../nea-over-time

# $ 2017-12-01 09:29:29 :: INFO :: Running simulation from SLiM input file "/tmp/tmpzvgjkb93"
# 2017-12-01 09:29:29 :: INFO :: Running simulation from SLiM input file "/tmp/tmpm77c43ut"
# 2017-12-01 09:29:29 :: INFO :: Running simulation from SLiM input file "/tmp/tmpzcza_60d"







# for region in protein_coding promoter utr3 tf_binding_site; do
# for h in 0.5; do
# for rep in `seq 4 10`; do
# python3 code/run_introgression.py \
#     --regions data/slim_coords/subset/${region}_regions.bed \
#     --sites data/slim_coords/subset/${region}_all_sites.bed \
#     --recomb-map data/slim_coords/subset/${region}_recomb_map.bed \
#     --mut-rate 1e-8 \
#     --dominance-coef $h \
#     --model constant \
#     --terminate-after 500 \
#     --output-prefix ${sim_dir}/subset_${region}_h_${h}_rep_${rep} \
#     --population-file ${burnin_dir}/subset_${region}_h_${h}.txt &
# done
# done
# done

# for region in protein_coding promoter utr3 tf_binding_site; do
# for h in 0.5; do
# for rep in `seq 4 10`; do
# python3 code/run_introgression.py \
#     --regions data/slim_coords/subset/${region}_regions.bed \
#     --sites data/slim_coords/subset/${region}_all_sites.bed \
#     --recomb-map data/slim_coords/subset/uniform_${region}_recomb_map.bed \
#     --mut-rate 1e-8 \
#     --dominance-coef $h \
#     --model constant \
#     --terminate-after 500 \
#     --output-prefix ${sim_dir}/uniform_subset_${region}_h_${h}_rep_${rep} \
#     --population-file ${burnin_dir}/uniform_subset_${region}_h_${h}.txt &
# done
# done
# done




