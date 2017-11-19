burnin_dir=data/burnins
sim_dir=data/simulations

region="protein_coding"
h=0.5

for rep in `seq 1 10`; do
python3 code/run_introgression.py \
    --regions data/slim_coords/${region}_regions.bed \
    --sites data/slim_coords/${region}_sites.bed \
    --recomb-map data/slim_coords/${region}_recomb_map.bed \
    --mut-rate 1e-8 \
    --dominance-coef $h \
    --model constant \
    --terminate-after 500 \
    --output-prefix ${sim_dir}/${region}_${h}_rep_${rep} \
    --population-file ${burnin_dir}/${region}_${h}.out.txt &
done

cd ../slim-neanderthal
for rep in `seq 1 10`; do
python3 run_introgression.py \
    --exon-coordinates /mnt/scratch/mp/slim-neanderthal/clean_data/exome_and_sites_exon_coordinates.txt \
    --exonic-sites /mnt/scratch/mp/slim-neanderthal/clean_data/admixture_array_coordinates_exonic.txt \
    --nonexonic-sites /mnt/scratch/mp/slim-neanderthal/clean_data/admixture_array_coordinates_nonexonic.txt \
    --recomb-map /mnt/scratch/mp/slim-neanderthal/clean_data/exome_and_sites_recombination_map.txt \
    --dominance-coef $h \
    --model constant \
    --terminate-after 500 \
    --save-trajectories \
    --output-prefix ../nea-over-time/data/simulations/old_${region}_${h}_rep_${rep} \
    --population-file /mnt/scratch/mp/slim-neanderthal/simulations/exome_and_sites__h_0.5__seed_6977220333793.txt &
done
cd ../nea-over-time
