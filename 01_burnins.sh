mkdir data/burnins

# chromosomes 1-5 only -------------------------------------------------------

mkdir data/slim_coords/subset

# generate chr1-3 subsets of the coordinate files - used to get faster burnins
for f in data/slim_coords/*.bed; do
    head -n 1 $f > data/slim_coords/subset/`basename ${f}`
    for chr in `seq 1 3`; do
        grep -w "^chr${chr}" $f
    done >> data/slim_coords/subset/`basename ${f}`
done

for region in protein_coding promoter utr3 tf_binding_site; do
for h in `seq 0.0 0.1 1.0`; do
    python3 code/run_mutation_accumulation.py \
        --regions data/slim_coords/subset/${region}_regions.bed \
        --sites data/slim_coords/subset/${region}_all_sites.bed \
        --recomb-map data/slim_coords/subset/${region}_recomb_map.bed \
        --mut-rate 1e-8 \
        --dominance-coef $h \
        --output data/burnins/subset_${region}_h_${h}.txt &
done
done

# uniform recombination rate - to test if the slower recombination rate really
# leads to a more steep decline - as observed in the old simulations
for region in protein_coding promoter utr3 tf_binding_site; do
for h in 0.5; do
    python3 code/run_mutation_accumulation.py \
        --regions data/slim_coords/subset/${region}_regions.bed \
        --sites data/slim_coords/subset/${region}_all_sites.bed \
        --recomb-map data/slim_coords/subset/uniform_${region}_recomb_map.bed \
        --mut-rate 1e-8 \
        --dominance-coef $h \
        --output data/burnins/uniform_subset_${region}_h_${h}.txt &
done
done



# all chromosomes ------------------------------------------------------------

for region in protein_coding promoter utr3 tf_binding_site; do
for h in 0.0 0.5 1.0; do
    python3 code/run_mutation_accumulation.py \
        --regions data/slim_coords/${region}_regions.bed \
        --sites data/slim_coords/${region}_all_sites.bed \
        --recomb-map data/slim_coords/${region}_recomb_map.bed \
        --mut-rate 1e-8 \
        --dominance-coef $h \
        --output data/burnins/${region}_h_${h}.txt &
done
done

for region in protein_coding promoter utr3 tf_binding_site; do
for h in 0.5; do
    python3 code/run_mutation_accumulation.py \
        --regions data/slim_coords/${region}_regions.bed \
        --sites data/slim_coords/${region}_all_sites.bed \
        --recomb-map data/slim_coords/uniform_${region}_recomb_map.bed \
        --mut-rate 1e-8 \
        --dominance-coef $h \
        --output data/burnins/uniform_${region}_h_${h}.txt &
done
done
