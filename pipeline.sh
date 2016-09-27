#!/bin/bash

input_dir=clean_data
sims_dir=simulations



# 
# Simulations of AMH and Neanderthal demography until the introgression event.
# Exome only, not including sites from the archaic admixture array.
#

exome_only_exon_coords=${input_dir}/exome_only_exon_coordinates.txt
exome_only_recomb_map=${input_dir}/exome_only_recombination_map.txt


h=0.5
python3 run_mutation_accumulation.py --exon-coordinates $exome_only_exon_coords \
				     --recomb-map $exome_only_recomb_map \
				     --neutral-spacing 10000 \
				     --dominance-coef $h \
				     --output-prefix ${sims_dir}/exome_only__h_${h}__ &

h=0.1
python3 run_mutation_accumulation.py --exon-coordinates $exome_only_exon_coords \
				     --recomb-map $exome_only_recomb_map \
				     --neutral-spacing 10000 \
				     --dominance-coef $h \
				     --output-prefix ${sims_dir}/exome_only__h_${h}__ &


h=0.0
python3 run_mutation_accumulation.py --exon-coordinates $exome_only_exon_coords \
				     --recomb-map $exome_only_recomb_map \
				     --neutral-spacing 10000 \
				     --dominance-coef $h \
				     --output-prefix ${sims_dir}/exome_only__h_${h}__ &



#
# Simulations of AMH and Neanderthal demography until the introgression event.
# Exome only, not including sites from the archaic admixture array.
# Using SLiM recombination map by Harris and Nielsen, 2016.
#

harris_config=slim/execute_mut_accum_partially_recessive.slim
harris_recomb_map=${input_dir}/harris_recombination_map.txt
harris_exon_coords=${input_dir}/harris_exon_coordinates.txt

# first extract the recombination map used by the H&K paper
# (converting 1-based coordinates of SliM 1.8 to 0-based coordinates
# of SLiM 2)
echo -e "slim_end\trecomb_rate" > $harris_recomb_map
awk -vOFS="\t" 'NR > 12 && NR < 390986 {print $1 - 1, $2}' $harris_config \
    >> $harris_recomb_map

# get the coordinates of the simulated genomic element
segment_length=`tail -1 $harris_recomb_map | cut -f1`
echo -e "slim_start\tslim_end" > $harris_exon_coords
echo -e "0\t${segment_length}" >> $harris_exon_coords

# then run the same mutation accumulation simulations as above
h=0.5
python3 run_mutation_accumulation.py --exon-coordinates $harris_exon_coords \
				     --recomb-map $harris_recomb_map \
				     --neutral-spacing 10000 \
				     --dominance-coef $h \
 				     --output-prefix ${sims_dir}/harris__exome_only__h_${h}__ &

h=0.1
python3 run_mutation_accumulation.py --exon-coordinates $harris_exon_coords \
				     --recomb-map $harris_recomb_map \
				     --neutral-spacing 10000 \
				     --dominance-coef $h \
				     --output-prefix ${sims_dir}/harris__exome_only__h_${h}__ &

h=0.0
python3 run_mutation_accumulation.py --exon-coordinates $harris_exon_coords \
				     --recomb-map $harris_recomb_map \
				     --neutral-spacing 10000 \
				     --dominance-coef $h \
				     --output-prefix ${sims_dir}/harris__exome_only__h_${h}__ &
