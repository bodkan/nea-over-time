input_dir=clean_data
sims_dir=simulations

exome_only_exon_coords=${input_dir}/exome_only_exon_coordinates.txt
exome_only_recomb_map=${input_dir}/exome_only_recombination_map.txt

#
# Simulations of AMH and Neanderthal demography until the introgression event.
# Exome only, not including sites from the archaic admixture array.
#
h=0.5
python3 run_mutation_accumulation.py --exon-coordinates $exome_only_exon_coords \
				     --recomb-map $exome_only_recomb_map \
				     --neutral-spacing 10000 \
				     --dominance-coef $h \
				     --output-prefix ${sims_dir}/exome_only__h_$h__ &

h=0.1
python3 run_mutation_accumulation.py --exon-coordinates $exome_only_exon_coords \
				     --recomb-map $exome_only_recomb_map \
				     --neutral-spacing 10000 \
				     --dominance-coef $h \
				     --output-prefix ${sims_dir}/exome_only__h_$h__ &

h=0.0
python3 run_mutation_accumulation.py --exon-coordinates $exome_only_exon_coords \
				     --recomb-map $exome_only_recomb_map \
				     --neutral-spacing 10000 \
				     --dominance-coef $h \
				     --output-prefix ${sims_dir}/exome_only__h_$h__ &
