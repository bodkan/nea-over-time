slim_dir := slim

recomb_map := $(slim_dir)/recombination_map.txt
exon_coords := $(slim_dir)/exon_coordinates.txt
sites_coords := $(slim_dir)/array_sites_coordinates.txt
recomb_map_notebook := notebooks/generate_concat_recomb_map_and_positions_of_exons_and_sites.ipynb

dominance_coefs := 0.5 0.1 0.0


dependencies:
	# generate the recombination map and coordinates of exons and array sites
	jupyter nbconvert $(recomb_map_notebook) --to notebook --execute --ExecutePreprocessor.timeout=-1 --output $(notdir $(recomb_map_notebook));


simulate_populations:
	for h in $(dominance_coefs); do \
		python3 run_mutation_accumulation.py --exon-coordinates $(exon_coords) \
						     --recomb-map $(recomb_map) \
						     --dominance-coef $${h} \
                                                     --output-prefix $(slim_dir)/populations__h_$${h}.txt & \
	done

