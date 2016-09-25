slim_dir := slim
sims_dir := simulations

recomb_map := $(slim_dir)/recombination_map.txt
exon_coords := $(slim_dir)/exon_coordinates.txt
sites_coords := $(slim_dir)/array_sites_coordinates.txt

dominance_coefs := 0.5 0.1 0.0


simulate_populations:
	for h in $(dominance_coefs); do \
		python3 run_mutation_accumulation.py --exon-coordinates $(exon_coords) \
						     --recomb-map $(recomb_map) \
						     --dominance-coef $${h} \
                                                     --output-prefix $(sims_dir)/populations__h_$${h}.txt & \
	done

