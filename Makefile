slim_dir := slim

dominance_coefs := 0.5 0.1 0.0

population_files := $(addsuffix .txt,$(addprefix $(clean_dir)/populations__h_,$(dominance_coefs)))

simulate_populations:
	for h in $(dominance_coefs); do \
		python3 run_mutation_accumulation.py --dominance-coef $${h} --output-prefix $(slim_dir)/populations__h_$${h}.txt --recomb-map $(slim_dir)/recomb_map__exome.txt & \
	done
