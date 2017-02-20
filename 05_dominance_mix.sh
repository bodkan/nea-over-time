# SLiM scripts were madually edited and simplified from the templates
# filled in by the wrapper scripts run_mutation_accumulation.py and
# run_introgression.py

#
# mutation accumulation simulations
#
slim -d prop_add=1.0 -d prop_partrec=0.0 -d prop_rec=0.0 -d segment_length=10000000 -d neutral_spacing=1000 slim/dominance_mix__mut_accum.slim &
slim -d prop_add=0.0 -d prop_partrec=1.0 -d prop_rec=0.0 -d segment_length=10000000 -d neutral_spacing=1000 slim/dominance_mix__mut_accum.slim &
slim -d prop_add=0.0 -d prop_partrec=0.0 -d prop_rec=1.0 -d segment_length=10000000 -d neutral_spacing=1000 slim/dominance_mix__mut_accum.slim &


slim -d prop_add=0.75 -d prop_partrec=0.0 -d prop_rec=0.25 -d segment_length=10000000 -d neutral_spacing=1000 slim/dominance_mix__mut_accum.slim &
slim -d prop_add=0.5  -d prop_partrec=0.0 -d prop_rec=0.5  -d segment_length=10000000 -d neutral_spacing=1000 slim/dominance_mix__mut_accum.slim &
slim -d prop_add=0.25 -d prop_partrec=0.0 -d prop_rec=0.75 -d segment_length=10000000 -d neutral_spacing=1000 slim/dominance_mix__mut_accum.slim &


#
# introgression simulations
#

NUM_REPLICATES=30

sim_replicates() {
    for rep in `seq 1 $NUM_REPLICATES`; do
	slim -d prop_add=$1 -d prop_partrec=0.0 -d prop_rec=$2 -d segment_length=10000000 -d neutral_spacing=1000 -d pop_size=$3 -d rep=$rep slim/dominance_mix__introgression.slim
    done
}

sim_replicates 1.0  0.0  2000 &
sim_replicates 0.25 0.75 2000 &
sim_replicates 0.5  0.5  2000 &
sim_replicates 0.75 0.25 2000 &
sim_replicates 0.0  1.0  2000 &


sim_replicates 1.0  0.0  10000 &
sim_replicates 0.25 0.75 10000 &
sim_replicates 0.5  0.5  10000 &
sim_replicates 0.75 0.25 10000 &
sim_replicates 0.0  1.0  10000 &
