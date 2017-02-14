# SLiM script were madually edited and simplified from the templates
# filled in by the wrapper script run_mutation_accumulation.py

slim -d prop_add=1.0 -d prop_partrec=0.0 -d prop_rec=0.0 -d segment_length=10000000 -d neutral_spacing=1000 slim/dominance_mix__mut_accum.slim &
slim -d prop_add=0.0 -d prop_partrec=1.0 -d prop_rec=0.0 -d segment_length=10000000 -d neutral_spacing=1000 slim/dominance_mix__mut_accum.slim &
slim -d prop_add=0.0 -d prop_partrec=0.0 -d prop_rec=1.0 -d segment_length=10000000 -d neutral_spacing=1000 slim/dominance_mix__mut_accum.slim &


slim -d prop_add=0.75 -d prop_partrec=0.0 -d prop_rec=0.25 -d segment_length=10000000 -d neutral_spacing=1000 slim/dominance_mix__mut_accum.slim &
slim -d prop_add=0.5 -d prop_partrec=0.0 -d prop_rec=0.5 -d segment_length=10000000 -d neutral_spacing=1000 slim/dominance_mix__mut_accum.slim &
slim -d prop_add=0.25 -d prop_partrec=0.0 -d prop_rec=0.75 -d segment_length=10000000 -d neutral_spacing=1000 slim/dominance_mix__mut_accum.slim &
