# SLiM scripts were madually edited and simplified from the templates
# filled in by the wrapper scripts run_mutation_accumulation.py and
# run_introgression.py

SEGMENT_LENGTH=10000000
NEUTRAL_SPACING=1000

#
# mutation accumulation simulations
#

mut_accum() {
    slim -d prop_add=$1                       \
	 -d prop_partrec=$2                   \
	 -d prop_rec=$3                       \
	 -d segment_length=$SEGMENT_LENGTH    \
	 -d neutral_spacing=$NEUTRAL_SPACING  \
	 -d recomb_rate=$4                    \
	 slim/dominance_mix__mut_accum.slim
}

# recombination rate 1e-8
mut_accum 1.0  0.0 0.0  1e-8 &
mut_accum 0.0  1.0 0.0  1e-8 &
mut_accum 0.0  0.0 1.0  1e-8 &
mut_accum 0.75 0.0 0.25 1e-8 &
mut_accum 0.5  0.0 0.5  1e-8 &
mut_accum 0.25 0.0 0.75 1e-8 &

# recombination rate 1e-7
mut_accum 1.0  0.0 0.0  1e-7 &
mut_accum 0.0  1.0 0.0  1e-7 &
mut_accum 0.0  0.0 1.0  1e-7 &
mut_accum 0.75 0.0 0.25 1e-7 &
mut_accum 0.5  0.0 0.5  1e-7 &
mut_accum 0.25 0.0 0.75 1e-7 &


# recombination rate 1e-6
mut_accum 1.0  0.0 0.0  1e-6 &
mut_accum 0.0  1.0 0.0  1e-6 &
mut_accum 0.0  0.0 1.0  1e-6 &
mut_accum 0.75 0.0 0.25 1e-6 &
mut_accum 0.5  0.0 0.5  1e-6 &
mut_accum 0.25 0.0 0.75 1e-6 &

#
# introgression simulations
#

NUM_REPLICATES=20

introgression() {
    slim -d prop_add=$1                       \
	 -d prop_partrec=$2                   \
	 -d prop_rec=$3                       \
	 -d segment_length=$SEGMENT_LENGTH    \
	 -d neutral_spacing=$NEUTRAL_SPACING  \
	 -d pop_size=$4                       \
	 -d recomb_rate=$5                    \
	 -d rep=$6                            \
	 slim/dominance_mix__introgression.slim
}

introgression_reps() {
    for rep in `seq 1 $NUM_REPLICATES`; do
        introgression $1 $2 $3 $4 $5 $rep
    done
}

introgression_reps 1.0  0.0 0.0  10000 1e-8 &
introgression_reps 0.25 0.0 0.75 10000 1e-8 &
introgression_reps 0.5  0.0 0.5  10000 1e-8 &
introgression_reps 0.75 0.0 0.25 10000 1e-8 &
introgression_reps 0.0  0.0 1.0  10000 1e-8 &
introgression_reps 0.0  1.0 0.0  10000 1e-8 &

introgression_reps 1.0  0.0 0.0  10000 1e-7 &
introgression_reps 0.25 0.0 0.75 10000 1e-7 &
introgression_reps 0.5  0.0 0.5  10000 1e-7 &
introgression_reps 0.75 0.0 0.25 10000 1e-7 &
introgression_reps 0.0  0.0 1.0  10000 1e-7 &
introgression_reps 0.0  1.0 0.0  10000 1e-7 &

introgression_reps 1.0  0.0 0.0  10000 1e-6 &
introgression_reps 0.25 0.0 0.75 10000 1e-6 &
introgression_reps 0.5  0.0 0.5  10000 1e-6 &
introgression_reps 0.75 0.0 0.25 10000 1e-6 &
introgression_reps 0.0  0.0 1.0  10000 1e-6 &
introgression_reps 0.0  1.0 0.0  10000 1e-6 &


#
# exonic structure simulations
#
slim -d prop_add=0.25 -d prop_rec=0.75 slim/dominance_mix__mut_accum__exons.slim &
slim -d prop_add=0.5  -d prop_rec=0.5  slim/dominance_mix__mut_accum__exons.slim &
slim -d prop_add=0.75 -d prop_rec=0.25 slim/dominance_mix__mut_accum__exons.slim &
