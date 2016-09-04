#!/usr/bin/env python3

import os
import random
import argparse
from string import Template
from tempfile import NamedTemporaryFile
import subprocess

import numpy as np


def years_to_gen(years, years_per_gen=25):
    """Convert years to generations."""
    return int(years / years_per_gen)


def parse_and_save(slim_input, output_file):
    '''Parse the list of neutral allele frequencies simulated by SLiM.'''
    prefix = '#OUT:\t'
    for line in slim_input.split('\n'):
        if line.startswith(prefix):
            print(line[len(prefix):], file=output_file, flush=True)


def count_finished_simulations(path):
    """Parse the table with SLiM desert simulation results and
    return the number of simulation results it already contains.
    """
    with open(path, 'r') as output_file:
        return sum(1 for _ in output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run SLiM simulation of Neanderthal deserts')

    parser.add_argument('--num-iters', type=int, required=True,
                        help='Number of iterations of the simulation')

    parser.add_argument('--segment-length', required=True,
                        help='Length of the simulated segment')
    parser.add_argument('--spacing', type=int, default=10000,
                        help='Number of bases between neutral markers')
    parser.add_argument('--recomb-rate', type=float, default=1e-8,
                        help='Recombination rate')

    parser.add_argument('--dominance-coef', type=float, required=True,
                        help='Dominance coefficient of deleterious mutations')

    parser.add_argument('--admixture-rate', type=float, default=0.1,
                        help='Neanderthal migration rate')
    parser.add_argument('--founder-size', type=int, default=1861,
                        help='Effective population size of the founding population')
    parser.add_argument('--afr-size', type=int, default=10000,
                        help='Effective population size of the African population')
    parser.add_argument('--nea-size', type=int, default=1000,
                        help='Effective population size of the Neanderthal population')

    parser.add_argument('--burnin', type=int, default=1000000,
                        help='Length of initial burnin [years]')
    parser.add_argument('--hum-nea-split', type=int, default=500000,
                        help='Split time of Neanderthals [years ago]')
    parser.add_argument('--out-of-africa', type=int, default=55000,
                        help='Out of Africa migration [years ago]')
    parser.add_argument('--admixture-start', type=int, default=55000,
                        help='Start of Neanderthal admixture [years ago]')
    parser.add_argument('--admixture-end', type=int,
                        help='End of Neanderthal admixture [years ago]')
    parser.add_argument('--eur-growth', type=int, default=23000,
                        help='Start of European growth after the split with Asians [years ago]')

    parser.add_argument('--output-file', required=True,
                        help='Where to save the output table with frequencies')

    args = parser.parse_args()

    slim_template = Template(open('neanderthal_admixture.slim', 'r').read())

    args.segment_length = int(float(args.segment_length))

    # place neutral mutations at regular intervals
    neutral_pos = [pos for pos in range(int(args.spacing / 2),
                                        args.segment_length,
                                        args.spacing)]

    # convert arguments specified in years BP to generations since the
    # start of the simulation
    burnin          = years_to_gen(args.burnin)
    hum_nea_split   = years_to_gen(args.hum_nea_split) 
    out_of_africa   = burnin + hum_nea_split - years_to_gen(args.out_of_africa)
    admixture_start = burnin + hum_nea_split - years_to_gen(args.admixture_start)
    if args.admixture_end:
        admixture_end   = burnin + hum_nea_split - years_to_gen(args.admixture_end)
    else:
        admixture_end = admixture_start + 1
    eur_growth      = burnin + hum_nea_split - years_to_gen(args.eur_growth)

    # values to fill in the SLiM template file
    mapping = {
        'segment_length' : args.segment_length,
        'spacing' : args.spacing,
        'neutral_pos' : 'c(' + ','.join(str(pos) for pos in neutral_pos) + ')',
        'neutral_count': len(neutral_pos),
        'dominance_coef' : args.dominance_coef,
        'recomb_rate' : args.recomb_rate,
        'founder_size' : args.founder_size,
        'afr_size' :  args.afr_size,
        'nea_size' : args.nea_size,
        'admixture_rate' : args.admixture_rate,
        'burnin' : burnin,
        'out_of_africa' : out_of_africa,
        'place_neutral' : admixture_start - 1,
        'admixture_start': admixture_start,
        'admixture_end' : admixture_end,
        'eur_growth' : eur_growth,
        'sim_length' : burnin + hum_nea_split,
    }
                
    # test if a given output file already contains some simulation
    # results and if it does, count how many
    already_finished = 0
    if os.path.isfile(args.output_file):
        already_finished = count_finished_simulations(args.output_file)

    with open(args.output_file, 'a') as output_file:
        # run the given number of SLiM iterations
        for i in range(1, args.num_iters - already_finished + 1):
            # fill in the SLiM template with simulation parameter values and
            # run SLiM with it as an input file
            with NamedTemporaryFile('w') as slim_file:
                print(slim_template.substitute(mapping),
                      file=slim_file, flush=True)

                slim_output = subprocess.run(['slim', slim_file.name],
                                             stdout=subprocess.PIPE,
                                             universal_newlines=True)
                parse_and_save(slim_output.stdout, output_file)
