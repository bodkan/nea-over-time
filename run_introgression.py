#!/usr/bin/env python3

import os
import argparse
from string import Template
from tempfile import NamedTemporaryFile
import subprocess
import logging

import pandas

logging.basicConfig(level=logging.INFO, format='%(asctime)s :: %(levelname)s :: %(message)s',
                    datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger()


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
    already_finished = 0

    if os.path.isfile(path):
        with open(path, 'r') as output_file:
            already_finished = sum(1 for _ in output_file)

    return already_finished


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run SLiM simulation of Neanderthal deserts')

    parser.add_argument('--num-iters', type=int, required=True,
                        help='Number of iterations of the simulation')
    parser.add_argument('--population-file', required=True,
                        help='File with simulated AMH and Neanderthal populations')

    parser.add_argument('--segment-length',
                        help='Length of the simulated segment')
    parser.add_argument('--spacing', type=int, default=10000,
                        help='Number of bases between neutral markers')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--recomb-rate', type=float, default=1e-8,
                        help='Recombination rate')
    group.add_argument('--recomb-map', metavar='FILE', help='Recombination map file')

    parser.add_argument('--dominance-coef', type=float, required=True,
                        help='Dominance coefficient of deleterious mutations')

    parser.add_argument('--dominance-coef', type=float, required=True,
                        help='Dominance coefficient of deleterious mutations')

    parser.add_argument('--admixture-rate', type=float, default=0.1,
                        help='Neanderthal migration rate')

    parser.add_argument('--founder-size', type=int, default=1861,
                        help='Effective population size of the founding population')

    parser.add_argument('--out-of-africa', type=int, default=55000,
                        help='Out of Africa migration [years ago] (start of the simulation)')
    parser.add_argument('--admixture-start', type=int, default=55000,
                        help='Start of Neanderthal admixture [years ago]')
    parser.add_argument('--admixture-length', type=int,
                        help='Duration of Neanderthal admixture [years]')
    parser.add_argument('--eur-growth', type=int, default=23000,
                        help='Start of European growth after the split with Asians [years ago]')

    parser.add_argument('--output-file', required=True,
                        help='Where to save the output table with frequencies')

    args = parser.parse_args()

    slim_template = Template(open('slim/init.slim', 'r').read() +
                             open('slim/admix.slim', 'r').read())

    # place neutral mutations at regular intervals
    neutral_pos = [pos for pos in range(int(args.spacing / 2),
                                        args.segment_length,
                                        args.spacing)]

    # convert arguments specified in years BP to generations since the
    # start of the simulation
    out_of_africa   = years_to_gen(args.out_of_africa)
    admixture_start = out_of_africa - years_to_gen(args.admixture_start) + 1
    if args.admixture_length:
        admixture_end = admixture_start + years_to_gen(args.admixture_length)
    else:
        admixture_end = admixture_start + 1
    eur_growth      = out_of_africa - years_to_gen(args.eur_growth)

    # specify recombination rate either as a fixed value or as
    # a pair of coordinates and recombination rates as required
    # by SLiM
    if args.recomb_map:
        recomb_map = pd.read_table(args.recomb_map,
                                   names=['interval_end', 'recomb_rate'])
        ends  = 'c(' + ','.join(str(i) for i in recomb_map.interval_end) + ')'
        rates = 'c(' + ','.join(str(i) for i in recomb_map.recomb_rate) + ')'

        args.recomb_rate = ends + ',\n' + rates
        args.segment_length = max(recomb_map.interval_end)

    if not args.segment_length:
        parser.error('Segment length has to be specified'
                     '(or taken from the recombination map)!')

    # values to fill in the SLiM template file
    mapping = {
        'population_file' : args.population_file,
        'segment_length' : int(float(args.segment_length)),
        'spacing' : args.spacing,
        'neutral_pos' : 'c(' + ','.join(str(pos) for pos in neutral_pos) + ')',
        'neutral_count': len(neutral_pos),
        'dominance_coef' : args.dominance_coef,
        'recomb_rate' : args.recomb_rate,
        'founder_size' : args.founder_size,
        'admixture_rate' : args.admixture_rate,
        'out_of_africa' : out_of_africa,
        'admixture_start': admixture_start,
        'admixture_end' : admixture_end,
        'eur_growth' : eur_growth,
        'sim_length' : out_of_africa,
    }
                
    # count the number of already finished simulations
    already_finished = count_finished_simulations(args.output_file)

    with open(args.output_file, 'a') as output_file:
        # run the given number of SLiM iterations
        for i in range(already_finished, args.num_iters):
            # fill in the SLiM template with simulation parameter values and
            # run SLiM with it as an input file
            with NamedTemporaryFile('w') as slim_file:
                print(slim_template.substitute(mapping),
                      file=slim_file, flush=True)
                logger.info('Simulation #{} (SLiM input file "{}")'.format(i + 1, slim_file.name))

                slim_output = subprocess.run(['slim', slim_file.name],
                                             stdout=subprocess.PIPE,
                                             universal_newlines=True)
                parse_and_save(slim_output.stdout, output_file)
