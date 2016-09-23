#!/usr/bin/env python3

import os
import argparse
from string import Template
from tempfile import NamedTemporaryFile
import subprocess
import logging

import pandas as pd

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
    parser.add_argument('--recomb-map', metavar='FILE', required=True,
                        help='Tab delimited text file with the recombination map')
    parser.add_argument('--informative-sites', metavar='FILE', required=True,
                        help='Positions of sites from the archaic admixture array')
    parser.add_argument('--exon-coordinates', metavar='FILE', required=True,
                        help='Tab delimited text file with exon coordinates')

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

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--output-file', metavar='FILE', help='Where to save'
                       ' the output table with frequencies at each informative position')
    group.add_argument('--dump-slim', metavar='FILE', help='Dump the SLiM config'
                        ' file without running the simulation')

    args = parser.parse_args()

    slim_template = Template(open('slim/init.slim', 'r').read() +
                             open('slim/admix.slim', 'r').read())

    # convert arguments specified in years BP to generations since the
    # start of the simulation
    out_of_africa   = years_to_gen(args.out_of_africa)
    admixture_start = out_of_africa - years_to_gen(args.admixture_start) + 1
    if args.admixture_length:
        admixture_end = admixture_start + years_to_gen(args.admixture_length)
    else:
        admixture_end = admixture_start + 1
    eur_growth      = out_of_africa - years_to_gen(args.eur_growth)

    # load the SLiM 0-based coordinates of exons
    exon_coords = pd.read_table(args.exon_coordinates, sep='\t')
    genomic_elements = '\n'.join('initializeGenomicElement(g1, {}, {});'.format(s, e)
                                 for s, e in zip(exon_coords.slim_start,
                                                 exon_coords.slim_end))

    # load the SLiM 0-based coordinates of recombination gaps
    recomb_map = pd.read_table(args.recomb_map)

    # read coordinates of sites from the archaic admixture array
    sites_coords = pd.read_table(args.informative_sites, names=['slim_start'])

    # values to fill in the SLiM template file
    mapping = {
        'population_file' : args.population_file,
        'recomb_ends'      : 'c(' + ','.join(str(i) for i in recomb_map.slim_end) + ')',
        'recomb_rates'     : 'c(' + ','.join(str(i) for i in recomb_map.recomb_rate) + ')',
        'genomic_elements' : genomic_elements,
        'neutral_pos'     : 'c(' + ','.join(str(pos) for pos in sites_coords.slim_start) + ')',
        'neutral_count'   : len(sites_coords),
        'dominance_coef'  : args.dominance_coef,
        'founder_size'    : args.founder_size,
        'admixture_rate'  : args.admixture_rate,
        'out_of_africa'   : out_of_africa,
        'admixture_start' : admixture_start,
        'admixture_end'   : admixture_end,
        'eur_growth'      : eur_growth,
        'sim_length'      : out_of_africa,
    }

    if args.dump_slim:
        with open(args.dump_slim, 'w') as slim_file:
            print(slim_template.substitute(mapping),
                  file=slim_file)
    else:
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
                    logger.info('Running simulation #{} (SLiM input file "{}")'.format(i + 1, slim_file.name))

                    slim_output = subprocess.run(['slim', slim_file.name],
                                                 stdout=subprocess.PIPE,
                                                 universal_newlines=True)
                    parse_and_save(slim_output.stdout, output_file)

                    logger.info('Simulation #{} done (return code = {})'.format(i + 1, slim_output.returncode))
