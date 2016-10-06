#!/usr/bin/env python3

import argparse
from string import Template
from tempfile import NamedTemporaryFile
import subprocess
import logging
import random

import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s :: %(levelname)s :: %(message)s',
                    datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger()


def years_to_gen(years, years_per_gen=25):
    """Convert years to generations."""
    return int(years / years_per_gen)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run SLiM simulation of Neanderthal deserts')

    parser.add_argument('--population-file', required=True,
                        help='File with simulated AMH and Neanderthal populations')

    parser.add_argument('--exon-coordinates', metavar='FILE', required=True,
                        help='Tab delimited text file with exon coordinates')
    parser.add_argument('--recomb-map', metavar='FILE', required=True,
                        help='Tab delimited text file with the recombination map')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--array-sites', metavar='FILE',
                        help='Positions of sites from the archaic admixture array')
    group.add_argument('--neutral-spacing', type=int, metavar='N',
                        help='Place evenly distributed neutral sites every N basepairs')

    parser.add_argument('--dominance-coef', type=float, required=True,
                        help='Dominance coefficient of deleterious mutations')

    parser.add_argument('--admixture-rate', type=float, default=0.1,
                        help='Neanderthal migration rate')

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
                       ' Neanderthal ancestries over time')
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

    if args.array_sites:
        # read coordinates of sites from the archaic admixture array
        sites_coords = pd.read_table(args.array_sites, names=['slim_start'])
    else:
        # place neutral mutations at regular interval
        sites_coords = pd.DataFrame({'slim_start': [pos for pos in range(int(args.neutral_spacing / 2),
                                                                         max(recomb_map.slim_end),
                                                                         args.neutral_spacing)]})

    # values to fill in the SLiM template file
    mapping = {
        'population_file' : args.population_file,
        'recomb_ends'      : 'c(' + ','.join(str(i) for i in recomb_map.slim_end) + ')',
        'recomb_rates'     : 'c(' + ','.join(str(i) for i in recomb_map.recomb_rate) + ')',
        'genomic_elements' : genomic_elements,
        'neutral_pos'     : 'c(' + ','.join(str(pos) for pos in sites_coords.slim_start) + ')',
        'neutral_count'   : len(sites_coords),
        'dominance_coef'  : args.dominance_coef,
        'admixture_rate'  : args.admixture_rate,
        'out_of_africa'   : out_of_africa,
        'admixture_start' : admixture_start,
        'admixture_end'   : admixture_end,
        'eur_growth'      : eur_growth,
        'sim_length'      : out_of_africa,
        'output_file'     : args.output_file
    }

    if args.dump_slim:
        with open(args.dump_slim, 'w') as slim_file:
            print(slim_template.substitute(mapping),
                  file=slim_file)
    else:
        # fill in the SLiM template with simulation parameter values and
        # use it as an input for SLiM
        with NamedTemporaryFile('w') as slim_file:
            print(slim_template.substitute(mapping),
                  file=slim_file, flush=True)
            logger.info('Running simulation from SLiM input file "{}"'.format(slim_file.name))

            slim_output = subprocess.run(['slim', '-s', str(random.randint(1, 10**13)),
                                          slim_file.name],
                                         universal_newlines=True)

            logger.info('Simulation from SLiM input file "{}" done (return code = {})'.format(slim_file.name, slim_output.returncode))
