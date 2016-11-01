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
    group.add_argument('--exonic-sites', metavar='FILE',
                        help='Positions of exonic sites from the archaic admixture array')
    group.add_argument('--neutral-spacing', type=int, metavar='N',
                        help='Place evenly distributed neutral sites every N basepairs')
    parser.add_argument('--nonexonic-sites', metavar='FILE',
                        help='Positions of non-exonic sites from the archaic admixture array')

    parser.add_argument('--dominance-coef', type=float, required=True,
                        help='Dominance coefficient of deleterious mutations')

    parser.add_argument('--admixture-rate', type=float, default=0.1,
                        help='Neanderthal migration rate')
    parser.add_argument('--admixture-time', type=int, default=55000,
                        help='Time of Neanderthal admixture [years ago]')

    parser.add_argument('--out-of-africa', type=int, default=55000,
                        help='Out of Africa migration [years ago] (start of the simulation)')
    parser.add_argument('--eur-growth', type=int, default=23000,
                        help='Start of European growth after the split with Asians [years ago]')

    parser.add_argument('--sampling-times', nargs='*', type=int, default=[],
                        help='List of timepoints (in years BP) at which to sample'
                        ' Neanderthal ancestry in a population')

    parser.add_argument('--save-nea-mutations', action='store_true',
                        help='Save the data about introgressed Neanderthal mutations'
                        ' in the founder population')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--output-prefix', metavar='FILE', help='Prefix of output files')
    group.add_argument('--dump-slim', metavar='FILE', help='Dump the SLiM config'
                        ' file without running the simulation')

    args = parser.parse_args()

    slim_template = Template(open('slim/introgression.slim', 'r').read())

    # convert arguments specified in years BP to generations since the
    # start of the simulation
    out_of_africa   = years_to_gen(args.out_of_africa) + 1
    admixture_time  = out_of_africa - years_to_gen(args.admixture_time) + 1
    eur_growth      = out_of_africa - years_to_gen(args.eur_growth)

    # load the SLiM 0-based coordinates of exons
    exon_coords = pd.read_table(args.exon_coordinates, sep='\t')
    genomic_elements = '\n'.join('initializeGenomicElement(g1, {}, {});'.format(s, e)
                                 for s, e in zip(exon_coords.slim_start,
                                                 exon_coords.slim_end))

    # convert sampling times from years BP to generations since
    # the start of the simulation
    sampling_times = [out_of_africa - years_to_gen(s) for s in args.sampling_times]

    # load the SLiM 0-based coordinates of recombination gaps
    recomb_map = pd.read_table(args.recomb_map)

    if args.exonic_sites:
        # read coordinates of sites from the archaic admixture array
        exonic_sites_coords = pd.read_table(args.exonic_sites, names=['slim_start'])
        if args.nonexonic_sites:
            nonexonic_sites_coords = pd.read_table(args.nonexonic_sites, names=['slim_start'])
    else:
        exonic_sites_coords = pd.DataFrame(
            {'slim_start': [pos for pos in range(0,
                                                 max(recomb_map.slim_end),
                                                 args.neutral_spacing)]
            })
        nonexonic_sites_coords = pd.DataFrame({'slim_start': []})

    # values to fill in the SLiM template file
    mapping = {
        'population_file' : args.population_file,
        'recomb_ends'      : 'c(' + ','.join(str(i) for i in recomb_map.slim_end) + ')',
        'recomb_rates'     : 'c(' + ','.join(str(i) for i in recomb_map.recomb_rate) + ')',
        'genomic_elements' : genomic_elements,
        'exonic_pos'       : 'c(' + ','.join(str(pos) for pos in exonic_sites_coords.slim_start) + ')',
        'nonexonic_pos'    : 'c(' + ','.join(str(pos) for pos in nonexonic_sites_coords.slim_start) + ')',
        'dominance_coef'  : args.dominance_coef,
        'admixture_rate'  : args.admixture_rate,
        'out_of_africa'   : out_of_africa,
        'prior_admixture' : admixture_time - 1,
        'admixture_time'  : admixture_time,
        'admixture_end'   : admixture_time + 1,
        'eur_growth'      : eur_growth,
        'sim_length'      : out_of_africa,
        'sampling_times'  : 'c(' + ','.join(str(i) for i in sampling_times) + ')',
        'output_prefix'   : args.output_prefix
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
