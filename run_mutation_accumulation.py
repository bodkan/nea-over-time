#!/usr/bin/env python3

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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate mutations in AMH and Neanderthals')

    parser.add_argument('--exon-coordinates', metavar='FILE', required=True,
                        help='Tab delimited text file with exon coordinates')
    parser.add_argument('--recomb-map', metavar='FILE', required=True,
                        help='Tab delimited text file with the recombination map')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--array-sites', metavar='FILE',
                        help='Positions of sites from the archaic admixture array')
    group.add_argument('--even-spacing', type=int, metavar='N',
                        help='Place evenly distributed neutral sites every N basepairs')

    parser.add_argument('--dominance-coef', type=float, required=True,
                        help='Dominance coefficient of deleterious mutations')

    parser.add_argument('--anc-size', type=int, default=10000,
                        help='Effective population size of the ancestral population')
    parser.add_argument('--nea-size', type=int, default=1000,
                        help='Effective population size of the Neanderthal population')
    parser.add_argument('--founder-size', type=int, default=1861,
                        help='Effective population size of the founding population')

    parser.add_argument('--hum-nea-split', type=int, default=500000,
                        help='Split time of Neanderthals [years ago]')
    parser.add_argument('--out-of-africa', type=int, default=55000,
                        help='Out of Africa migration [years ago]')

    parser.add_argument('--output-prefix', required=True,
                        help='Where to save the populations simulated by SLiM')

    parser.add_argument('--dump-slim', metavar='FILE', help='Dump the SLiM config'
                        ' file without running the simulation')

    args = parser.parse_args()

    slim_template = Template(open('slim/init.slim', 'r').read() +
                             open('slim/burnin.slim', 'r').read())

    # convert arguments specified in years BP to generations since the
    # start of the simulation
    burnin          = 5 * args.anc_size  # loosely based on SLiM 1.8 manual
    hum_nea_split   = years_to_gen(args.hum_nea_split)
    out_of_africa   = burnin + hum_nea_split - years_to_gen(args.out_of_africa)

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
        sites_coords = pd.DataFrame({'slim_start': [pos for pos in range(0,
                                                                         max(recomb_map.slim_end),
                                                                         args.even_spacing)]})

    # values to fill in the SLiM template file
    mapping = {
        'recomb_ends'      : 'c(' + ','.join(str(i) for i in recomb_map.slim_end) + ')',
        'recomb_rates'     : 'c(' + ','.join(str(i) for i in recomb_map.recomb_rate) + ')',
        'genomic_elements' : genomic_elements,
        'neutral_pos'      : 'c(' + ','.join(str(pos) for pos in sites_coords.slim_start) + ')',
        'dominance_coef'   : args.dominance_coef,
        'founder_size'     : args.founder_size,
        'anc_size'         : args.anc_size,
        'nea_size'         : args.nea_size,
        'burnin'           : burnin,
        'out_of_africa'    : out_of_africa,
        'output_prefix'    : args.output_prefix
    }

    # fill in the SLiM template with simulation parameter values and
    # run SLiM with it as an input file
    slim_file = open(args.dump_slim, 'w') if args.dump_slim else NamedTemporaryFile('w')

    print(slim_template.substitute(mapping),
          file=slim_file, flush=True)

    # run the simulation if not in the debugging mode
    if not args.dump_slim:
        logger.info('Simulating populations (SLiM input {})'.format(slim_file.name))
        slim_output = subprocess.run(['slim', slim_file.name])
        logger.info('Simulation using "{}" done (returned {})'.format(slim_file.name, slim_output.returncode))

    slim_file.close()
