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

    parser.add_argument('--segment-length',
                        help='Length of the simulated segment')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--recomb-rate', type=float, help='Recombination rate')
    group.add_argument('--recomb-map', metavar='FILE', help='Recombination map file')

    parser.add_argument('--dominance-coef', type=float, required=True,
                        help='Dominance coefficient of deleterious mutations')

    parser.add_argument('--afr-size', type=int, default=10000,
                        help='Effective population size of the ancestral population')
    parser.add_argument('--nea-size', type=int, default=1000,
                        help='Effective population size of the Neanderthal population')
    parser.add_argument('--founder-size', type=int, default=1861,
                        help='Effective population size of the founding population')

    parser.add_argument('--burnin', type=int, default=1100000,
                        help='Length of initial burnin [years]')
    parser.add_argument('--hum-nea-split', type=int, default=400000,
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
    burnin          = years_to_gen(args.burnin)
    hum_nea_split   = years_to_gen(args.hum_nea_split)
    out_of_africa   = burnin + hum_nea_split - years_to_gen(args.out_of_africa)

    # specify recombination rate either as a fixed value or as
    # a pair of coordinates and recombination rates as required
    # by SLiM
    if args.recomb_map:
        recomb_map = pd.read_table(args.recomb_map,
                                   names=['interval_end', 'recomb_rate'])
        ends  = 'c(' + ','.join(str(i) for i in recomb_map.interval_end) + ')'
        rates = 'c(' + ','.join(str(i) for i in recomb_map.recomb_rate) + ')'

        args.recomb_rate = rates + ',\n' + ends
        args.segment_length = max(recomb_map.interval_end)

    if not args.segment_length:
        parser.error('Segment length has to be specified'
                     '(or taken from the recombination map)!')

    # values to fill in the SLiM template file
    mapping = {
        'segment_length' : int(float(args.segment_length)),
        'dominance_coef' : args.dominance_coef,
        'recomb_rate'    : args.recomb_rate,
        'founder_size'   : args.founder_size,
        'afr_size'       : args.afr_size,
        'nea_size'       : args.nea_size,
        'burnin'         : burnin,
        'out_of_africa'  : out_of_africa,
        'output_prefix'  : args.output_prefix
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
