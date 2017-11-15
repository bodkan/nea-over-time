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
    parser.add_argument('--mut-rate', metavar='FILE', default=7e-9,
                        help='Mutation rate')

    group = parser.add_argument_group()
    group.add_argument('--modify-at', type=int, help='At what time to modify selection'
                       ' coefficients [generations]')
    group.add_argument('--modify-what', type=str, help='Which mutation type to modify?', default='m0')
    group.add_argument('--modify-fraction', type=float, help='What fraction of selection'
                        ' coefficients to modify?')
    group.add_argument('--modify-count', type=int, help='What number of mutations to modify?')
    group.add_argument('--multiply-s', type=float, help='Multiply each selection coefficient')
    group.add_argument('--fix-s', type=float, help='Fix a selection coefficient to single value')

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
                        help='Time of the Neanderthal admixture [years ago]')
    parser.add_argument('--admixture-end', type=int, default=55000,
                        help='End of the Neanderthal admixture [years ago]')

    parser.add_argument('--out-of-africa', type=int, default=55000,
                        help='Out of Africa migration [years ago] (start of the simulation)')

    parser.add_argument('--model', type=str, required=True, choices=['constant', 'gravel', 'linear'],
                        help='Demographic model to use in the simulation.')
    parser.add_argument('--founder-size', type=int)

    parser.add_argument('--dilution', nargs=3, default=None,
                        help='Tuple of length 3, specifying the rate, start and end'
                        ' of the dilution from a population with no Neanderthal admixture')

    parser.add_argument('--sampling-times', nargs='*', type=int, default=[],
                        help='List of timepoints (in years BP) at which to sample'
                        ' Neanderthal ancestry in a population')

    parser.add_argument('--terminate-after', type=int, help='How many generations to simulate?'
                        ' (stop after "--out-of-africa" years by default)')

    parser.add_argument('--save-trajectories', action='store_true',
                        help='Save table with Neanderthal ancestry statistics')
    parser.add_argument('--save-mutations', action='store_true',
                        help='Save the data about deleterious mutations over'
                        ' the course of the simulation')

    parser.add_argument('--output-prefix', metavar='FILE', help='Prefix of output files')
    parser.add_argument('--dump-slim', metavar='FILE', help='Dump the SLiM config'
                        ' file without running the simulation')

    args = parser.parse_args()

    if not args.save_trajectories and not args.save_mutations:
        parser.error('At least one output option must be chosen (either saving'
                     ' Nea. trajectories or whole SLiM output files).')

    if args.modify_fraction and args.modify_count:
        parser.error('Both fraction and count of mutations to modify cannot be specified')

    if args.founder_size and args.model != 'constant':
        parser.error('Founder population size can be specified only for the constant'
                     ' population size model.')

    # create the SLiM template file for a specified demographic model
    slim_template = Template(
        open('slim/introgression.slim', 'r').read() +
        open('slim/' + args.model + '.slim', 'r').read() +
        (open('slim/dilution.slim', 'r').read() if args.dilution else '')
    )

    # convert arguments specified in years BP to generations since the
    # start of the simulation
    out_of_africa   = years_to_gen(args.out_of_africa) + 1
    admixture_time  = out_of_africa - years_to_gen(args.admixture_time) + 1
    admixture_end   = admixture_time if not args.admixture_end else out_of_africa - years_to_gen(args.admixture_end) + 1

    if args.modify_at and args.modify_at <= admixture_time:
        parser.error('Time of mutation modification must occur at least one generation'
                     ' after admixture')

    # set the appropriate growth rate and effective population size of the non-African
    # population after the out of Africa migration
    if args.model == 'constant':
        founder_size = 10000 if not args.founder_size else args.founder_size
        exp_growth = -1
    elif args.model == 'gravel':
        founder_size = 1861
        exp_growth = out_of_africa - years_to_gen(23000)
    elif args.model == 'linear':
        founder_size = 1861
        exp_growth = out_of_africa - years_to_gen(10000)

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
        'mut_rate'         : float(args.mut_rate),
        'modify_at'        : args.modify_at if args.modify_at else admixture_time + 1,
        'modify_what'      : args.modify_what,
        'modify_fraction'  : args.modify_fraction if args.modify_fraction is not None else 'F',
        'modify_count'     : args.modify_count if args.modify_count is not None else 'F',
        'multiply_s'       : args.multiply_s if args.multiply_s is not None else 'F',
        'fix_s'            : args.fix_s if args.fix_s is not None else 'F',
        'genomic_elements' : genomic_elements,
        'exonic_pos'       : 'c(' + ','.join(str(pos) for pos in exonic_sites_coords.slim_start) + ')',
        'nonexonic_pos'    : 'c(' + ','.join(str(pos) for pos in nonexonic_sites_coords.slim_start) + ')',
        'dominance_coef'  : args.dominance_coef,
        'founder_size'    : founder_size,
        'admixture_rate'  : args.admixture_rate,
        'admixture_time'  : admixture_time,
        'admixture_end'   : admixture_end,
        'before_admixture' : admixture_time - 1,
        'after_admixture' : admixture_end + 1,
        'exp_growth'      : exp_growth,
        'sim_length'      : out_of_africa,
        'terminate_after' : args.terminate_after + 1 if args.terminate_after else min(out_of_africa, (sampling_times[-1] + admixture_time) if args.sampling_times else 1e9),
        'sampling_times'  : 'c(' + ','.join(str(i) for i in sampling_times) + ')',
        'sim_length'        : out_of_africa,
        'simulate_dilution' : 'T' if args.dilution else 'F',
        'save_trajectories' : 'T' if args.save_trajectories else 'F',
        'save_mutations'  : 'T' if args.save_mutations else 'F',
        'output_prefix'   : args.output_prefix
    }

    # if simulating dilution, add the necessary parameters
    if args.dilution:
        mapping.update({
            'dilution_rate'  : float(args.dilution[0]),
            'dilution_start' : out_of_africa - years_to_gen(int(args.dilution[1])) + 1,
            'dilution_end'   : out_of_africa - years_to_gen(int(args.dilution[2]))
        })

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
