#!/usr/bin/env python3

import argparse
from string import Template
from tempfile import NamedTemporaryFile
import subprocess
import logging
import random
import glob

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s :: %(levelname)s :: %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger()


def years_to_gen(years, years_per_gen=25):
    """Convert years to generations."""
    return int(years / years_per_gen)


def slim_vector(xs):
    """Convert a list of numbers to a SLiM code creating the same list."""
    return "c(" + ",".join(str(x) for x in xs) + ")"

def chrom_subset(df, chrom):
    """Subset SLiM coordinates in a dataframe to a single chromosome."""
    df = df.query("chrom == '" + chrom + "'").reset_index(drop=True).copy()
    chrom_start = df.slim_start[0]
    df.slim_start = df.slim_start - chrom_start
    df.slim_end = df.slim_end - chrom_start
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run SLiM simulation of Neanderthal deserts")

    parser.add_argument("--population-file", required=True,
                        help="File with simulated AMH and Neanderthal populations")

    parser.add_argument("--regions", metavar="FILE", required=True,
                        help="Table with SLiM-based coordinates of functional regions "
                        "(0-based SLiM start position, 0-based SLiM end position")
    parser.add_argument("--sites", metavar="FILE", required=True,
                        help="Table with 0-based SLiM coordinates of informative sites")
    parser.add_argument("--recomb-map", metavar="FILE", required=True,
                        help="Table  with the SLiM recombination map (0-based SLiM end "
                        "position, recombination rate)")
    parser.add_argument("--chrom", help="Simulate just one chromosome ('chrN')")

    parser.add_argument("--mut-rate", metavar="FILE", type=float, default=0.0,
                        help="Mutation rate in the simulated region")
    parser.add_argument("--dominance-coef", type=float, required=True,
                        help="Dominance coefficient of deleterious mutations")

    group = parser.add_argument_group()
    group.add_argument("--modify-at", type=int, help="At what time to modify selection"
                       " coefficients [generations]")
    group.add_argument("--modify-what", type=str, help="Which mutation type to modify?", default="m0")
    group.add_argument("--modify-fraction", type=float, help="What fraction of selection"
                        " coefficients to modify?")
    group.add_argument("--modify-count", type=int, help="What number of mutations to modify?")
    group.add_argument("--multiply-s", type=float, help="Multiply each selection coefficient")
    group.add_argument("--fix-s", type=float, help="Fix a selection coefficient to single value")
    group.add_argument("--force-neutral", action="store_true", help="Simulate neutrality "
                       "(set all deleterious mutations to s=0)")

    parser.add_argument("--admixture-source", type=str, default="p2",
                        help="SLiM ID of an admixture source population (p2 = Nea, p4 = Den)")
    parser.add_argument("--admixture-rate", type=float, default=0.1,
                        help="Neanderthal migration rate")
    parser.add_argument("--admixture-time", type=int, default=55000,
                        help="Time of the Neanderthal admixture [years ago]")
    parser.add_argument("--admixture-end", type=int, default=55000,
                        help="End of the Neanderthal admixture [years ago]")

    parser.add_argument("--out-of-africa", type=int, default=55000,
                        help="Out of Africa migration [years ago] (start of the simulation)")

    parser.add_argument("--model", type=str, required=True, choices=["constant", "gravel", "linear"],
                        help="Demographic model to use in the simulation.")
    parser.add_argument("--founder-size", type=int)

    parser.add_argument("--dilution", nargs=3, default=None,
                        help="Tuple of length 3, specifying the rate, start and end "
                        "of the dilution from a population with no Neanderthal admixture")

    parser.add_argument("-sampling-times", nargs="*", type=int, default=[],
                        help="List of timepoints (in years BP) at which to save the simulation state")

    parser.add_argument("--terminate-after", type=int, help="How many generations to simulate? "
                        "(stop after '--out-of-africa' years by default)")

    parser.add_argument("--output-prefix", metavar="FILE", help="Prefix of output VCF files")
    parser.add_argument("--vcf-times", nargs="+", type=int, help="Generation times (in generations)"
                        " for VCF output", default=[])
    parser.add_argument("--vcf-sample", type=int, help="How many individuals to sample?")
    parser.add_argument("--dump-slim", metavar="FILE", help="Dump the generated SLiM config "
                        "file without running the simulation")

    args = parser.parse_args()

    if args.modify_fraction and args.modify_count:
        parser.error("Both fraction and count of mutations to modify cannot be specified")

    if args.founder_size and args.model != "constant":
        parser.error("Founder population size can be specified only for the constant"
                     " population size model.")

    # create the SLiM template file for a specified demographic model
    slim_template = Template(
        open("code/slim/introgression.slim", "r").read() +
        open("code/slim/" + args.model + ".slim", "r").read() +
        (open("code/slim/dilution.slim", "r").read() if args.dilution else "")
    )

    # convert arguments specified in years BP to generations since the
    # start of the simulation
    out_of_africa   = years_to_gen(args.out_of_africa) + 1
    admixture_time  = out_of_africa - years_to_gen(args.admixture_time) + 1
    admixture_end   = admixture_time if not args.admixture_end else out_of_africa - years_to_gen(args.admixture_end) + 1

    if args.modify_at and args.modify_at <= admixture_time:
        parser.error("Time of mutation modification must occur at least one generation"
                     " after admixture")

    # set the appropriate growth rate and effective population size of the non-African
    # population after the out of Africa migration
    if args.model == "constant":
        founder_size = 10000 if not args.founder_size else args.founder_size
        exp_growth = -1
    elif args.model == "gravel":
        founder_size = 1861
        exp_growth = out_of_africa - years_to_gen(23000)
    elif args.model == "linear":
        founder_size = 1861
        exp_growth = out_of_africa - years_to_gen(10000)

    # convert sampling times from years BP to generations since
    # the start of the simulation
    sampling_times = [out_of_africa - years_to_gen(t) for t in args.sampling_times]

    # load the SLiM 0-based coordinates of regions to simulate
    region_coords = pd.read_table(args.regions, sep="\t")

    # load the SLiM 0-based coordinates of recombination breaks
    recomb_map = pd.read_table(args.recomb_map)

    # read coordinates of sites from the archaic admixture array
    sites_coords = pd.read_table(args.sites)

    if args.chrom:
        region_coords = chrom_subset(region_coords, args.chrom)
        recomb_map = chrom_subset(recomb_map, args.chrom)
        sites_coords = chrom_subset(sites_coords, args.chrom)

    if args.multiply_s is not None:
        modifier = args.multiply_s
    elif args.fix_s is not None:
        modifier = args.fix_s
    else:
        modifier = "NULL"

    if args.terminate_after:
        terminate_after = args.terminate_after + 1
    else:
        end = sampling_times[-1] + admixture_time if args.sampling_times else 1e9
        terminate_after = min(out_of_africa, end)

    # values to fill in the SLiM template file
    mapping = {
        "population_file"   : args.population_file,
        "recomb_ends"       : slim_vector(recomb_map.slim_end),
        "recomb_rates"      : slim_vector(recomb_map.recomb_rate),
        "genomic_elements"  : "\n".join("initializeGenomicElement(g1, {}, {});".format(s, e)
                                        for s, e in zip(region_coords.slim_start,
                                                        region_coords.slim_end)),
        "mut_rate"          : args.mut_rate,
        "dominance_coef"    : args.dominance_coef,
        "positions"         : slim_vector(sites_coords.slim_start),

        "modify_at"         : args.modify_at if args.modify_at else admixture_time + 1,
        "modify_what"       : args.modify_what,
        "modify_fraction"   : args.modify_fraction if args.modify_fraction is not None else "F",
        "modify_count"      : args.modify_count if args.modify_count is not None else "F",
        "multiply_s"        : "T" if args.multiply_s is not None else "F",
        "fix_s"             : "T" if args.fix_s is not None else "F",
        "modifier"          : modifier,
        "neutral_sim"       : "T" if args.force_neutral else "F",

        "founder_size"      : founder_size,
        "admixture_source"  : args.admixture_source,
        "admixture_rate"    : args.admixture_rate,
        "admixture_time"    : admixture_time,
        "admixture_end"     : admixture_end,
        "before_admixture"  : admixture_time - 1,
        "after_admixture"   : admixture_end + 1,
        "exp_growth"        : exp_growth,
        "sim_length"        : out_of_africa,
        "terminate_after"   : terminate_after,
        "sampling_times"    : slim_vector(sampling_times),
        "sim_length"        : out_of_africa,
        "simulate_dilution" : "T" if args.dilution else "F",
        "output_prefix"     : args.output_prefix,
        "vcf_times"         : slim_vector(args.vcf_times),
        "vcf_sample"        : args.vcf_sample
    }

    # if simulating dilution, add the necessary parameters
    if args.dilution:
        mapping.update({
            "dilution_rate"  : float(args.dilution[0]),
            "dilution_start" : out_of_africa - years_to_gen(int(args.dilution[1])) + 1,
            "dilution_end"   : out_of_africa - years_to_gen(int(args.dilution[2]))
        })

    if args.dump_slim:
        with open(args.dump_slim, "w") as slim_file:
            print(slim_template.substitute(mapping), file=slim_file)
    else:
        # fill in the SLiM template with simulation parameter values and
        # use it as an input for SLiM
        with NamedTemporaryFile("w") as slim_file:
            print(slim_template.substitute(mapping), file=slim_file, flush=True)
            logger.info("Running simulation from SLiM input file '{}'".format(slim_file.name))

            slim_output = subprocess.run(["slim", "-s", str(random.randint(1, 10**13)),
                                          slim_file.name],
                                         universal_newlines=True)

            logger.info("Simulation from SLiM input file '{}' done (return code = {})".format(slim_file.name, slim_output.returncode))
