#!/usr/bin/env python3

import argparse
from string import Template
from tempfile import NamedTemporaryFile
import subprocess
import logging
import random

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

def generate_elements(regions, id, h, mean, shape):
    mut_type = "initializeMutationType(\"m{}\", {}, \"g\", {}, {});\n\n".format(id, h, mean, shape)
    elem_type = "initializeGenomicElementType(\"g{}\", m{}, 1.0);\n\n".format(id, id)
    elements = "\n".join("initializeGenomicElement(g{}, {}, {});".format(id, s, e)
                        for s, e in zip(regions.slim_start, regions.slim_end))
    return mut_type + elem_type + elements

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate mutations in AMH and Neanderthals")

    parser.add_argument("--regions", metavar="FILE", required=True,
                        help="Table with SLiM-based coordinates of functional regions "
                             "(0-based SLiM start position, 0-based SLiM end position")
    parser.add_argument("--sites", metavar="FILE", required=True,
                        help="Table with 0-based SLiM coordinates of informative sites")
    parser.add_argument("--recomb-map", metavar="FILE", required=True,
                        help="Table  with the SLiM recombination map (0-based SLiM end "
                             "position, recombination rate)")
    parser.add_argument("--chrom", help="Simulate just one chromosome ('chrN')")

    parser.add_argument("--mut-rate", metavar="FILE", type=float, required=True,
                        help="Mutation rate in the simulated region")
    parser.add_argument("--dominance-coef", type=float, default=0.5,
                        help="Dominance coefficient of deleterious mutations")

    parser.add_argument("--anc-size", type=int, default=10000,
                        help="Effective population size of the ancestral population")
    parser.add_argument("--nea-size", type=int, default=1000,
                        help="Effective population size of the Neanderthal population")
    parser.add_argument("--founder-size", type=int, default=1861,
                        help="Effective population size of the founding non-African population")

    parser.add_argument("--hum-nea-split", type=int, default=500000,
                        help="Split time between modern humans and Neanderthals [years ago]")
    parser.add_argument("--nea-den-split", type=int, default=None,
                        help="Split time between modern humans and Neanderthals [years ago]")
    parser.add_argument("--out-of-africa", type=int, default=55000,
                        help="Out of Africa migration [years ago]")

    parser.add_argument("--output", required=True, help="Where to save the simulation output file")

    parser.add_argument("--dump-slim", metavar="FILE", help="Dump the generated SLiM config"
                        " file without running the simulation")

    args = parser.parse_args()

    slim_template = Template(
        open("code/slim/mut_accum.slim", "r").read() +
        (open("code/slim/denisovans.slim", "r").read() if args.nea_den_split else "")
    )

    # convert arguments specified in years BP to generations since the
    # start of the simulation
    burnin          = 7 * args.anc_size  # SLiM 1.8 manual recommends 5xNe
    hum_nea_split   = years_to_gen(args.hum_nea_split)
    out_of_africa   = burnin + hum_nea_split - years_to_gen(args.out_of_africa)
    if args.nea_den_split:
        nea_den_split = burnin + hum_nea_split - years_to_gen(args.nea_den_split)

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
    
    # generate only one class of mutation types and genomic elements
    if not any(c.startswith("bin_") for c in region_coords.columns):
        genomic_elements = generate_elements(regions = region_coords,
                                             id = 0, h = args.dominance_coef,
                                             mean=-0.043, shape=0.23)
    else: # generate several different classes of mutations/genomic elements
        # extract the variable name (h_scale or s_scale)
        var = region_coords.filter(like="bin_", axis="columns").columns[0]
        # get the list of multipliers
        multipliers = list(sorted(region_coords[var].unique()))
        n = len(multipliers)
        dom, sel = (multipliers, [1] * n) if var == "bin_h" else ([args.dominance_coef] * n, multipliers)
        # generate SLiM code for generation mutations and genomic elements
        genomic_elements = "\n\n".join(
            generate_elements(regions = region_coords.query("{} == {}".format(var, h if var == "bin_h" else s)),
                              id=10 + i, h=h, mean=-0.043 * s, shape=0.23)
            for i, h, s in zip(range(n), dom, sel)
        )

    # values to fill in the SLiM template file
    mapping = {
        "recomb_ends"      : slim_vector(recomb_map.slim_end),
        "recomb_rates"     : slim_vector(recomb_map.recomb_rate),
        "genomic_elements" : genomic_elements,
        "mut_rate"         : args.mut_rate,
        "dominance_coef"   : args.dominance_coef,
        "positions"        : slim_vector(sites_coords.slim_start),
        "founder_size"     : args.founder_size,
        "anc_size"         : args.anc_size,
        "nea_size"         : args.nea_size,
        "burnin"           : burnin,
        "nea_den_split"    : nea_den_split if args.nea_den_split else -1,
        "out_of_africa"    : out_of_africa,
        "output"           : args.output
    }

    # fill in the SLiM template with simulation parameter values
    slim_file = open(args.dump_slim, "w") if args.dump_slim else NamedTemporaryFile("w")
    print(slim_template.substitute(mapping), file=slim_file, flush=True)

    # run the simulation if not in the debugging mode
    if not args.dump_slim:
        logger.info("Simulating populations (SLiM input {})".format(slim_file.name))
        slim_output = subprocess.run(["slim", "-s", str(random.randint(1, 10**13)), slim_file.name])
        logger.info("Simulation using '{}' done (returned {})".format(slim_file.name, slim_output.returncode))

    slim_file.close()
