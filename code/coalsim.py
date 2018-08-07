#!/usr/bin/env python3

import argparse
import random
from itertools import repeat
from collections import defaultdict

import msprime as msp
import pandas as pd
import numpy as np


def years_to_gen(years, gen_time=25):
    return int(years / gen_time)


def get_all_snps(ts, samples):
    """Extract all simulated SNPs from the tree sequence replicates."""
    if isinstance(ts, msp.TreeSequence):
        ts = [ts]

    return pd.concat(
        pd.DataFrame(
            t.genotype_matrix(),
            columns=samples,
            index=(v.position for v in t.variants())
        ) for t in ts
    )


def get_nea_snps(snps):
    """Filter all SNPs to fixed differences between Africans and Neandertals."""
    sample_names = snps.columns

    yoruba = list(sample_names[sample_names.str.startswith("yoruba")])
    nea = list(sample_names[sample_names.str.startswith("nea")])

    afr_freq = calc_freq(snps, yoruba)
    nea_freq = calc_freq(snps, nea)

    return snps.loc[((afr_freq == 0) & (nea_freq == 1)) | ((afr_freq == 1) & (nea_freq == 0))]


def get_true_snps(ts, all_snps, pop_params):
    """Examine the coalescent trees and find those SNPs that truly originated in
    the Neandertals."""
    n_nea = len(pop_params["nea"]["t_sample"])
    nea_names = [f"nea{i}" for i in range(n_nea)]

    mut_pops = assign_pops(ts)
    nea_snps = all_snps[mut_pops == pop_params["nea"]["id"]]
    fixed_nea_snps = nea_snps[nea_snps[nea_names].sum(axis=1) == n_nea]

    return fixed_nea_snps


def assign_times(ts):
    """Randomly assign time of origin to each mutation.
    Inspired by an example from msprime's documentation."""
    rng = random.Random()
    mut_times = np.zeros(ts.num_mutations)
    for tree in ts.trees():
        for mutation in tree.mutations():
            a = tree.time(mutation.node)
            b = tree.time(tree.parent(mutation.node))
            mut_times[mutation.id] = rng.uniform(a, b)
    return mut_times


def gather_migrations(ts):
    """Gather all migrations in each node in a tree."""
    node_migrations = defaultdict(list)
    for migration in ts.migrations():
        node_migrations[migration.node].append(migration)
    return node_migrations


def assign_pops(ts):
    """Assign population of origin to each mutation.
    Inspired by an example from msprime's documentation."""
    mut_times = assign_times(ts)
    node_migrations = gather_migrations(ts)

    pop_assign = np.repeat(-1, ts.num_mutations)
    for tree in ts.trees():
        for site in tree.sites():
            for mut in site.mutations:
                pop_assign[mut.id] = tree.population(mut.node)
                for mig in node_migrations[mut.node]:
                    if mig.left <= site.position < mig.right:
                        if mig.time < mut_times[mut.id]:
                            assert pop_assign[mut.id] == mig.source
                            pop_assign[mut.id] = mig.dest
    return pop_assign


def pop_samples(pop, pop_params):
    """Generate list of string names of population samples (i.e. ["eur0", "eur1",
    ..., "eurN"] etc.)"""
    names = generate_names(pop_params)
    return [name for name in names if name.startswith(pop)]


def generate_names(pop_params):
    """Generate list of all simulated sample names."""
    sample_names = []
    for p in pop_params:
        n_pop = len(pop_params[p]["t_sample"])
        sample_names.extend([f"{p}{i}" for i in range(n_pop)])
    return sample_names


def define_samples(pop_params):
    """Generate list of sample definitions form msprime."""
    sample_names = []
    for i, pop in enumerate(pop_params):
        times = [years_to_gen(t) for t in pop_params[pop]["t_sample"]]
        sample_names.extend([msp.Sample(population=i, time=t) for t in times])
    return sample_names


def define_popconfig(pop_params):
    """Generate list of population configurations for msprime."""
    return [msp.PopulationConfiguration(
        initial_size=pop_params[p]["Ne"]) for p in pop_params]


def sample_ages(ages):
    """Generate DataFrame with European sample ages."""
    df = pd.DataFrame({
        "name": [f"eur{i}" for i in range(len(ages))],
        "age": ages
    })
    return df


def calc_freq(snps, sample_names=None):
    """Calculate the frequencies of SNPs."""
    if not sample_names:
        sample_names = list(snps.columns)

    if not isinstance(sample_names, list):
        sample_names = [sample_names]

    allele_counts = snps.loc[:, sample_names].sum(axis=1)
    return allele_counts / len(sample_names)


def f4(snps, w, x, y, z):
    w_freq = calc_freq(snps, w)
    x_freq = calc_freq(snps, x)
    y_freq = calc_freq(snps, y)
    z_freq = calc_freq(snps, z)

    return ((w_freq - x_freq) * (y_freq - z_freq)).sum()


def f4_ratio(snps, x, a, b, c, o):
    return f4(snps, a, o, x, c) / f4(snps, a, o, b, c)


def d(snps, w, x, y, z):
    w_freq = calc_freq(snps, w)
    x_freq = calc_freq(snps, x)
    y_freq = calc_freq(snps, y)
    z_freq = calc_freq(snps, z)

    nom = ((w_freq - x_freq) * (y_freq - z_freq)).sum()
    denom = ((w_freq + x_freq - 2 * w_freq * x_freq) * \
             (y_freq + z_freq - 2 * y_freq * z_freq)).sum()

    return nom / denom


def run_sim(admix_params, pop_params, migr_params, *,
            seq_len=10_000, num_replicates=None,
            mut_rate=1e-8, debug=False):
    CHIMP, YORUBA, DIN, NEA, EUR, EAS = [pop_params[p]["id"] for p in pop_params]

    # population split times
    t_split_eur, t_split_eas, t_split_dinka, t_split_nea, t_split_ch = \
        [years_to_gen(pop_params[p]["t_split"])
         for p in ["eur", "eas", "dinka", "nea", "chimp"]]

    t_admix = years_to_gen(admix_params["t_admix"])

    demography = [
        # EUR-AFR gene-flow
        msp.MigrationRateChange(time=0, rate=migr_params["eur_to_dinka"],
                                matrix_index=(DIN, EUR)),
        msp.MigrationRateChange(time=0, rate=migr_params["eur_to_yoruba"],
                                matrix_index=(YORUBA, EUR)),
        msp.MigrationRateChange(time=0, rate=migr_params["dinka_to_eur"],
                                matrix_index=(EUR, DIN)),

        # end of EUR-AFR gene-flow backwards in time
        msp.MigrationRateChange(time=years_to_gen(migr_params["t"]), rate=0.0,
                                matrix_index=(DIN, EUR)),
        msp.MigrationRateChange(time=years_to_gen(migr_params["t"]), rate=0.0,
                                matrix_index=(YORUBA, EUR)),
        msp.MigrationRateChange(time=years_to_gen(migr_params["t"]), rate=0.0,
                                matrix_index=(EUR, DIN)),

        # EUR-EAS split
        msp.MassMigration(time=t_split_eas, source=EAS, destination=EUR, proportion=1.0),

        # Neanderthal admixture
        msp.MassMigration(time=t_admix, source=EUR, dest=NEA, proportion=admix_params["rate"]),

        # population size during the bottleneck
        msp.PopulationParametersChange(
            time=t_admix,
            initial_size=pop_params["eur"]["bottle_Ne"],
            population_id=EUR),

        # out of Africa migration
        msp.MassMigration(time=t_split_eur, source=EUR, destination=DIN, proportion=1.0),

        # YORUBA-DIN split
        msp.MassMigration(time=t_split_dinka, source=DIN, destination=YORUBA, proportion=1.0),

        # Neanderthal split
        msp.MassMigration(time=t_split_nea, source=NEA, destination=YORUBA, proportion=1.0),

        # chimpanzee split
        msp.MassMigration(time=t_split_ch, source=YORUBA, destination=CHIMP, proportion=1.0),
    ]

    samples = define_samples(pop_params)
    pop_config = define_popconfig(pop_params)

    # effective population sizes
    Ne0 = 10000

    if debug:
        msp.DemographyDebugger(
            Ne=Ne0,
            population_configurations=pop_config,
            demographic_events=demography
        ).print_history()
    else:
        return msp.simulate(
            length=seq_len,
            Ne=Ne0,
            mutation_rate=mut_rate,
            recombination_rate=1e-8,
            samples=samples,
            population_configurations=pop_config,
            demographic_events=demography,
            num_replicates=num_replicates,
            record_migrations=True
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--time", help="Start of EUR-AFR gene-flow [years BP]", type=int, required=True)
    parser.add_argument("--eur-to-dinka", help="EUR -> Dinka migration rate", type=float, required=True)
    parser.add_argument("--eur-to-yoruba", help="EUR -> Yoruba migration rate", type=float, required=True)
    parser.add_argument("--dinka-to-eur", help="Dinka -> EUR migration rate", type=float, required=True)
    parser.add_argument("--nea-rate", help="Neandertal admixture rate", type=float, required=True)
    parser.add_argument("--hap-length", help="Length of simulated haplotypes", type=int, required=True)
    parser.add_argument("--hap-num", help="Number of simulated haplotypes", type=int)
    parser.add_argument("--eur-ages", nargs="+", type=int, help="Ages of non-African samples [years BP]", required=True)
    parser.add_argument("--ascertainment", help="Ascertainment scheme", choices=["all", "nea", "yoruba", "dinka", "eur", "eas"],
                        default="all")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--dump-snps", action="store_true", help="Save all SNPs into a table")
    group.add_argument("--calc-stats", action="store_true", help="Calculate the admixture statistics")
    parser.add_argument("--output-file", metavar="FILE", help="Where to save the table of results", required=True)

    # args = parser.parse_args("--output-file asd.txt --dump-snps --hap-length 5_000_000 --eur-ages 0 0 0 --time 0 --output-file asd.txt --eur-to-afr 0 --afr-to-eur 0 --nea-rate 0.05".split())
    args = parser.parse_args()

    # prepare all simulation parameters
    samples = sample_ages(args.eur_ages)

    admix_params = {
        "rate": args.nea_rate,
        "t_admix": 55000
    }

    pop_params = {
        "chimp": {"id": 0, "Ne": 10000, "t_sample": 1 * [0], "t_split": 6_000_000},
        "yoruba": {"id": 1, "Ne": 10000, "t_sample": 50 * [0]},
        "dinka": {"id": 2, "Ne": 10000, "t_sample": 2 * [0], "t_split": 150_000},
        "nea": {"id": 3, "Ne": 1000, "t_sample": 4 * [80000], "t_split": 500_000},
        "eur": {"id": 4, "Ne": 10000, "t_sample": list(samples.age) + 2 * [0], "t_split": 60_000, "bottle_Ne": 2000},
        "eas": {"id": 5, "Ne": 10000, "t_sample": 2 * [0], "t_split": 40_000}
    }

    migr_params = {
        "t": args.time,
        "dinka_to_eur": args.dinka_to_eur,
        "eur_to_yoruba": args.eur_to_yoruba,
        "eur_to_dinka": args.eur_to_dinka
    }

    # simulate the data
    ts = run_sim(admix_params,
                 pop_params,
                 migr_params,
                 seq_len=args.hap_length,
                 num_replicates=None)

    # process the simulations into different tables of SNPs
    all_snps = get_all_snps(ts, generate_names(pop_params))

    if args.dump_snps:
        all_snps.to_csv(args.output_file, sep="\t", index=False)
    else:
        admix_snps = get_nea_snps(all_snps)
        true_snps = get_true_snps(ts, all_snps, pop_params)

        if args.ascertainment != "all":
            last_index = len(pop_params[args.ascertainment]["t_sample"]) - 1
            all_snps = all_snps.query(f"{args.ascertainment}{last_index - 1} != {args.ascertainment}{last_index}")

        # calculate admixture statistics and bind them into a DataFrame
        dinka = ["dinka0", "dinka1"]
        yoruba = ["yoruba0", "yoruba1"]
        altai = ["nea0", "nea1"]
        vindija = ["nea2", "nea3"]

        stats = defaultdict(list)
        for s in samples.name:
            stats["true_prop"].append(true_snps[s].mean())
            stats["admix_prop"].append((admix_snps[s] == admix_snps.nea0).mean())
            stats["direct_f4"].append(f4_ratio(all_snps, s, a=altai, b=vindija, c=yoruba, o="chimp0"))
            stats["indirect_f4"].append(1 - f4_ratio(all_snps, s, a=yoruba, b=dinka, c=altai + vindija, o="chimp0"))
            stats["d"].append(d(all_snps, w="eur0", x=s, y=dinka, z="chimp0") if s != "eur0" else None)
        stats_df = pd.DataFrame(stats)
        stats_df["name"] = samples.name

        final_df = pd.merge(samples, stats_df, on="name")

        final_df.to_csv(args.output_file, sep="\t", index=False)
