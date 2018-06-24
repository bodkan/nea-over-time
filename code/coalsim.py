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
    sample_names = snps.columns

    yri = list(sample_names[sample_names.str.startswith("yri")])
    nea = list(sample_names[sample_names.str.startswith("nea")])
    dinka = list(sample_names[sample_names.str.startswith("dinka")])

    afr_freq = calc_freq(snps, yri + dinka)
    nea_freq = calc_freq(snps, nea)

    return snps.loc[((afr_freq == 0) & (nea_freq == 1)) | ((afr_freq == 1) & (nea_freq == 0))]


def generate_names(pop_params):
    sample_names = []
    for p in pop_params:
        n_pop = len(pop_params[p]["t_sample"])
        sample_names.extend([f"{p}{i}" for i in range(n_pop)])
    return sample_names


def define_samples(pop_params):
    sample_names = []
    for i, pop in enumerate(pop_params):
        times = [years_to_gen(t) for t in pop_params[pop]["t_sample"]]
        sample_names.extend([msp.Sample(population=i, time=t) for t in times])
    return sample_names


def define_popconfig(pop_params):
    return [msp.PopulationConfiguration(
        initial_size=pop_params[p]["Ne"]) for p in pop_params]


def sample_ages(emh_ages_path):
    emh_ages = pd.read_csv(emh_ages_path, sep=" ", names=["name", "age"]).age
    emh = pd.DataFrame({"name": [f"eur{i}" for i in range(len(emh_ages))],
                        "age": emh_ages})

    n_eur = 20
    eur = pd.DataFrame({
        "name": [f"eur{i}" for i in range(len(emh), len(emh) + n_eur)],
        "age": list(repeat(0, n_eur))
    })

    samples = pd.concat([emh, eur]).reset_index(drop=True)
    samples["post_admixture"] = 55000 - samples.age

    return samples


def get_sfs(freq, n_bins=100):
    freq_bins = pd.cut(freq,
                       bins=[i / n_bins for i in range(n_bins + 1)],
                       labels=range(1, n_bins + 1),
                       include_lowest=True) \
        .value_counts()
    return pd.DataFrame({"bin": freq_bins.index,
                         "count": freq_bins,
                         "prop": freq_bins / sum(freq_bins)}).sort_values("bin")


def calc_freq(snps, sample_names=None):
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
    CHIMP, YRI, DIN, NEA, EUR = [pop_params[p]["id"] for p in pop_params]

    # population split times
    t_split_eur, t_split_dinka, t_split_nea, t_split_ch = \
        [years_to_gen(pop_params[p]["t_split"])
         for p in ["eur", "dinka", "nea", "chimp"]]

    t_admix = years_to_gen(admix_params["t_admix"])

    demography = [
        # EUR-AFR gene-flow
        msp.MigrationRateChange(time=0, rate=migr_params["eur_to_afr"],
                                matrix_index=(DIN, EUR)),
        msp.MigrationRateChange(time=0, rate=2 * migr_params["eur_to_afr"],
                                matrix_index=(YRI, EUR)),
        msp.MigrationRateChange(time=0, rate=migr_params["afr_to_eur"],
                                matrix_index=(EUR, DIN)),

        # end of EUR-AFR gene-flow backwards in time
        msp.MigrationRateChange(time=years_to_gen(migr_params["t"]), rate=0.0,
                                matrix_index=(DIN, EUR)),
        msp.MigrationRateChange(time=years_to_gen(migr_params["t"]), rate=0.0,
                                matrix_index=(YRI, EUR)),
        msp.MigrationRateChange(time=years_to_gen(migr_params["t"]), rate=0.0,
                                matrix_index=(EUR, DIN)),

        # Neanderthal admixture
        msp.MassMigration(time=t_admix, source=EUR, dest=NEA,
                          proportion=admix_params["rate"]),

        # population size during the bottleneck
        msp.PopulationParametersChange(
            time=t_admix,
            initial_size=pop_params["eur"]["bottle_Ne"],
            population_id=EUR),

        # out of Africa migration
        msp.MassMigration(time=t_split_eur, source=EUR, destination=DIN,
                          proportion=1.0),

        # YRI-DIN split
        msp.MassMigration(time=t_split_dinka, source=DIN, destination=YRI,
                          proportion=1.0),

        # Neanderthal split
        msp.MassMigration(time=t_split_nea, source=NEA, destination=YRI,
                          proportion=1.0),

        # chimpanzee split
        msp.MassMigration(time=t_split_ch, source=YRI, destination=CHIMP,
                          proportion=1.0),
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


def assign_times(ts):
    """Randomly assign time of origin to each mutation."""
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
    """Assign population of origin to each mutation."""
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


def get_pop(idx, pops):
    return pops[pops.idx == idx]["pop"].values[0]


def get_true_snps(ts, all_snps, pop_params):
    n_nea = len(pop_params["nea"]["t_sample"])
    nea_names = [f"nea{i}" for i in range(n_nea)]

    mut_pops = assign_pops(ts)
    nea_snps = all_snps[mut_pops == pop_params["nea"]["id"]]
    fixed_nea_snps = nea_snps[nea_snps[nea_names].sum(axis=1) == n_nea]

    pops = pd.DataFrame({"pop": mut_pops, "idx": all_snps.index})
    print(pd.Series([get_pop(i, pops) for i in fixed_nea_snps.index]).value_counts())

    return fixed_nea_snps


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--t", help="Time of gene-flow", required=True)
    parser.add_argument("--eur-to-afr", help="EUR -> AFR migration rate", required=True)
    parser.add_argument("--afr-to-eur", help="AFR -> EUR migration rate", required=True)
    parser.add_argument("--afr-to-eur", help="AFR -> EUR migration rate", required=True)
    parser.add_argument("--hap-length", help="Length of simulated haplotypes", required=True)
    parser.add_argument("--hap-num", help="Number of simulated haplotypes", required=True)
    parser.parse_args()

    samples = sample_ages("data/emh_ages.txt")

    admix_params = {
        "rate": 0.015,
        "t_admix": 55000
    }

    pop_params = {
        "chimp": {"id": 0, "Ne": 10000, "t_sample": 1 * [0], "t_split": 6_000_000},
        "yri": {"id": 1, "Ne": 10000, "t_sample": 5 * [0]},
        "dinka": {"id": 2, "Ne": 10000, "t_sample": 5 * [0], "t_split": 150_000},
        "nea": {"id": 3, "Ne": 1000, "t_sample": 4 * [70000], "t_split": 500_000},
        "eur": {"id": 4, "Ne": 10000, "t_sample": samples.age, "t_split": 60_000,
                "bottle_Ne": 2000}
    }

    migr_params = {"t": 0, "afr_to_eur": 0, "eur_to_afr": 0}

    ts = run_sim(admix_params, pop_params, migr_params, seq_len=100_000_000)

    all_snps = get_all_snps(ts, generate_names(pop_params))
    fix_snps = get_nea_snps(all_snps)
    true_snps = get_true_snps(ts, all_snps, pop_params)

    pd.Series((fix_snps[s] == fix_snps.nea0).mean() for s in samples.name).mean()
    pd.Series((true_snps[s] == true_snps.nea0).mean() for s in samples.name).mean()







#
#
#
#
# fix = set(fix_snps.index)
# true = set(true_snps.index)
#
# mut_pops = assign_pops(ts)
#
#
#
#
# fix_pops = pd.Series([get_pop(i) for i in fix])
# true_pops = pd.Series([get_pop(i) for i in true])
#
# # pd.DataFrame((m.source, m.dest) for m in ts.migrations()).drop_duplicates()
#
#
#
#

