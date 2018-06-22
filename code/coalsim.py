import random
from itertools import repeat
from collections import defaultdict

import msprime as msp
import pandas as pd
import numpy as np
from intervaltree import IntervalTree, Interval


def years_to_gen(years, gen_time=25):
    return int(years / gen_time)


def nea_snps(snps):
    samples = snps.columns

    afr_freq = calc_freq(snps, list(samples[samples.str.startswith("yri")]) \
                             + list(samples[samples.str.startswith("dinka")]))
    nea_freq = calc_freq(snps, list(samples[samples.str.startswith("nea")]))

    return snps.loc[((afr_freq == 0) & (nea_freq == 1)) | ((afr_freq == 1) & (nea_freq == 0))]


def get_snps(tree_seq, samples):
    if isinstance(tree_seq, msp.TreeSequence):
        tree_seq = [tree_seq]
    
    return pd.concat(
        pd.DataFrame(
            ts.genotype_matrix(),
            columns=samples,
            index = (int(v.position) for v in ts.variants())
        ) for ts in tree_seq
    )


def generate_names(pop_params):
    samples = []
    for p in pop_params:
        samples.extend([p + str(i) for i in range(len(pop_params[p]["t_sample"]))])
    return samples


def define_samples(pop_params):
    samples = []
    for i, pop in enumerate(pop_params):
        times = [years_to_gen(t) for t in pop_params[pop]["t_sample"]]
        samples.extend([msp.Sample(population=i, time=t) for t in times])
    return samples


def define_popconfig(pop_params):
    return [msp.PopulationConfiguration(
        initial_size=pop_params[p]["Ne"]) for p in pop_params]


def sample_ages(emh_ages_path, n_eur=20):
    emh_ages = pd.read_csv(emh_ages_path, sep=" ", names=["name", "age"]).age

    emh = pd.DataFrame({"name" : ["eur{}".format(i) for i in range(len(emh_ages))],
                        "age" : emh_ages})

    eur = pd.DataFrame(
        {"name" : ["eur{}".format(i) for i in range(len(emh), len(emh) + n_eur)],
                   "age" : list(repeat(0, n_eur))})

    samples = pd.concat([emh, eur]).reset_index(drop=True)
    samples["post_admixture"] = 55000 - samples.age

    return samples


def get_sfs(freq, n_bins=100):
    freq_bins = pd.cut(freq,
                       bins = [i / n_bins for i in range(n_bins + 1)],
                       labels=range(1, n_bins + 1),
                       include_lowest=True) \
        .value_counts()
    return pd.DataFrame({"bin": freq_bins.index,
                         "count": freq_bins,
                         "prop": freq_bins / sum(freq_bins)}).sort_values("bin")


def calc_freq(snps, samples=None):
    if not samples:
        samples = list(snps.columns)

    if not isinstance(samples, list):
        samples = [samples]

    allele_counts = snps.loc[:, samples].sum(axis = 1)
    return allele_counts / len(samples)


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

    nom =  ((w_freq - x_freq) * (y_freq - z_freq)).sum()
    denom = ((w_freq + x_freq - 2 * w_freq * x_freq) * \
             (y_freq + z_freq - 2 * y_freq * z_freq)).sum()

    return nom / denom


def run_sim(admix_params, pop_params, migr_params, *, seq_len=10_000, mut_rate=1e-8,
            debug=False, num_replicates=None):
    """Generic function for simulating Neandertal introgression."""
    CHIMP, YRI, DIN, NEA, EUR = [pop_params[p]["id"] for p in pop_params]

    # Neandertal admixture parameters
    admix_duration = years_to_gen(admix_params["duration"])
    admix_rate = admix_params["rate"] / admix_duration
    admix_start = years_to_gen(admix_params["t_admix"])
    admix_end = admix_start - admix_duration

    # population split times
    t_split_eur, t_split_dinka, t_split_nea, t_split_ch = \
        [years_to_gen(pop_params[p]["t_split"])
         for p in ["eur", "dinka", "nea", "chimp"]]

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
        msp.MigrationRateChange(time=admix_end, rate=admix_rate,
                                matrix_index=(EUR, NEA)),
        msp.MigrationRateChange(time=admix_start, rate=0,
                                matrix_index=(EUR, NEA)),

        # population size during the bottleneck
        msp.PopulationParametersChange(
            time=admix_start,
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


samples = sample_ages("data/emh_ages.txt")

admix_params = {
    "rate"     : 0.03,
    "t_admix"  : 55000,
    "duration" : 1000
}

pop_params = {
    "chimp" : {"id" : 0, "Ne" : 10000, "t_sample" : 1  * [0],     "t_split" : 6_000_000},
    "yri"   : {"id" : 1, "Ne" : 10000, "t_sample" : 5 * [0]},
    "dinka" : {"id" : 2, "Ne" : 10000, "t_sample" : 5 * [0],     "t_split" : 150_000},
    "nea"   : {"id" : 3, "Ne" : 1000,  "t_sample" : 4  * [70000], "t_split" : 500_000},
    "eur"   : {"id" : 4, "Ne" : 10000, "t_sample" : samples.age,  "t_split" : 60_000,    "bottle_Ne" : 2000}
}

migr_params = { "t" : 0, "afr_to_eur" : 0, "eur_to_afr" : 0 }

ts = run_sim(admix_params, pop_params, migr_params, seq_len=100_000_000)





def simulate_mutation_times(ts, random_seed=None):
    rng = random.Random(random_seed)
    mutation_time = np.zeros(ts.num_mutations)
    for tree in ts.trees():
        for mutation in tree.mutations():
            a = tree.time(mutation.node)
            b = tree.time(tree.parent(mutation.node))
            mutation_time[mutation.id] = rng.uniform(a, b)
    return mutation_time



def get_mutation_population(ts, mutation_time):
    node_migrations = collections.defaultdict(list)
    for migration in ts.migrations():
        node_migrations[migration.node].append(migration)

    mutation_population = np.zeros(ts.num_mutations, dtype=int)
    for tree in ts.trees():
        for site in tree.sites():
            for mutation in site.mutations:                
                mutation_population[mutation.id] = tree.population(mutation.node)
                for mig in node_migrations[mutation.node]:
                    # Stepping through all migations will be inefficient for large 
                    # simulations. Should use an interval tree (e.g. 
                    # https://pypi.python.org/pypi/intervaltree) to find all 
                    # intervals intersecting with site.position.
                    if mig.left <= site.position < mig.right:
                        # Note that we assume that we see the migration records in 
                        # increasing order of time!
                        if mig.time < mutation_time[mutation.id]:
                            assert mutation_population[mutation.id] == mig.source
                            mutation_population[mutation.id] = mig.dest
    return mutation_population

mutation_time = simulate_mutation_times(ts)
mutation_population = get_mutation_population(ts, mutation_time)



















all_snps = get_snps(ts, generate_names(pop_params))
info_snps = nea_snps(all_snps)
pd.Series((info_snps[s] == info_snps.nea0).mean() for s in samples.name).mean()


tree = IntervalTree().from_tuples([(1, 20, "a"), (30, 40, "b"), (80, 90, "c")])
sites = [1, 20, 5, 7, 18, 28, 36, 39, 75, 83, 120]




def mut_assign_pops(ts, mutation_time):
    node_migrations = defaultdict(list)
    for migration in ts.migrations():
        node_migrations[migration.node].append(migration)

    mutation_population = np.zeros(ts.num_mutations, dtype=int)
    for tree in ts.trees():
        for site in tree.sites():
            for mutation in site.mutations:
                mutation_population[mutation.id] = tree.population(mutation.node)
                migration_edges = IntervalTree().from_tuples(
                    (m.left, m.right, m) for m in node_migrations[mutation.node])
                for _, _, mig in migration_edges[site.position]:
                    if mig.time < mutation_time[mutation.id]:
                        if mutation_population[mutation.id] != mig.source:
                            return tree, site, mutation, migration_edges, node_migrations
                        mutation_population[mutation.id] = mig.dest

    return pd.Series(mutation_source), pd.Series(mutation_dest)



mutation_time = mut_assign_times(ts)
tree, site, mutation, migration_edges, node_migrations = mut_assign_pops(ts, mutation_time)

mutation_population = get_mutation_population(ts, mutation_time)



# true_nea = all_snps[list(mutation_population == pop_params["nea"]["id"])]

for mig in node_migrations[mutation.node]:
    if mig.left <= site.position < mig.right:
        print(mig)


for mig in migration_edges[site.position]:
    print(mig.data)





for mig in node_migrations[mutation.node]:
    if mig.left <= site.position < mig.right:
        if mig.time < mutation_time[mutation.id]:
            print(mig, mutation_population[mutation.id] == mig.source)
            assert mutation_population[mutation.id] == mig.source

for _, _, mig in migration_edges[site.position]:
    if mig.time < mutation_time[mutation.id]:

        print(mig, mutation_population[mutation.id] == mig.source)

