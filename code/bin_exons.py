import random
import pandas as pd
import numpy as np

random.seed(123456789)

# load the SLiM 0-based coordinates of regions to simulate
region_coords = pd.read_table("data/slim_coords/exon_regions.bed", sep="\t")

h_scale = [round(x, 1) for x in np.linspace(0, 1, 11)]
s_scale = [round(x + 0.5, 1) for x in np.linspace(0, 1, 11)]

h_regions = region_coords.copy()
h_regions["bin_h"] = [random.choice(h_scale) for _ in range(len(region_coords))]
h_regions.to_csv("data/slim_coords/bin_h_exon_regions.bed", sep="\t")

s_regions = region_coords.copy()
s_regions["bin_s"] = [random.choice(s_scale) for _ in range(len(region_coords))]
s_regions.to_csv("data/slim_coords/bin_s_exon_regions.bed", sep="\t")
