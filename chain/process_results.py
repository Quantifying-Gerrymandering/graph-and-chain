import geopandas as gpd
from gerrychain import Graph, Partition
import json
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

if __name__ == "__main__":
    shapefile = "../data/shapefile_with_islands/shapefile_with_islands.shp"
    dem_votes = []
    rep_votes = []

    for file_id in range(74):
        partition_assignment_file = f"results/chain-final-partitions/spanning_tree_final_partition{file_id}_recom.json"

        with open(partition_assignment_file) as f:
            assignment = json.load(f)
        assignment = {int(k): v for k, v in assignment.items()}

        gdf = gpd.read_file(shapefile)
        gdf.drop([f"district_{i}" for i in range(1, 10)], axis=1, inplace=True)
        gdf.drop([f"district{i}" for i in range(10, 91)], axis=1, inplace=True)

        # gdf.info()

        gdf["district_i"] = gdf.index.map(assignment)
        district_to_index_map = gdf.groupby("district_i").apply(lambda x: x.index.tolist()).to_dict()

        dem_vote = defaultdict(int)
        rep_vote = defaultdict(int)
        for district_id, parts in district_to_index_map.items():
            for part in parts:
                dem_vote[district_id] += gdf.iloc[part]["PRES_Dem"]
                rep_vote[district_id] += gdf.iloc[part]["PRES_Rep"]
        print(f"partiition: {file_id}")
        print(sum(dem_vote[i] > rep_vote[i] for i in range(52)))
        print(sum(dem_vote[i] < rep_vote[i] for i in range(52)))
        dem_votes.append(int(sum(dem_vote[i] > rep_vote[i] for i in range(52))))
        rep_votes.append(int(sum(dem_vote[i] < rep_vote[i] for i in range(52))))

    print(dem_votes)
    print(rep_votes)

    count, bins, _ = plt.hist(dem_votes, bins=12, density=True, alpha=0.6, color='b', edgecolor='black')
    mu, sigma = np.mean(dem_votes), np.std(dem_votes)
    x = np.linspace(40, 52, 100)
    pdf = norm.pdf(x, mu, sigma)
    plt.plot(x, pdf, 'r', linewidth=2, label=f'Normal Fit ($\\mu$={mu:.2f}, $\\sigma$={sigma:.2f})')

    plt.xticks(np.arange(min(dem_votes), max(dem_votes) + 1, 1))

    # Labels and title
    plt.xlabel('Value')
    plt.ylabel('Density')
    plt.title('Dem Seats')
    plt.legend()
    plt.show()

