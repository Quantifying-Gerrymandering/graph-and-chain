import geopandas as gpd
from gerrychain import Graph, Partition
import json
from collections import defaultdict
import pickle

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

if __name__ == "__main__":
    shapefile = "../data/shapefile_with_islands/shapefile_with_islands.shp"
    dem_votes = []
    rep_votes = []

    for file_id in range(95):
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
                dem_vote[district_id] += int(gdf.iloc[part]["PRES_Dem"])
                rep_vote[district_id] += int(gdf.iloc[part]["PRES_Rep"])
        print(f"partiition: {file_id}")
        dem_votes.append(dem_vote)
        rep_votes.append(rep_vote)

    print(dem_votes)
    print(rep_votes)
    with open('dem_votes.pkl', 'wb') as fp:
        pickle.dump(dem_votes, fp)
    with open('rep_votes.pkl', 'wb') as fp:
        pickle.dump(rep_votes, fp)