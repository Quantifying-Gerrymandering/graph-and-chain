import geopandas as gpd
from gerrychain import Graph, Partition
import networkx as nx
import matplotlib.pyplot as plt
import random
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
import numpy as np
from collections import deque
import heapq
from collections import defaultdict

from gerrychain.updaters import Tally, cut_edges

from gerrychain import MarkovChain
from gerrychain.constraints import single_flip_contiguous
from gerrychain.proposals import propose_random_flip
from gerrychain.accept import always_accept

# Load Graph
# graph = Graph.from_json("./graphs/shapefile_with_islands.json")
graph = Graph.from_json("../graph/dual-graph.json")
print(nx.is_connected(graph))
# print(len(graph))

# Read shapefile
gdf = gpd.read_file("../data/shapefile_with_islands/shapefile_with_islands.shp")

PARTITIONS = 52


def grow_districts(graph, num_districts, gdf):
    """
    Randomly select num_districts starting nodes and grow districts with preference for compactness.
    Uses a priority queue to grow from nodes closest to district centers.
    """
    
    nodes = list(graph.nodes)
    seeds = random.sample(nodes, num_districts)
    
    district_assignment = {node: None for node in graph.nodes}
    
    # Use priority queues instead of deques
    boundaries = {i: [] for i in range(num_districts)}
    centers = {i: gdf.iloc[seeds[i]].geometry.centroid for i in range(num_districts)}
    
    for i, seed in enumerate(seeds):
        district_assignment[seed] = i
        gdf.loc[seed, 'district_id'] = i
        # Add neighbors to priority queue with distance as priority
        for neighbor in graph.neighbors(seed):
            if district_assignment[neighbor] is None:
                dist = centers[i].distance(gdf.iloc[neighbor].geometry.centroid)
                heapq.heappush(boundaries[i], (dist, neighbor))

    # Expand districts using priority queue
    print("Expanding")
    while any(boundaries.values()):
        for district, boundary in boundaries.items():
            if not boundary:
                continue

            # Get closest unassigned neighbor
            while boundary:
                _, node = heapq.heappop(boundary)
                if district_assignment[node] is None:
                    break
            else:
                continue
            
            district_assignment[node] = district
            gdf.loc[node, 'district_id'] = district
            
            # Add unassigned neighbors to priority queue
            for neighbor in graph.neighbors(node):
                if district_assignment[neighbor] is None:
                    dist = centers[district].distance(gdf.iloc[neighbor].geometry.centroid)
                    heapq.heappush(boundaries[district], (dist, neighbor))

    return district_assignment

# def grow_districts(graph, num_districts, gdf):
#     """
#     Randomly select num_districts starting nodes and grow until all nodes are covered.
#     """
    
#     nodes = list(graph.nodes)
#     seeds = random.sample(nodes, num_districts)
    
#     district_assignment = {node: None for node in graph.nodes}
    
#     boundaries = {i: deque() for i in range(num_districts)}
    
#     for i, seed in enumerate(seeds):
#         district_assignment[seed] = i
#         gdf.loc[seed, 'district_id'] = i
#         boundaries[i].append(seed)

#     # Expand districts using BFS
#     print("Expanding")
#     while any(boundaries.values()):
#         for district, boundary in boundaries.items():
#             if not boundary:
#                 continue

#             node = boundary.popleft()
            
#             for neighbor in graph.neighbors(node):
#                 if district_assignment[neighbor] is None:
#                     district_assignment[neighbor] = district
#                     gdf.loc[neighbor, 'district_id'] = district
#                     boundary.append(neighbor)

#     return district_assignment

# def grow_districts(graph, num_districts, gdf):
#     """
#     Randomly select num_districts starting nodes and grow districts with preference for compactness.
#     Keeps trying until population deviation is under 70%.
#     """
    
#     while True:  # Keep trying until we get a valid partition
#         nodes = list(graph.nodes)
#         seeds = random.sample(nodes, num_districts)
        
#         district_assignment = {node: None for node in graph.nodes}
        
#         # Use priority queues instead of deques
#         boundaries = {i: [] for i in range(num_districts)}
#         centers = {i: gdf.iloc[seeds[i]].geometry.centroid for i in range(num_districts)}
        
#         for i, seed in enumerate(seeds):
#             district_assignment[seed] = i
#             gdf.loc[seed, 'district_id'] = i
#             # Add neighbors to priority queue with distance as priority
#             for neighbor in graph.neighbors(seed):
#                 if district_assignment[neighbor] is None:
#                     dist = centers[i].distance(gdf.iloc[neighbor].geometry.centroid)
#                     heapq.heappush(boundaries[i], (dist, neighbor))

#         # Expand districts using priority queue
#         print("Expanding")
#         while any(boundaries.values()):
#             for district, boundary in boundaries.items():
#                 if not boundary:
#                     continue

#                 # Get closest unassigned neighbor
#                 while boundary:
#                     _, node = heapq.heappop(boundary)
#                     if district_assignment[node] is None:
#                         break
#                 else:
#                     continue
                
#                 district_assignment[node] = district
#                 gdf.loc[node, 'district_id'] = district
                
#                 # Add unassigned neighbors to priority queue
#                 for neighbor in graph.neighbors(node):
#                     if district_assignment[neighbor] is None:
#                         dist = centers[district].distance(gdf.iloc[neighbor].geometry.centroid)
#                         heapq.heappush(boundaries[district], (dist, neighbor))

#         # Check population deviation using graph node attributes
#         populations = defaultdict(int)
#         for node, district in district_assignment.items():
#             populations[district] += graph.nodes[node]["CENS_Total"]  # Using Census Total Population
        
#         total_pop = sum(populations.values())
#         ideal_pop = total_pop / num_districts
#         max_dev = max(abs(pop - ideal_pop) / ideal_pop for pop in populations.values())
        
#         if max_dev <= 0.70:  # 70% threshold
#             print(f"Found valid partition with {max_dev*100:.1f}% max deviation")
#             return district_assignment
#         else:
#             print(f"Retrying - max deviation was {max_dev*100:.1f}%")


# district_assignment = grow_districts(graph, PARTITIONS, gdf)
# if gdf['district_id'].isna().sum() > 0:
#     print(f"Warning: There are {gdf['district_id'].isna().sum()} NaN values in 'district_id' column.")
#     print([gdf.index[gdf['district_id'].isna()]])
# print("Unique districts:", np.sort(gdf['district_id'].unique()))

# # Assign the 'district_id' attribute to all nodes in the graph
# # nx.set_node_attributes(graph, district_assignment, name='district_id')
# # Ensure 'district_id' is assigned correctly in gdf
# for node in graph.nodes:
#     graph.nodes[node]["district_id"] = district_assignment[node]
# graph.to_json("./graphs/shapefile_with_islands.json")
# gdf.to_file("./shapefile_with_islands/shapefile_with_islands.shp", driver="ESRI Shapefile")

# # Plot
# base_colors = cm.tab20.colors * 3
# shuffled_colors = np.random.permutation(base_colors[:52])
# random_cmap = ListedColormap(shuffled_colors)
# gdf.plot("district_id", 
#          cmap=random_cmap,
#          edgecolor="black",
#          linewidth=0.1)  
# plt.show()