from pathlib import Path
import pandas as pd
import geopandas as gpd
from gerrychain import Graph
import networkx as nx
import matplotlib.pyplot as plt

gdf = gpd.read_file("./data/ca-cbgs-graphing.zip") # synced shapefile generated using script from sync-cbg-shape-with-voting repo

# graph = Graph.from_geodataframe(gdf, adjacency="queen")
graph = Graph.from_geodataframe(gdf, ignore_errors=False)

print(len(graph.nodes)) # 25584 CBGs that are non-water (of the 25607 CBGs total in California)

def add_edge(b1, b2):
    b1_index = gdf.index[gdf["GEOID20"] == b1].to_list()[0]
    b2_index = gdf.index[gdf["GEOID20"] == b2].to_list()[0]
    graph.add_edges_from([(b1_index, b2_index)])

add_edge("060759804011", "060750604002") # Farallon Islands - mainland San Francisco
add_edge("060750615071", "060750179031") # Treasure Island - mainland San Francisco
add_edge("060750179032", "060750101011") # Alcatraz Island - mainland San Francisco
add_edge("060839801001", "061110025003") # Santa Cruz Island - Ventura (connected by ferry)
add_edge("061110036181", "061110036183") # Anacapa Island - Oxnard (connected by ferry)
add_edge("061119800001", "061110036183") # San Nicolas Island - mainland
add_edge("060375990001", "060375760011") # Santa Catalina Island - Long Beach (connected by ferry)
add_edge("060375991001", "060375991002") # San Clemente Island - Santa Catalina Island
add_edge("060590995145", "060590995144") # Trinidad Island, Huntington Beach
add_edge("060590995143", "060590995141") # Huntington Beach
add_edge("060590629001", "060590635001") # Newport Beach
add_edge("060590630062", "060590630051") # Newport Beach
add_edge("060730109001", "060730050001") # Coronado - San Diego

graph.to_json("./graph/dual-graph.json")

positions = {node: (row.geometry.centroid.x, row.geometry.centroid.y) for node, row in gdf.iterrows()}
nx.draw(graph, pos=positions, node_size=10, edge_color="blue")
plt.show()