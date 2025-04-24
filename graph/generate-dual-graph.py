from pathlib import Path
import pandas as pd
import geopandas as gpd
from gerrychain import Graph
import networkx as nx
import matplotlib.pyplot as plt

root = (Path(__file__).parent).parent
shapefile_path = root / "data" / "cacbg20" / "cacbg20.shp"
gdf = gpd.read_file(shapefile_path) # shapefile with population and election data as attributes

graph = Graph.from_geodataframe(gdf, ignore_errors=False)

print(len(graph.nodes)) # 25584 CBGs that are non-water (of the 25607 CBGs total in California)

def add_edge(b1, b2):
    b1_index = gdf.index[gdf["GEOID20"] == b1].to_list()[0]
    b2_index = gdf.index[gdf["GEOID20"] == b2].to_list()[0]
    graph.add_edges_from([(b1_index, b2_index)])

def remove_edge(b1, b2):
    b1_index = gdf.index[gdf["GEOID20"] == b1].to_list()[0]
    b2_index = gdf.index[gdf["GEOID20"] == b2].to_list()[0]
    graph.remove_edges_from([(b1_index, b2_index)])

# Edges to remove
remove_edge("060750179031", "060411242003") # due to San Francisco's boundaries, CBG 060750179031 is not contiguous 
remove_edge("060750179031", "060014287002") # due to San Francisco's boundaries, CBG 060750179031 is not contiguous
remove_edge("060750601002", "060411302041") # two CBGs in the bay area that border each other on water only
remove_edge("060730099021", "060730113001") # San Diego - Coronado (not on bridge)
remove_edge("060730099021", "060730111001") # San Diego - Coronado (not on bridge)


# Edges to add
add_edge("060750179031", "060750615072") # Treasure Island - San Francisco
add_edge("060750179031", "060014017003") # Treasure Island - Oakland
add_edge("060750179032", "060750101011") # Alcatraz Island - San Francisco
add_edge("060759804011", "060759803001") # Farallon Islands - San Francisco
add_edge("060839801001", "061110025003") # Santa Cruz Island - Ventura (connected by ferry)
add_edge("061110036181", "061110036183") # Anacapa Island - Oxnard (connected by ferry)
add_edge("061119800001", "060839801001") # San Nicolas Island - Santa Barbara Island
add_edge("060839801001", "061110036183") # Santa Barbara Island - Oxnard (connected by ferry)
add_edge("060375991001", "060375991002") # San Clemente Island - Santa Catalina Island
add_edge("060375991002", "060379800311") # Santa Catalina Island - San Pedro, Los Angeles
add_edge("060375991002", "060375760011") # Santa Catalina Island - Long Beach
add_edge("060375991002", "060590628001") # Santa Catalina Island - Newport Beach

graph_filepath = root / "graph" / "dual-graph.json"
graph.to_json(graph_filepath)

positions = {node: (row.geometry.centroid.x, row.geometry.centroid.y) for node, row in gdf.iterrows()}
nx.draw(graph, pos=positions, node_size=10, edge_color="blue")
plt.show()