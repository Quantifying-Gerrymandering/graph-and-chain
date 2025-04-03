from pathlib import Path
import pandas as pd
import geopandas as gpd
from gerrychain import Graph
import networkx as nx
import matplotlib.pyplot as plt

root = (Path(__file__).parent).parent
shapefile_path = root / "data" / "ca-cbgs-graphing.zip"
gdf = gpd.read_file(shapefile_path) # synced shapefile generated using script from sync-cbg-shape-with-voting repo


# graph = Graph.from_geodataframe(gdf, adjacency="queen")
graph = Graph.from_geodataframe(gdf, ignore_errors=False)

print(len(graph.nodes)) # 25584 CBGs that are non-water (of the 25607 CBGs total in California)

def add_edge(b1, b2):
    b1_index = gdf.index[gdf["GEOID20"] == b1].to_list()[0]
    b2_index = gdf.index[gdf["GEOID20"] == b2].to_list()[0]
    graph.add_edges_from([(b1_index, b2_index)])

add_edge("060759804011", "060750604002") # Farallon Islands - mainland San Francisco
add_edge("060750615072", "060750179031") # Treasure Island - mainland San Francisco (Bay Bridge)
add_edge("060014017003", "060750179031") # Treasure Island - Oakland (Bay Bridge)
add_edge("060750179032", "060750101011") # Alcatraz Island - mainland San Francisco
add_edge("060750601002", "060411311001") # Golden Gate Bridge (San Francisco - Marin County)
add_edge("060014061003", "060014271001") # Alameda - Oakland
add_edge("060014061002", "060014271003") # Alameda - Oakland
add_edge("060133570001", "060952506042") # Carquinez Strait (Crockett - Vallejo)
add_edge("060133200011", "060952521021") # Benicia - Martinez bridge
add_edge("060411122022", "060133780001") # Richmond - San Rafael bridge
add_edge("060816082002", "060014371011") # San Mateo - Hayward bridge
add_edge("060816118001", "060014443031") # Dunbarton bridge
add_edge("060839801001", "061110025003") # Santa Cruz Island - Ventura (connected by ferry)""
add_edge("061110036181", "061110036183") # Anacapa Island - Oxnard (connected by ferry)
add_edge("061119800001", "061110036183") # San Nicolas Island - mainland
add_edge("060375990001", "060375760011") # Santa Catalina Island - Long Beach (connected by ferry)
add_edge("060375991001", "060375991002") # San Clemente Island - Santa Catalina Island
add_edge("060590995145", "060590995144") # Trinidad Island, Huntington Beach
add_edge("060590995143", "060590995141") # Huntington Beach
add_edge("060590629001", "060590635001") # Newport Beach
add_edge("060590630062", "060590630051") # Newport Beach
add_edge("060730109001", "060730050001") # Coronado - San Diego


graph_filepath = root / "graph" / "dual-graph.json"
graph.to_json(graph_filepath)

positions = {node: (row.geometry.centroid.x, row.geometry.centroid.y) for node, row in gdf.iterrows()}
nx.draw(graph, pos=positions, node_size=10, edge_color="blue")
plt.show()