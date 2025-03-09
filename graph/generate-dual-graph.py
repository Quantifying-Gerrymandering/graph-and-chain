# 2019 cbg map to use for comparison: https://databasin.org/datasets/b6359a64b2fa4d19a8b38ff0c348f2d1/ 

from pathlib import Path
import pandas as pd
import geopandas as gpd
from gerrychain import Graph
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

## Read in new data file without water CBGs
# gdf = gpd.read_file("./data/final_shapefile_without_water/final_shapefile_without_water.shp") # synced shapefile generated using script from sync-cbg-shape-with-voting repo
gdf = gpd.read_file("../data/shapefile_with_islands/shapefile_with_islands.shp")

# graph = Graph.from_geodataframe(gdf, adjacency="queen")
graph = Graph.from_geodataframe(gdf, ignore_errors=False)

print(len(graph.nodes)) # 25584 CBGs that are non-water (of the 25607 CBGs total in California)

island_indices = [1316, 19333, 14151, 23207, 4685, 12109, 880, 13040, 9110, 1079, 11609]
island_geoids = [gdf.iloc[i]["GEOID20"] for i in island_indices]
print(island_geoids)


"""
Add the following edges:
060750604002 (port) - 060759804011 (farallon island) --- DONE
060750101011 (pier 39) - 060750179032 (alcatraz island) --- DONE
060750615072 - 060750179031 --- DONE
060839801001 (channel island 0) - 061110025003 (ventura harbor) --- DONE
061110036181 (channel island 1) - 061110036172 --- DONE
061119800001 (channel island 2) - 061110036172 --- DONE
060375991001 (channel island 3) - 060375991002 (channel island 4) --- DONE
060375990001 (channel island 5) - 060375760011 --- DONE
060375990001 (channel island 5) - 060379800311 --- DONE
060375991002 (channel island 4main) - 060375990002 (channel island 4main_2) # don't need
060375991002 (channel island 4main) - 060375990003 (channel island 4main_3) # don't need
060375991002 (channel island 4main) - 060375990004 (channel island 4main_4) # don't need
060730050001 - 060730110001 --- DONE

060014272005 - 060014060001 --- DONE
060730101091 - 060730102011 --- DONE
060590995143 - 060590995141 --- DONE
060590995145 - 060590995141 --- DONE
060375775011 - 060375776041 --- DONE

060590630051 - 060590630061 --- DONE
060590635001 - 060590629001 --- DONE
060375991001 - 060375990001 --- DONE
"""

# Add edge  between farallon islands and port
farallon_island_geoid = "060759804011"
farallon_island_index = gdf.index[gdf["GEOID20"] == farallon_island_geoid].tolist()[0]
port_geoid = "060750604002"
port_index = gdf.index[gdf["GEOID20"] == port_geoid].tolist()[0]
graph.add_edges_from([(farallon_island_index, port_index)])

# Add edge between alcatraz and  pier 39
alcatraz_island_geoid = "060750179032"
alcatraz_island_index = gdf.index[gdf["GEOID20"] == alcatraz_island_geoid].tolist()[0]
pier_39_geoid = "060750101011"
pier_39_index = gdf.index[gdf["GEOID20"] == pier_39_geoid].tolist()[0]
graph.add_edges_from([(alcatraz_island_index, pier_39_index)])

# Add edge between 060750615072 and 060750179031
connection1_geoid = "060750615072"
connection2_geoid = "060750179031"
connection1_index = gdf.index[gdf["GEOID20"] == connection1_geoid].tolist()[0]
connection2_index = gdf.index[gdf["GEOID20"] == connection2_geoid].tolist()[0]
graph.add_edges_from([(connection1_index, connection2_index)])

## Add edge between channel islands and their respective ports
channel_islands_geoids = ["060839801001", "061110036181", "061119800001", "060375991001", "060375991002", "060375990001", "060375990002", "060375990003", "060375990004"] # retrieved using mapshaper.org, last digit is the block group number
channel_islands_indices = gdf.index[gdf["GEOID20"].isin(channel_islands_geoids)].tolist()
# Add edge between channel island 0 and ventura harbor (geoid: 061110025003)
ventura_harbor_geoid = "061110025003"
ventura_harbor_index = gdf.index[gdf["GEOID20"] == ventura_harbor_geoid].tolist()[0]
graph.add_edges_from([(channel_islands_indices[0], ventura_harbor_index)])
# Add edges between channel island 1 and 2 and their shared port (geoid: 061110036172)
channel_island_1and2_port_geoid = "061110036172"
channel_island_1and2_port_index = gdf.index[gdf["GEOID20"] == channel_island_1and2_port_geoid].tolist()[0]
graph.add_edges_from([(channel_islands_indices[1], channel_island_1and2_port_index)])
graph.add_edges_from([(channel_islands_indices[2], channel_island_1and2_port_index)])
# Add edge between channel island 3 and channel island 4port (san clemente to avalon)
graph.add_edges_from([(channel_islands_indices[3], channel_islands_indices[5])])
graph.add_edges_from([(channel_islands_indices[3], channel_islands_indices[4])])
# Add edges between channel island 5 and its two ports (geoids: 060375760011 and 060379800311)
channel_island_5_port1_geoid = "060375760011"
channel_island_5_port1_index = gdf.index[gdf["GEOID20"] == channel_island_5_port1_geoid].tolist()[0]
channel_island_5_port2_geoid = "060379800311"
channel_island_5_port2_index = gdf.index[gdf["GEOID20"] == channel_island_5_port2_geoid].tolist()[0]
graph.add_edges_from([(channel_islands_indices[5], channel_island_5_port1_index)])
graph.add_edges_from([(channel_islands_indices[5], channel_island_5_port2_index)])
graph.add_edges_from([(channel_islands_indices[5], ventura_harbor_index)])

# Add edge between 060730050001 and 060730110001
connection3_geoid = "060730050001"
connection4_geoid = "060730110001"
connection3_index = gdf.index[gdf["GEOID20"] == connection3_geoid].tolist()[0]
connection4_index = gdf.index[gdf["GEOID20"] == connection4_geoid].tolist()[0]
graph.add_edges_from([(connection3_index, connection4_index)])

# Add edge between 060014272005 and 060014060001
connection5_geoid = "060014272005"
connection6_geoid = "060014060001"
connection5_index = gdf.index[gdf["GEOID20"] == connection5_geoid].tolist()[0]
connection6_index = gdf.index[gdf["GEOID20"] == connection6_geoid].tolist()[0]
graph.add_edges_from([(connection5_index, connection6_index)])

# Add edge between 060730101091 and 060730102011
connection7_geoid = "060730101091"
connection8_geoid = "060730102011"
connection7_index = gdf.index[gdf["GEOID20"] == connection7_geoid].tolist()[0]
connection8_index = gdf.index[gdf["GEOID20"] == connection8_geoid].tolist()[0]
graph.add_edges_from([(connection7_index, connection8_index)])

# Add edge between 060590995143 and 060590995141
connection9_geoid = "060590995143"
connection10_geoid = "060590995141"
connection11_geoid = "060590995145"
connection9_index = gdf.index[gdf["GEOID20"] == connection9_geoid].tolist()[0]
connection10_index = gdf.index[gdf["GEOID20"] == connection10_geoid].tolist()[0]
connection11_index = gdf.index[gdf["GEOID20"] == connection11_geoid].tolist()[0]
graph.add_edges_from([(connection9_index, connection10_index)])
graph.add_edges_from([(connection11_index, connection10_index)])

# Add edge between 060375775011 and 060375776041
connection12_geoid = "060375775011"
connection13_geoid = "060375776041"
connection12_index = gdf.index[gdf["GEOID20"] == connection12_geoid].tolist()[0]
connection13_index = gdf.index[gdf["GEOID20"] == connection13_geoid].tolist()[0]
graph.add_edges_from([(connection12_index, connection13_index)])

# Add edge between 060590630051 and 060590630061
connection14_geoid = "060590630051"
connection15_geoid = "060590630061"
connection14_index = gdf.index[gdf["GEOID20"] == connection14_geoid].tolist()[0]
connection15_index = gdf.index[gdf["GEOID20"] == connection15_geoid].tolist()[0]
graph.add_edges_from([(connection14_index, connection15_index)])

# Add edge between 060590635001 and 060590629001
connection16_geoid = "060590635001"
connection17_geoid = "060590629001"
connection16_index = gdf.index[gdf["GEOID20"] == connection16_geoid].tolist()[0]
connection17_index = gdf.index[gdf["GEOID20"] == connection17_geoid].tolist()[0]
graph.add_edges_from([(connection16_index, connection17_index)])

# After all edges are added, find unconnected components
components = list(nx.connected_components(graph))
print(f"\nFound {len(components)} connected components")

# After adding an edge, add this debug code:
def verify_edge(graph, geoid1, geoid2, gdf):
    """Verify that an edge exists between two nodes."""
    index1 = gdf.index[gdf["GEOID20"] == geoid1].tolist()[0]
    index2 = gdf.index[gdf["GEOID20"] == geoid2].tolist()[0]
    
    print(f"\nVerifying edge between {geoid1} and {geoid2}")
    print(f"Node indices: {index1}, {index2}")
    print(f"Edge exists: {graph.has_edge(index1, index2)}")
    print(f"Nodes exist: {index1 in graph.nodes} and {index2 in graph.nodes}")
    return index1, index2

# Example usage after adding an edge:
connection1_geoid = "060839801001"
connection2_geoid = "061110025003"
idx1, idx2 = verify_edge(graph, connection1_geoid, connection2_geoid, gdf)
graph.add_edges_from([(idx1, idx2)])
verify_edge(graph, connection1_geoid, connection2_geoid, gdf)

connection3_geoid = "060375991001"
connection4_geoid = "060375990001"
idx3, idx4 = verify_edge(graph, connection3_geoid, connection4_geoid, gdf)
graph.add_edges_from([(idx3, idx4)])
verify_edge(graph, connection3_geoid, connection4_geoid, gdf)

# Add edge between disconnected parts of District 4
connection5_geoid = "060990008061"
connection6_geoid = "060952506042"
idx5, idx6 = verify_edge(graph, connection5_geoid, connection6_geoid, gdf)
graph.add_edges_from([(idx5, idx6)])
verify_edge(graph, connection5_geoid, connection6_geoid, gdf)

# Add edge between disconnected parts of District 7
connection1_geoid = "060990009082"
connection2_geoid = "060710103002"
idx1, idx2 = verify_edge(graph, connection1_geoid, connection2_geoid, gdf)
graph.add_edges_from([(idx1, idx2)])
verify_edge(graph, connection1_geoid, connection2_geoid, gdf)

# Add edge between disconnected parts of District 9
connection1_geoid = "060770021001"
connection2_geoid = "060750132002"
idx1, idx2 = verify_edge(graph, connection1_geoid, connection2_geoid, gdf)
graph.add_edges_from([(idx1, idx2)])
verify_edge(graph, connection1_geoid, connection2_geoid, gdf)

# Add edge between disconnected parts of District 22
connection1_geoid = "060759804011"
connection2_geoid = "060750171011"
idx1, idx2 = verify_edge(graph, connection1_geoid, connection2_geoid, gdf)
graph.add_edges_from([(idx1, idx2)])
verify_edge(graph, connection1_geoid, connection2_geoid, gdf)

# Add edge between disconnected parts of District 25
connection1_geoid = "060990005063"
connection2_geoid = "060050003011"
idx1, idx2 = verify_edge(graph, connection1_geoid, connection2_geoid, gdf)
graph.add_edges_from([(idx1, idx2)])
verify_edge(graph, connection1_geoid, connection2_geoid, gdf)

# Add edge between disconnected parts of District 27
connection1_geoid = "060770051332"
connection2_geoid = "060372371012"
idx1, idx2 = verify_edge(graph, connection1_geoid, connection2_geoid, gdf)
graph.add_edges_from([(idx1, idx2)])
verify_edge(graph, connection1_geoid, connection2_geoid, gdf)

# Add edge between disconnected parts of District 39
connection1_geoid = "060990019002"
connection2_geoid = "060650408092"
idx1, idx2 = verify_edge(graph, connection1_geoid, connection2_geoid, gdf)
graph.add_edges_from([(idx1, idx2)])
verify_edge(graph, connection1_geoid, connection2_geoid, gdf)

# Add edge between disconnected parts of District 43
connection1_geoid = "060990004042"
connection2_geoid = "060190041002"
idx1, idx2 = verify_edge(graph, connection1_geoid, connection2_geoid, gdf)
graph.add_edges_from([(idx1, idx2)])
verify_edge(graph, connection1_geoid, connection2_geoid, gdf)
connection1_geoid = "060770022012"
connection2_geoid = "060530118022"
idx1, idx2 = verify_edge(graph, connection1_geoid, connection2_geoid, gdf)
graph.add_edges_from([(idx1, idx2)])
verify_edge(graph, connection1_geoid, connection2_geoid, gdf)

# Add edge between disconnected parts of District 44
connection1_geoid = "060990010022"
connection2_geoid = "060371434013"
idx1, idx2 = verify_edge(graph, connection1_geoid, connection2_geoid, gdf)
graph.add_edges_from([(idx1, idx2)])
verify_edge(graph, connection1_geoid, connection2_geoid, gdf)

# After all edges are added
components = list(nx.connected_components(graph))
print(f"\nNumber of components after adding edges: {len(components)}")

positions = {node: (row.geometry.centroid.x, row.geometry.centroid.y) 
             for node, row in gdf.iterrows()}
nx.draw(graph, pos=positions, node_size=10, edge_color="blue")
plt.show()

if len(components) > 1:
    # sort components in decreasing size
    components = sorted(components, key=len, reverse=True)
    
    print(f"\nMain component size: {len(components[0])} nodes")
    print("Disconnected components:")
    
    for i, component in enumerate(components[1:], 1):
        print(f"\nComponent {i} ({len(component)} nodes):")

        geoids = [gdf.iloc[node]["GEOID20"] for node in component]
        print("GEOIDs:", geoids)
        
        # find closest nodes in main component
        main_component = components[0]
        min_distance = float('inf')
        closest_pair = None
        
        for node1 in component:
            point1 = gdf.iloc[node1].geometry.centroid
            for node2 in main_component:
                point2 = gdf.iloc[node2].geometry.centroid
                dist = point1.distance(point2)
                if dist < min_distance:
                    min_distance = dist
                    closest_pair = (node1, node2)
        
        if closest_pair:
            print(f"Closest connection:")
            print(f"  Node 1: {gdf.iloc[closest_pair[0]]['GEOID20']}")
            print(f"  Node 2: {gdf.iloc[closest_pair[1]]['GEOID20']}")
            print(f"  Distance: {min_distance}")

print(nx.is_connected(graph))

graph.to_json("./dual-graph.json")

loaded_graph = Graph.from_json("./dual-graph.json")
print("\nVerifying saved graph:")
print(f"Number of nodes: {len(loaded_graph.nodes)}")
print(f"Number of edges: {len(loaded_graph.edges)}")
print(f"Loaded graph is connected: {nx.is_connected(loaded_graph)}")

# Visualize components
# plt.figure(figsize=(15, 15))
# positions = {node: (row.geometry.centroid.x, row.geometry.centroid.y) 
#             for node, row in gdf.iterrows()}

# colors = plt.cm.rainbow(np.linspace(0, 1, len(components)))
# for component, color in zip(components, colors):
#     nx.draw_networkx_nodes(graph, pos=positions, 
#                           nodelist=list(component),
#                           node_color=[color], 
#                           node_size=50,
#                           alpha=0.6)
# nx.draw_networkx_edges(graph, pos=positions, edge_color="gray", alpha=0.3)
# plt.title(f"Connected Components ({len(components)} total)")
# plt.savefig("connected_components.png", dpi=300, bbox_inches='tight')
# plt.show()