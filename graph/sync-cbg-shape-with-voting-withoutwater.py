import os
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
from gerrychain import Graph, Partition
import networkx as nx

shapefile_path = "../data/tl_2020_06_bg20/tl_2020_06_bg20.shp"
cbg_path = "../data/2020-pres-ca-cbg.csv"

gdf = gpd.read_file(shapefile_path)[["GEOID20", "geometry", "STATEFP20", "COUNTYFP20"]]
cbgdf = pd.read_csv(cbg_path, 
                    usecols=["GEOID20", "Name", "T_20_CENS_ADJ_Total", "E_20_PRES_Total", "E_20_PRES_Dem", "E_20_PRES_Rep"], 
                    dtype={"GEOID20": str})

gdf["FIPS"] = gdf["STATEFP20"] + gdf["COUNTYFP20"]
gdf.drop(["STATEFP20", "COUNTYFP20"], axis=1, inplace=True)
cbgdf = cbgdf.rename(columns={"T_20_CENS_ADJ_Total": "CENS_Total",
                            "E_20_PRES_Total": "PRES_Total",
                            "E_20_PRES_Dem": "PRES_Dem",
                            "E_20_PRES_Rep": "PRES_Rep"})
gdf["geometry"] = gdf["geometry"].buffer(0)
gdf = pd.merge(gdf, cbgdf, on="GEOID20", how="inner")

gdf = gdf.to_crs(3310) # Reproject to projected coordinate system (3310 = California Albers)

 # Remove water CBGs
gdf = gdf[gdf["Name"] != "Block Group 0"]
print("Filtered out water CBGs successfully")

# connect islands to land nodes by ferry, all channel islands connect to ventura harbor in ventura county (geoid: 061110025003)
# https://geocoding.geo.census.gov/geocoder/geographies/onelineaddress?address=2950%20Pierpont%20Blvd%2C%20Ventura%2C%20CA%2093001&benchmark=4&vintage=4
# https://www.nps.gov/chis/planyourvisit/island-transportation.htm 

# GEOIDs for islands (first one is near the bay, remaining five are the channel islands)
island_geoids = ["060759804011", "060839801001", "061110036181", "061119800001", "060375991002", "060375991001", "060375990001", "060375990002", "060375990003", "060375990004"] # retrieved using mapshaper.org, last digit is the block group number
island_indices = gdf.index[gdf["GEOID20"].isin(island_geoids)].tolist()
print("Island indices: ", island_indices)

print("index of error cbg:", gdf.index[gdf["GEOID20"] == "060530114002"].tolist()[0]) # 060530114002

# # Remove islands from shapefile
# gdf = gdf.drop(index=island_indices)
# print("Removed islands successfully")

gdf.info()

# Save the merged dataframe to a shapefile
if not os.path.exists("../data/shapefile_with_islands"):
    os.makedirs("../data/shapefile_with_islands")
gdf.to_file("../data/shapefile_with_islands/shapefile_with_islands.shp", driver="ESRI Shapefile")

# Plot the gdf colored by FIPS code
gdf.plot("FIPS", edgecolor="black", linewidth=0.1)

# Convert gdf to dual graph 
graph = Graph.from_geodataframe(gdf, ignore_errors=False)
print("number of nodes: ", len(graph.nodes))

plt.show()
