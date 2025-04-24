# graph-and-chain
Graph and chain generation for Quantifying Gerrymandering project.

**'graph' folder**:
This folder contains the dual graph of the shapefile and the code needed to generate it. Running "generate-dual-graph.py" will create "dual-graph.json" - a file showing the borders between California census block groups. Some borders are added and removed to take into account land connections, bridges, and ferries.

**'data' folder**:
This folder contains the shapefiles, population, and election data for 2020 census block groups (CBGs) in California.

"tl_2020_06_bg20" folder contains the TIGER/Line shapefile for 2020 California CBGs (.shp), its shape index (.shx), projection (.prj), attributes (.dbf), and the code page for the attributes (.cpg).

"cacbg20" contains a shapefile containing all non-water CBGs in California, its prison-adjusted population as of the 2020 census, the population of each ethnic group (Non-Hispanic White, Hispanic/Latino, Black, Asian, Native American, Pacific Islander), and the results of the 2020 presidential election and 2022 gubernatorial election in each CBG. 

The populations and 2020 election results are sourced from Dave's Redistricting. 2022 election results are sourced from the "2022_ca_gov" repository in this GitHub organization, which was synthesized through data from the Statewide Database. The file is in a format that can be uploaded to Dave's Redistricting as a custom election dataset.

**'chain' folder**:
This folder contains the code used to generate random districtings and evaluate their "expected" partisan distribution, which is used to evaluate whether a district is gerrymandered. If a districting plan contains an abnormally large number of Democratic or Republican districts that cannot be explained by random chance, the plan is considered gerrymandered.