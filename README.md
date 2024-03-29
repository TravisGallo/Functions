# Functions
This is a repository of various functions that I have created during my work that may be useful to others.

1. `2017-08-24_urban_metric_fnc.R`: An R script created by Travis Gallo and Mason Fidino to manually calculate a metric of urbanization based off the manuscripts Czúni et al 2012 "Estimation of Urbanization Using Visual Features of Satellite Images" & 
Liker et al. "Lean birds in the city: body size and condition of house spar- rows along the urbanization gradient".
Creating square polygons was modified from "Create A Square Buffer Around a Plot Centroid in R" by Leah A. Wasser, Natalie Robinson, Sarah Elmendorf, Megan A. Jones

2. `2017-08-24_extract_landcover_parallel.R`: An R script created by Travis Gallo with help from Mason Fidino that extracts the proportion of landcover types within a given buffer size. Here we used a high resolution (1-m) landcover raster for the Chicago area. It also works with the NLCD landcover data with a few tweaks. This script is modified from Maxwell Joseph's "Spatial data extraction around buffered points in R". I simply made it work in parallel and made some tweaks to manage memory issues that we were facing since this project involved a large raster layer and many points. Obviously, there are a few project specific tweaks within the function.

3. `2018-05-07_costDistance_mod.R`: An R script created by Travis Gallo that modifies the costDistance function in the gdistance packaged so that one can calculate pairwise distances in parrallel. Thinking about dispersal distances, you can also set a cutoff distance to only calculate distances between neighborhing patches within a given distance.

4. `visualizing_priors`: a small little model to visualize priors that are not that intuitive to code up in base R. Like the half-Cauchy in the example.
