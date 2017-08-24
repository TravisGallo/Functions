# Functions
This is a repository of various functions that I have created during my work that may be useful to others.

1. 2017-08-24_urban_metric_fnc.R: An R script created by Travis Gallo and Mason Fidino to manually calculate a metric of urbanization based off the manuscripts Czúni et al 2012 "Estimation of Urbanization Using Visual Features of Satellite Images" & 
Liker et al. "Lean birds in the city: body size and condition of house spar- rows along the urbanization gradient".
Creating square polygons was modified from "Create A Square Buffer Around a Plot Centroid in R" by Leah A. Wasser, Natalie Robinson, Sarah Elmendorf, Megan A. Jones

2. 2017-08-24_extract_landcover_parallel.R: An R script created by Travis Gallo with help from Mason Fidino that extracts the proportion of landcover types within a given buffer size. Here we used a high resolution (1-m) landcover raster for the Chicago area. It also works with the NLCD landcover data with a few tweaks. This script is modified from Maxwell Joseph's "Spatial data extraction around buffered points in R". I simply made it work in parallel and made some tweaks to manage memory issues that we were facing since this project involved a large raster layer and many points. Obviously, there are a few project specific tweaks within the function.

