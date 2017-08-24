# An R script created by Travis Gallo and Mason Fidino to manually calculate a metric of urbanization based off
# the manuscripts Czúni et al 2012 "Estimation of Urbanization Using Visual Features of Satellite Images" & 
# Liker et al. "Lean birds in the city: body size and condition of house spar- rows along the urbanization gradient"
# Creating square polygons was modified from "Create A Square Buffer Around a Plot Centroid in R" by
# Leah A. Wasser, Natalie Robinson, Sarah Elmendorf, Megan A. Jones

# Last updated: 2017-08-24 by Travis Gallo

library(rgdal)
library(sp)
library(sf)
library(raster)


urbanization_score <- function (scale){
  
  # load sites (centroids)
  sites <- read.csv("2017-08-23_all_sites.csv", stringsAsFactors = FALSE)
  
  # set the radius of what will be our square buffer
  radius <- scale/2 # radius in meters
  
  # define the buffer edges based upon the buffer radius. 
  yPlus <- sites$UTM_N+radius
  xPlus <- sites$UTM_E+radius
  yMinus <- sites$UTM_N-radius
  xMinus <- sites$UTM_E-radius
  
  # calculate polygon coordinates for each buffer centroid. 
  square=cbind(xMinus,yPlus,  # NW corner
               xPlus, yPlus,  # NE corner
               xPlus,yMinus,  # SE corner
               xMinus,yMinus, # SW corner
               xMinus,yPlus)  # NW corner again - this closes the polygon
  
  # create unique ID for each polygon
  # was getting unique id error in SpatialPolygons() below when using site names
  ID=seq(0,nrow(sites),1)
  
  # create lists to be populated
  a <- vector('list', length(2))
  a_grids <- vector("list", length(2))
  
  # loop through each site and create a polygon
  for (i in 1:nrow(sites)) {
    # make it an polygon object
    a[[i]]<-Polygons(list(Polygon(matrix(square[i, ], ncol=2, byrow=TRUE))), ID[i])
    # create a grid across each polygon feature
    # grid size will be proportional to the size of the buffer and divisible by 100
    a_grids[[i]]<-as(st_make_grid(st_sfc(st_polygon(list(matrix(square[i, ], ncol=2, byrow=TRUE))), crs="+init=epsg:26916"), 
                                  n = c(scale/100, scale/100), what = 'polygons'), "Spatial")
  }
  
  # convert buffers to SpatialPolygon and assign CRS
  polys<-SpatialPolygons(a,proj4string=CRS(as.character("+init=epsg:26916")))
  
  # extract data for each cell in grid
  # load landcover raster
  rm <- raster("chicago_landcover_1m_merged.tif")
  
  # create list to hold proportion tables
  # create matrix to hold values
  prop_cell <- vector('list', length(a_grids[[1]]))
  class_mat <- matrix(0,nrow=length(a_grids), ncol=5)
  
  # function to calculate Rb
  doRb <- function(y=NULL, build=5){
    if(sum(build == as.numeric(names(y))) > 0){
      w5 <- which(as.numeric(names(y)) == build)
      answer <- ifelse(y[w5] > 0.5, 2, 1)
    } else {
      answer <- 0}
    return(answer)
  }
  
  # function to caculate Rv
  doRv <- function(y=NULL, veg=c(1,2)){
    if(sum(veg %in% as.numeric(names(y))) > 0){
      w5 <- which(as.numeric(names(y)) %in% veg)
      answer <- ifelse(sum(y[w5]) > 0.5, 2, 1)
    } else {
      answer <- 0}
    return(answer)
  }
  
  # function to calculate Pr
  doPr <- function(y=NULL, rd=6){
    answer <- ifelse(rd %in% as.numeric(names(y)), 1, 0)
    return(answer)
  }
  
  # extract info from raster
  for (k in 1:length(a_grids)){
    for(i in 1:length(a_grids[[1]])){
      prop_cell[[i]] <- prop.table(table(extract(rm, a_grids[[k]][i])))
    }
    # run through functions to caculate metrics and build matrix
    class_mat[k,1] <- mean(sapply(prop_cell,doRb))
    class_mat[k,2] <- sum(sapply(prop_cell,doRb) == 2)
    class_mat[k,3] <- mean(sapply(prop_cell,doRv))
    class_mat[k,4] <- sum(sapply(prop_cell,doRv) == 2)
    class_mat[k,5] <- sum(sapply(prop_cell,doPr) ==1)
  }
  
  # caculate first principle component for each site
  urban_score <- princomp(class_mat)
  
  # data frame with info from sites
  data <- cbind(sites,urban_score$scores[,1])
  colnames(data) <- c("site","UTM_E","UTM_N","score")
  
  return(list(classifications=class_mat,
              sq_buffers = polys,
              site_data=data))
}