# An R script created by Travis Gallo with help from Mason Fidino
# This script is modified from Maxwell Josephs "Spatial data extraction around buffered points in R"
# I simply made it work in parallel and made some tweaks to manage memory issues 
# to help processing when a project involves large raster layers or many points

# Last updated 2017-08-24 by Travis Gallo

# load packages
package_load<-function(packages = NULL, quiet=TRUE, verbose=FALSE, warn.conflicts=FALSE){
  
  # download required packages if they're not already
  pkgsToDownload<- packages[!(packages  %in% installed.packages()[,"Package"])]
  if(length(pkgsToDownload)>0)
    install.packages(pkgsToDownload, repos="http://cran.us.r-project.org", quiet=quiet, verbose=verbose)
  
  # then load them
  for(i in 1:length(packages))
    require(packages[i], character.only=T, quietly=quiet, warn.conflicts=warn.conflicts)
}

# the packages needed for this analysis
packs <- c("dplyr", "reshape2", "runjags", "mcmcplots",
           "runjags", 'parallel','raster','rgdal','stringr',
           'foreach','doParallel')

package_load(packs)


# calculate the proportion of each landcover variable within a buffer

extract_bufferData <- function (buff) {
  
  # load fine scale landcover raster
  rm <- raster("~/Documents/GIS/LandCover2010ChicagoRegion/merge_NC_landcover.tif")
  
  # load points
  points <- readOGR(dsn=path.expand("~/Documents/GIS/Bat Project/Grid Points"),layer="Grid_100m_reduced")
  points$id <- as.character(points$site) #change "site" to whatever your sites are named in the SpatialPointsDataFrame
  
  #setup parallel backend to use many processors
  cl <- makeCluster(7) #set number of cores
  registerDoParallel(cl) # register backend
  
  # loop through each point
  dataLong <- foreach(i=1:length(points),.combine=rbind, .packages=c("raster", "rgdal", "stringr")) %dopar% {
    # extract land cover data for each point, given buffer size
    Landcover <- prop.table(table(extract(rm, points[i,], buffer=buff)))
    if(length(Landcover)==0) {
      Landcover <- NA
      names(Landcover) <- "BLANK"
    }
    # summarize each site's data by proportion of each cover type
    # convert to data frame
    data.frame(id = points[i,]$id,
               cover = names(Landcover),
               percent = as.numeric(Landcover)
    )
  }
  # stop cluster
  stopCluster(cl)
  
  # reshape data
  mydf_reshape <- reshape(dataLong,idvar="id",timevar="cover", direction="wide")
  mydf_reshape$id <- as.character(mydf_reshape$id)
  
  # remove dataLong to save memory
  rm(dataLong)
  
  # NA's to 0 
  mydf_reshape[is.na(mydf_reshape)] <- 0
  if(sum(grep("BLANK",colnames(mydf_reshape))) > 0){
    mydf_reshape <- mydf_reshape[,-grep("BLANK",colnames(mydf_reshape))]
  }
  
  # create cover name column
  # order by site and order columns 1-7
  df <- mydf_reshape[order((as.numeric(mydf_reshape$id))),order(names(mydf_reshape))]
  
  #setup parallel backend
  cl <- makeCluster(7) #set number of cores
  
  # 5,6,& 7 combined are impervious cover
  df$imp <- 0
  df$imp <- parApply(cl,df[,6:8],1,sum)
  
  #stop cluster
  stopCluster(cl)
  
  # remove sites column and proportions of ID #3,6 & 7 because we dont want bare ground, roads and paved areas
  df <- df[,c(1:3,9)]
  colnames(df) <- c("id",paste("tree_", buff, sep=""),paste("openveg_", buff, sep=""), paste("imp_", buff, sep=""))
  
  return(df)
  
}