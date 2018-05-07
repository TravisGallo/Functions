### Modification to gdistance::costDistance to run in parrallel and ###
### to calculate distances only between neighboring cells within some given cutoff ###
### Modified from Jacob van Etten's gdistance package ###

costDistance_mod <- function(x, fromCoords, toCoords, dist.cutoff, n.cores) {
  
  ## create a list of sites only within dispersal distance (Euclidian) from each site
  toCoords.loc <- vector("list", nrow(toCoords)) # location of each retained coordinate in the distance matrix
  toCoords.list <- vector("list", nrow(toCoords)) # list of coordinates within dispersal distance of each individual point (list in order of points)
  for(j in 1:nrow(toCoords)){
    toCoords.loc[[j]] <- which(sqrt((fromCoords[j,1]-fromCoords[,1])^2 + (fromCoords[j,2]-fromCoords[,2])^2) < dist.cutoff)
    toCoords.list[[j]] <- toCoords[toCoords.loc[[j]],]
  }
  
  # indicates which cells the fromCoords are in
  fromCells <- cellFromXY(x, fromCoords) 
  
  # have to set up toCells a bit different
  toCells <- vector("list", nrow(fromCoords)) # list to hold the toCells for each individual fromCell
  # indicate the toCells, but keep them within the list to hold true to the dispersal cutoff
  for(i in 1:length(toCoords.list)){
    toCells[[i]] <- cellFromXY(x, toCoords.list[[i]])
  }
  
  # make the transition matrix an object that igraph can read
  y <- gdistance::transitionMatrix(tr1CorrC)
  # create an adjacencyGraph to get edge weights
  if(isSymmetric(y)) {m <- "undirected"} else{m <- "directed"}
  adjacencyGraph <- graph.adjacency(y, mode=m, weighted=TRUE)
  # reclassify edge weights as cost
  E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight
  
  # parallel through each fromCell and its respective list of toCells
  cl <- makeCluster(n.cores) # setup parallel backend to use x number of cores
  registerDoParallel(cl) # register backend
  costDist <- foreach(i=1:length(fromCells), .packages=c("gdistance")) %dopar% {
    shortestPaths <- shortest.paths(adjacencyGraph, v=fromCells[i], to=toCells[[i]], mode="out", algorithm="dijkstra")
    rownames(shortestPaths) <- rownames(fromCoords)[i]
    colnames(shortestPaths) <- rownames(toCoords.list[[i]])
    return(shortestPaths)
  }
  stopCluster(cl) #stop cluster
  
  # put back together in a similar matrix as original function
  D_mat <- matrix(NA, nrow(fromCoords), nrow(fromCoords))
  for(i in 1:nrow(fromCoords)){
    D_mat[toCoords.loc[[i]],i] <- costDist[[i]]
  }
  rownames(D_mat) <- rownames(fromCoords)
  colnames(D_mat) <- rownames(fromCoords)
  
  return(D_mat)
}