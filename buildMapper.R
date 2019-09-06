buildMapper = function(completeData,catVars = FALSE){ 
  
  if(catVars){
    distanceValues = as.matrix(daisy(completeData,metric = "gower"))
  }else{
    distanceValues = as.matrix(daisy(completeData,metric = "euclidean"))
  }
  
  filterValues = numeric(dim(completeData)[1])
  points_in_vertex = list()
  
  for(i in 1:length(filterValues)){
    filterValues[i] = sum(exp(-((distanceValues[i,] ** 2)) / (2 * sd((distanceValues[i,])) ** 2)))
  }
   
  cent = floor(sqrt(length(filterValues) /2))
  kclust = kmeans(filterValues,centers = cent)
  for(i in 1:cent)
  {
    start = length(points_in_vertex)
    if(length(which(kclust$cluster == i)) > 1)
    {
    levelDist = as.dist(distanceValues[which(kclust$cluster == i),which(kclust$cluster ==i)])
    levelMaxDist = max(levelDist)
    levelClusterOutput = hclust(levelDist,method="single")
    heights = levelClusterOutput$height
    cut = cluster_cutoff_at_first_empty_bin(heights,levelMaxDist,10)
    clusterIndices = as.vector(cutree(levelClusterOutput,h=cut))
    for(j in 1:max(clusterIndices))
    {
      points_in_vertex[[start + j]] = which(kclust$cluster == i)[which(clusterIndices == j)]
    }
    }else{
      points_in_vertex[[start + 1]] == which(kclust$cluster == i)
    }
  }
  return(points_in_vertex)
}
