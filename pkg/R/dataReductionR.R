`dataReductionR` <- function(arrayData, startData, repstartData){
  rpsel <- logicalMask(arrayData, startData)
  
  clusterNonMember <- startData[!rpsel, ]
  clusterMember <- startData[rpsel, ]
  repclusterMember <- repstartData[rpsel, ]
  repclusterNonMember <- repstartData[!rpsel, ]
  clusterdata <- list(clusterMember, clusterNonMember, repclusterMember, repclusterNonMember)
  names(clusterdata) <- c("clusterMember", "clusterNonMember", "repclusterMember ", "repclusterNonMember")
  
  return(clusterdata)
}

