isoClust3 <- function(directionDoseSpecificMeans, directionObservedData, trend, lambda, phi, includeObserved){

      DRdata <- as.data.frame(directionDoseSpecificMeans)
      nDose <- repDose(doseData = trend)
      repData <- directionObservedData
      incClusterData <- if (includeObserved)
        OCDMR(repData = repData, DRdata = DRdata, alpha = 1.2, lambda = lambda,
            phi = phi, nDose = nDose)
      else 
        OCDM(DRdata = DRdata, alpha = alpha, lambda = lambda, phi = phi)
    
      incClusterRowNames <- incClusterData$clusterRowNames
      incClusterNumber <-  incClusterData$clusterNrow
      incCluster <- clusteredData(parentData = DRdata, clusterRowNames = incClusterRowNames,
          clusterNumber = incClusterNumber, monoDir = "g")
      return(incCluster)
}