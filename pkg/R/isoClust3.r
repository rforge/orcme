isoClust3<- function(directionDoseSpecificMeans, directionObservedData,Trend, lambda, phi, includeObserved){
      alpha <- 1.2
      directionDoseSpecificMeans <- data.frame(directionDoseSpecificMeans)
      nDose <- repDose(doseData = Trend)
      DRdata <-  directionDoseSpecificMeans
      repData <- directionObservedData
        if (includeObserved){
            incClusterData <- OCDMR(repData=repData ,DRdata = DRdata,alpha = alpha,lambda = lambda,phi = phi,nDose=nDose)
            incClusterRowNames <- incClusterData$clusterRowNames
            incClusterNumber <-  incClusterData$clusterNrow
            incCluster <- clusteredData(parentData=DRdata, clusterRowNames = incClusterRowNames,clusterNumber = incClusterNumber,monoDir = 'g')
            return(incCluster)
        } else {
          incClusterData <- OCDM(DRdata = DRdata,alpha = alpha,lambda = lambda,phi = phi )
          incClusterRowNames <- incClusterData$clusterRowNames
          incClusterNumber <-   incClusterData$clusterNrow
          incCluster <- clusteredData(parentData=DRdata,clusterRowNames = incClusterRowNames,clusterNumber = incClusterNumber,monoDir = 'g')
          return(incCluster)
        }
}