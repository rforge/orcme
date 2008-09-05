#' Order-Restricted Clustering
#' 
#' TODO
#' @param expression object containing the expression values
#'        from the experiment
#' @param dose object indicating the doses for each sample in
#'        the experiment
#' @param ... further arguments to be passed to the methods,
#'        currently none are used
#' @return returns a factor indicating cluster membership for all genes
#'        in the dataset; the order of the values reflects the order of
#'        the genes (rows) in the original expression data 
#' @author Tobias Verbeke
#' @export
setGeneric("isoClust", function(expression, dose, ...){
    standardGeneric("isoClust")
})    

#' isoClust Method for matrix and factor
#' 
#' @param geneData numeric matrix of expression values; each row represents
#'                 the expression in each of the samples (columns) 
#' @param doseData factor indicating for each sample the dose that was
#'                 administered
#' @param alpha 
#' @param lambda 
#' @param phi 
#' @param isoDir 
#' @param includeObserved 
#' @returnType factor 
#' @return factor indicating cluster membership for each gene; the order
#'         used is the initial order of the genes
#' @export
setMethod("isoClust",
    signature = c("matrix", "factor"),
    function(geneData, doseData, alpha, lambda, phi, isoDir, includeObserved){
      geneData <- data.frame(geneData)        
      arrayMean <- isoMean(geneData, doseData) 
      dirData <- monotoneDirection(geneData = geneData, arrayMean = arrayMean,
          doseData = doseData)
      incData <- data.frame(dirData$incData)
      decData <- data.frame(dirData$decData)
      repincData <-data.frame(dirData$obsincData)
      repdecData <-data.frame(dirData$obsdecData) 
      nDose <- repDose(doseData = doseData)
      if (includeObserved){
        if ( isoDir == 'up' ){
          incClusterData <- OCDMR(repData=repincData ,DRdata = incData,alpha = alpha,lambda = lambda,phi = phi,nDose=nDose)
          incClusterRowNames <- incClusterData$clusterRowNames
          incClusterNumber <-  incClusterData$clusterNrow
          incCluster <- clusteredData(parentData=incData,clusterRowNames = incClusterRowNames,clusterNumber = incClusterNumber,monoDir = 'u')
          return(incCluster)
        } else if( isoDir == 'down' ){
          decClusterData <- OCDMR(repData=repdecData ,DRdata = decData,alpha = alpha,lambda = lambda,phi = phi,nDose=nDose )
          decClusterRowNames <- decClusterData$clusterRowNames
          decClusterNumber <-   decClusterData$clusterNrow
          decCluster <- clusteredData(parentData=decData,clusterRowNames = decClusterRowNames,clusterNumber = decClusterNumber,monoDir = 'd')
          return(decCluster)
        } else {
          incClusterData <- OCDMR(repData=repincData ,DRdata = incData,alpha = alpha,lambda = lambda,phi = phi,nDose=nDose )
          incClusterRowNames <- incClusterData$clusterRowNames
          incClusterNumber <-   incClusterData$clusterNrow
          incCluster <- clusteredData(parentData=incData, clusterRowNames = incClusterRowNames,clusterNumber = incClusterNumber,monoDir = 'u')
          
          decClusterData <- OCDMR(repData=repdecData ,DRdata = decData,alpha = alpha,lambda = lambda,phi = phi,nDose=nDose )
          decClusterRowNames <- decClusterData$clusterRowNames
          decClusterNumber <-   decClusterData$clusterNrow
          decCluster <- clusteredData(parentData=incData,clusterRowNames = decClusterRowNames,clusterNumber = decClusterNumber,monoDir = 'd')
          isoMeansClusters <- rbind(incCluster,decCluster)
          return(isoMeansClusters)
        }
      } else{
        if (isoDir == 'up'){
          incClusterData <- OCDM(DRdata = incData,alpha = alpha,lambda = lambda,phi = phi )
          incClusterRowNames <- incClusterData$clusterRowNames
          incClusterNumber <-   incClusterData$clusterNrow
          incCluster <- clusteredData(parentData=incData,clusterRowNames = incClusterRowNames,clusterNumber = incClusterNumber,monoDir = 'u')
          return(incCluster)
        } else if(isoDir == 'down' ){
          decClusterData <- OCDM(DRdata = decData,alpha = alpha,lambda = lambda,phi = phi )
          decClusterRowNames <- decClusterData$clusterRowNames
          decClusterNumber <-   decClusterData$clusterNrow
          decCluster <- clusteredData(parentData=decData,clusterRowNames = decClusterRowNames,clusterNumber = decClusterNumber,monoDir = 'd')
          return(decCluster)
        } else {
          incClusterData <- OCDM(DRdata=incData, alpha=alpha, lambda=lambda, phi=phi)
          incClusterRowNames <- incClusterData$clusterRowNames
          incClusterNumber <-   incClusterData$clusterNrow
          incCluster <- clusteredData(parentData=incData, clusterRowNames = incClusterRowNames,clusterNumber = incClusterNumber,monoDir = 'u')
          decClusterData <- OCDM(DRdata = decData, alpha = alpha, lambda = lambda, phi = phi)
          decClusterRowNames <- decClusterData$clusterRowNames
          decClusterNumber <-   decClusterData$clusterNrow
          decCluster <- clusteredData(parentData=decData,clusterRowNames = decClusterRowNames,clusterNumber = decClusterNumber,monoDir = 'd')
          
          isoMeansClusters <- rbind(incCluster, decCluster)
          return(isoMeansClusters)
        }
        
      }
})

