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

setMethod("isoClust",
    signature = c("matrix", "factor"),
    function(expression, dose, alpha, lambda, phi, isoDir, includeObserved, doseResponse){
      expression <- data.frame(expression) 
      nDose <- repDose(doseData = dose)
      
      if(doseResponse) {     
        dirData <- monotoneDirection(geneData = expression,
            doseData = dose)
        incData <- data.frame(dirData$incData)
        decData <- data.frame(dirData$decData)
        repincData <-data.frame(dirData$obsincData)
        repdecData <-data.frame(dirData$obsdecData) 
        if (includeObserved){
          if ( isoDir == 'up' ){
            incClusterData <- OCDMR(repData=repincData ,DRdata = incData,alpha = alpha,lambda = lambda,phi = phi,nDose=nDose)
            incClusterRowNames <- incClusterData$clusterRowNames
            incClusterNumber <-  incClusterData$clusterNrow
            incCluster <- clusteredData(parentData=incData, clusterRowNames = incClusterRowNames,clusterNumber = incClusterNumber,monoDir = 'u')
            return(incCluster)
          } else if( isoDir == 'down' ){
            decClusterData <- OCDMR(repData=repdecData, DRdata = decData,alpha = alpha,lambda = lambda,phi = phi,nDose=nDose )
            decClusterRowNames <- decClusterData$clusterRowNames
            decClusterNumber <-   decClusterData$clusterNrow
            decCluster <- clusteredData(parentData=decData, clusterRowNames = decClusterRowNames,clusterNumber = decClusterNumber,monoDir = 'd')
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
          decCluster <- clusteredData(parentData=decData, 
              clusterRowNames=decClusterRowNames, clusterNumber = decClusterNumber,
              monoDir = 'd')
          
          isoMeansClusters <- rbind(incCluster, decCluster)
          return(isoMeansClusters)
        }
      }
    }else{
    if (includeObserved){
      repincData <- expression
      incData <- matrix(NA,nrow=nrow(expression),ncol=length(unique(dose)))
      for(i in 1:nrow(expression)){
          gene.i <- as.vector(expression[i,])
          incData[i,] <- tapply(gene.i,dose,mean)
      }
      incClusterData <- OCDMR(repData=repincData ,DRdata = incData,alpha = alpha,lambda = lambda,phi = phi,nDose=nDose)
      incClusterRowNames <- incClusterData$clusterRowNames
      incClusterNumber <-  incClusterData$clusterNrow
      incCluster <- clusteredData(parentData=incData, clusterRowNames = incClusterRowNames,clusterNumber = incClusterNumber,monoDir = 'u')
      return(incCluster)
    }else{
     incData <- expression
     incClusterData <- OCDM(DRdata = incData,alpha = alpha,lambda = lambda,phi = phi )
     incClusterRowNames <- incClusterData$clusterRowNames
     incClusterNumber <-   incClusterData$clusterNrow
     incCluster <- clusteredData(parentData=incData,clusterRowNames = incClusterRowNames,clusterNumber = incClusterNumber,monoDir = 'u')
     return(incCluster)
} 
}
}    
)

#' isoClust Method for matrix and factor
#' 
#' @param expression numeric matrix of expression values; each row represents
#'                 the expression in each of the samples (columns) 
#' @param dose factor indicating for each sample the dose that was
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
    signature = c("matrix", "numeric"),
    function(expression, dose, alpha, lambda, phi, isoDir, includeObserved,doseResponse){
      expression <- data.frame(expression) 
      nDose <- repDose(doseData = dose) 
      if(doseResponse) {     
        dirData <- monotoneDirection(geneData = expression,
            doseData = dose)
        incData <- data.frame(dirData$incData)
        decData <- data.frame(dirData$decData)
        repincData <-data.frame(dirData$obsincData)
        repdecData <-data.frame(dirData$obsdecData) 
        if (includeObserved){
          if ( isoDir == 'up' ){
            incClusterData <- OCDMR(repData=repincData ,DRdata = incData,alpha = alpha,lambda = lambda,phi = phi,nDose=nDose)
            incClusterRowNames <- incClusterData$clusterRowNames
            incClusterNumber <-  incClusterData$clusterNrow
            incCluster <- clusteredData(parentData=incData, clusterRowNames = incClusterRowNames,clusterNumber = incClusterNumber,monoDir = 'u')
            return(incCluster)
          } else if( isoDir == 'down' ){
            decClusterData <- OCDMR(repData=repdecData, DRdata = decData,alpha = alpha,lambda = lambda,phi = phi,nDose=nDose )
            decClusterRowNames <- decClusterData$clusterRowNames
            decClusterNumber <-   decClusterData$clusterNrow
            decCluster <- clusteredData(parentData=decData, clusterRowNames = decClusterRowNames,clusterNumber = decClusterNumber,monoDir = 'd')
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
          decCluster <- clusteredData(parentData=decData, 
              clusterRowNames=decClusterRowNames, clusterNumber = decClusterNumber,
              monoDir = 'd')
          
          isoMeansClusters <- rbind(incCluster, decCluster)
          return(isoMeansClusters)
        }
      }
    }else{
    if (includeObserved){
      repincData <- expression
      incData <- matrix(NA,nrow=nrow(expression),ncol=length(unique(dose)))
      for(i in 1:nrow(expression)){
          gene.i <- as.vector(expression[i,])
          incData[i,] <- tapply(gene.i,dose,mean)
      }
      incClusterData <- OCDMR(repData=repincData ,DRdata = incData,alpha = alpha,lambda = lambda,phi = phi,nDose=nDose)
      incClusterRowNames <- incClusterData$clusterRowNames
      incClusterNumber <-  incClusterData$clusterNrow
      incCluster <- clusteredData(parentData=incData, clusterRowNames = incClusterRowNames,clusterNumber = incClusterNumber,monoDir = 'u')
      return(incCluster)
    }else{
     incData <- expression
     incClusterData <- OCDM(DRdata = incData,alpha = alpha,lambda = lambda,phi = phi )
     incClusterRowNames <- incClusterData$clusterRowNames
     incClusterNumber <-   incClusterData$clusterNrow
     incCluster <- clusteredData(parentData=incData,clusterRowNames = incClusterRowNames,clusterNumber = incClusterNumber,monoDir = 'u')
     return(incCluster)
    }
  
  } 
})

#' 
#' @param expression ExpressionSet 
#' @param dose character string indicating the column name of the 
#'        phenoData that contains the dose data
#' @param alpha 
#' @param lambda 
#' @param phi 
#' @param isoDir 
#' @param includeObserved 
#' @returnType 
#' @return 
#' @author tobias
#' @export
setMethod("isoClust",
    signature = c("ExpressionSet", "character"),
    function(expression, dose, alpha, lambda, phi, isoDir, includeObserved,doseResponse){
    if (length(dose) != 1)
        stop("'dose' should be a character vector of length one indicating the name of the phenoData column that contains the dose data")
      dose <- pData(expression)[,dose]
      expression <- exprs(expression) 
      nDose <- repDose(doseData = dose) 
      if(doseResponse) {     
        dirData <- monotoneDirection(geneData = expression,
            doseData = dose)
        incData <- data.frame(dirData$incData)
        decData <- data.frame(dirData$decData)
        repincData <-data.frame(dirData$obsincData)
        repdecData <-data.frame(dirData$obsdecData) 
        if (includeObserved){
          if ( isoDir == 'up' ){
            incClusterData <- OCDMR(repData=repincData ,DRdata = incData,alpha = alpha,lambda = lambda,phi = phi,nDose=nDose)
            incClusterRowNames <- incClusterData$clusterRowNames
            incClusterNumber <-  incClusterData$clusterNrow
            incCluster <- clusteredData(parentData=incData, clusterRowNames = incClusterRowNames,clusterNumber = incClusterNumber,monoDir = 'u')
            return(incCluster)
          } else if( isoDir == 'down' ){
            decClusterData <- OCDMR(repData=repdecData, DRdata = decData,alpha = alpha,lambda = lambda,phi = phi,nDose=nDose )
            decClusterRowNames <- decClusterData$clusterRowNames
            decClusterNumber <-   decClusterData$clusterNrow
            decCluster <- clusteredData(parentData=decData, clusterRowNames = decClusterRowNames,clusterNumber = decClusterNumber,monoDir = 'd')
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
          decCluster <- clusteredData(parentData=decData, 
              clusterRowNames=decClusterRowNames, clusterNumber = decClusterNumber,
              monoDir = 'd')
          
          isoMeansClusters <- rbind(incCluster, decCluster)
          return(isoMeansClusters)
        }
      }
    }else{
    if (includeObserved){
      repincData <- expression
      incData <- matrix(NA,nrow=nrow(expression),ncol=length(unique(dose)))
      for(i in 1:nrow(expression)){
          gene.i <- as.vector(expression[i,])
          incData[i,] <- tapply(gene.i,dose,mean)
      }
      incClusterData <- OCDMR(repData=repincData ,DRdata = incData,alpha = alpha,lambda = lambda,phi = phi,nDose=nDose)
      incClusterRowNames <- incClusterData$clusterRowNames
      incClusterNumber <-  incClusterData$clusterNrow
      incCluster <- clusteredData(parentData=incData, clusterRowNames = incClusterRowNames,clusterNumber = incClusterNumber,monoDir = 'u')
      return(incCluster)
    }else{
     incData <- expression
     incClusterData <- OCDM(DRdata = incData,alpha = alpha,lambda = lambda,phi = phi )
     incClusterRowNames <- incClusterData$clusterRowNames
     incClusterNumber <-   incClusterData$clusterNrow
     incCluster <- clusteredData(parentData=incData,clusterRowNames = incClusterRowNames,clusterNumber = incClusterNumber,monoDir = 'u')
     return(incCluster)
    }
  
  } 
})
