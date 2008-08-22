`multipleNodeDeletion` <-
function(arrayData,delta,alpha,phi){
  mndData <- meanSquaredResidue(arrayData)
  overallHvalue <- mndData$overallHvalue 
  geneResiduals <- mndData$geneResiduals 
  kkeepData <- arrayData
  while (overallHvalue > delta){
    if (nrow(arrayData) > 100){
      mkeepData <- arrayData
      geneBound <- alpha * overallHvalue 
      geneResGTbound <- geneResiduals  > geneBound
      arrayData <- arrayData[!geneResGTbound,]
      mndData <- meanSquaredResidue(arrayData) #update parameters
      overallHvalue <- mndData$overallHvalue 
      geneResiduals <- mndData$geneResiduals 
    } else {
      sngReturnData <- singleNodeDeletion(arrayData = arrayData,delta = delta,phi = phi) 
      return (sngReturnData)
      break
    }
  }
  return(arrayData)
}

