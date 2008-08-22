`meanSquaredResidue` <- function(arrayData){
  overallMean <- sum(arrayData) / ( nrow(arrayData) * ncol(arrayData) )
  overallMeanMat <- matrix(overallMean , nrow = nrow(arrayData), ncol = ncol(arrayData))
  
  geneMeans <- rowSums(arrayData) / ncol(arrayData)
  geneMeansMat <- matrix(geneMeans, nrow = nrow(arrayData),ncol = ncol(arrayData))
  
  doseMeans <- colSums(arrayData) / nrow(arrayData)
  doseMeansMat <- t(matrix(t(doseMeans ), ncol = nrow(arrayData),nrow = ncol(arrayData)))
  geneDoseResiduals <- (arrayData - geneMeansMat  - doseMeansMat 
  + overallMeanMat)^2
  
  geneResiduals <- rowSums(geneDoseResiduals) / ncol(geneDoseResiduals)
  overallHvalue <- sum(geneDoseResiduals ) / (nrow(geneDoseResiduals) * ncol(geneDoseResiduals))
  
  mseData <- list(geneResiduals ,overallHvalue )
  names(mseData ) <- c("geneResiduals ","overallHvalue ")
  
  return(mseData )
}

