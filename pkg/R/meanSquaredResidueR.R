`meanSquaredResidueR` <- function(arrayData,arrayrepData,nDose){
  InitData <- arrayData
  overallMean <- sum(InitData)  /( nrow(InitData) * ncol(InitData) )
  overallMeanMat <- matrix(overallMean ,nrow = nrow(arrayrepData),ncol = ncol(arrayrepData))
  
  geneMeans <- rowSums(InitData) / ncol(InitData)
  geneMeansMat <- matrix(geneMeans ,nrow = nrow(arrayrepData),ncol = ncol(arrayrepData))
  
  doseMeans <- colSums(InitData) / nrow(InitData)
  
  doseMeansRep <- list()
  for(k in 1:length(doseMeans )){
  doseMeansRep[[k]] <- rep(doseMeans[k],nDose[k])
  }
  
  doseMeansRep <- unlist(doseMeansRep)
  doseMeansMat <- t(matrix(t(doseMeansRep ),ncol = nrow(arrayrepData),nrow = ncol(arrayrepData)))
  
  
  geneDoseResiduals <- (arrayrepData - geneMeansMat  - doseMeansMat 
  + overallMeanMat)^2
  
  geneResiduals <- rowSums(geneDoseResiduals) / ncol(geneDoseResiduals)
  overallHvalue <- sum(geneDoseResiduals ) / (nrow(geneDoseResiduals) * ncol(geneDoseResiduals))
  
  mseData <- list(geneResiduals ,overallHvalue )
  
  names(mseData ) <- c("geneResiduals ","overallHvalue ")
  
  return(mseData)
}

