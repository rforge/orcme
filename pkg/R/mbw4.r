
MBW4<- function(arrayData){
overallMean <- sum(arrayData)  /( nrow(arrayData) * ncol(arrayData) )
overallMeanMat <- matrix(overallMean ,nrow = nrow(arrayData),ncol = ncol(arrayData))

geneMeans <- rowSums(arrayData) / ncol(arrayData)
geneMeansMat <- matrix(geneMeans ,nrow = nrow(arrayData),ncol = ncol(arrayData))

doseMeans <- colSums(arrayData) / nrow(arrayData)
doseMeansMat <- t(matrix(t(doseMeans ),ncol = nrow(arrayData),nrow = ncol(arrayData)))
geneDoseResiduals <- (arrayData - geneMeansMat  - doseMeansMat
+ overallMeanMat)


return(geneDoseResiduals )

}