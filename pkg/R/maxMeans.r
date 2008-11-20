maxMeans <-  function(timeData,timegeneData){
tmpMatrix <- matrix(NA,nrow=nrow(timegeneData),ncol=length(unique(timeData)))
for(i in 1:nrow(timegeneData)){
  tmpMatrix[i,] <- tapply(timegeneData[i,],timeData,mean)
}
 tmpMatrix <- as.data.frame(tmpMatrix)
 row.names(tmpMatrix) <- row.names(timegeneData)
 return(tmpMatrix)
}