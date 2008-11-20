`monotoneDirection` <- function(geneData, doseData){
  geneDataMat <- as.matrix(geneData)
  H1 <- isominmax(compData=geneDataMat,doseData=doseData)
  gdir <- as.vector(H1$Direction)
  arrayMean <- H1$isomeans
  incData <- arrayMean [gdir == "up",]
  decData <- arrayMean [gdir == "dn",]
  geneincData <- geneData [gdir == "up",]
  genedecData <- geneData [gdir == "dn",]
  dirData <- list(gdir,incData,decData,geneincData ,genedecData,arrayMean)
  names(dirData)<-c('direction','incData','decData','obsincData','obsdecData','arrayMean')
  
  return(dirData)
}

