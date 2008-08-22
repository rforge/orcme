`monotoneDirection` <-
function(geneData, arrayMean, doseData){
  geneDataMat <- as.matrix(geneData)
  H1 <- IsoGenem(doseData, geneDataMat)
  gdir <- as.vector(H1$direction)
  incData <- arrayMean [gdir == "u",]
  decData <- arrayMean [gdir == "d",]
  geneincData <- geneData [gdir == "u",]
  genedecData <- geneData [gdir == "d",]
  dirData <- list(incData,decData,geneincData ,genedecData)
  names(dirData)<-c('incData','decData','obsincData','obsdecData')
  return(dirData)

}

