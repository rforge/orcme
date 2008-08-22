`addRows` <-
function(startData,Rdata1) {

addData1 <- meanSquaredResidue(Rdata1)
addData2 <- meanSquaredResidue(startData)
geneAddition <- addData2$geneResiduals <= addData1$overallHvalue 
addData <- startData[geneAddition, ]
Cdata <- rbind(Rdata1,addData)
return(Cdata)
}

