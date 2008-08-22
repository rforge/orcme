`addRowsR` <-
function(startData,repstartData,Rdata1,repRdata1,nDose) {

addData1 <- meanSquaredResidueR(arrayData=Rdata1,arrayrepData=repRdata1,nDose=nDose)
addData2 <- meanSquaredResidueR(arrayData=startData,arrayrepData=repstartData,nDose=nDose)

geneAddition <- addData2$geneResiduals <= addData1$overallHvalue 

addData <- startData[geneAddition  == T, ]

Cdata <- rbind(Rdata1,addData)

return(Cdata)

}

