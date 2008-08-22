`multipleNodeDeletionR` <-
function(arrayData,arrayrepData,delta,nDose,alpha,phi){

mndData <-meanSquaredResidueR(arrayData=arrayData,arrayrepData=arrayrepData,nDose=nDose)
overallHvalue <- mndData$overallHvalue 
geneResiduals <- mndData$geneResiduals 

while(overallHvalue > delta){

mkeepData <- arrayData

geneBound <- alpha * overallHvalue 
geneResGTbound <- geneResiduals  > geneBound
arrayData <- arrayData[geneResGTbound ==F,]
arrayrepData<- arrayrepData[geneResGTbound ==F,]


#update parameters

mndData <-meanSquaredResidueR(arrayData=arrayData,arrayrepData=arrayrepData,nDose=nDose)
overallHvalue <- mndData$overallHvalue 
geneResiduals <- mndData$geneResiduals 

if((sum(dim(mkeepData )==dim(arrayData))==2)|nrow(arrayData) <= 100){

sngReturnData <-singleNodeDeletionR(arrayData=arrayData ,arrayrepData=arrayrepData,nDose=nDose,delta=delta,phi=phi) 

return (sngReturnData )
break
}
}

return(arrayData)

}

