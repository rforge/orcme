`singleNodeDeletion` <-
function(arrayData,delta,phi,kkeepData){

keepData <- arrayData
sndData <-meanSquaredResidue(arrayData)
overallHvalue <- sndData$overallHvalue 
geneResiduals <- sndData$geneResiduals 

    
while(overallHvalue > delta)  
{
     if(dim(arrayData)[1] <=  phi){  # if number of gene less or equal phi
return (arrayData)
break
}
geneResidualsMax <- max( geneResiduals )
if(sum(geneResiduals == geneResidualsMax) == nrow(arrayData)) { ## if genes have the same generesidual values
                    return(arrayData)
                            break
}
arrayData <- arrayData[geneResiduals < geneResidualsMax , ]
sndData <- meanSquaredResidue(arrayData)
overallHvalue <- sndData$overallHvalue 
geneResiduals <- sndData$geneResiduals 
}

return(arrayData)

}

