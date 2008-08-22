`deltaInitR` <- function(arrayData, arrayrepData, nDose, lambda){
  initialHvalue <- meanSquaredResidueR(arrayData=arrayData, arrayrepData=arrayrepData, nDose=nDose)
  delta <-lambda * initialHvalue$overallHvalue  
  
  return(delta)
}

