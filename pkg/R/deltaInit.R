`deltaInit` <- function(arrayData,lambda){
 
  initialHvalue <- meanSquaredResidue(arrayData)
  delta <- lambda * initialHvalue$overallHvalue  

  return(delta)
}

