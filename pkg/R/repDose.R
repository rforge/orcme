`repDose` <- function(doseData){                
  uniqueDose <- unique(doseData)
  nDoseReplicate <- list()
  for(j in seq(along = uniqueDose)){
    nDoseReplicate[[j]] <- length(doseData[doseData == uniqueDose[j]])
  }
  nDoseReplicate <- unlist(nDoseReplicate)
  return(nDoseReplicate)
}

