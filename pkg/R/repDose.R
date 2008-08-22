`repDose` <-
function(doseData){                
    uniqueDose <- unique(doseData )
  nDoseReplicate <- list()
 for(j in 1:length(uniqueDose )){
    nDoseReplicate[[j]] <- length(doseData[doseData ==uniqueDose[j]])
  }
nDoseReplicate <- unlist(nDoseReplicate)
return(nDoseReplicate)
 }

