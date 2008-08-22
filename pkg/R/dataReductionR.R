`dataReductionR` <-
function(arrayData, startData,repstartData){

rpsel <- logicalMask(arrayData,startData)

clusterNonMember <- startData[rpsel == F,]
clusterMember <- startData[rpsel==T,]
repclusterMember <- repstartData[rpsel==T,]
repclusterNonMember <- repstartData[rpsel == F,]
clusterdata <- list(clusterMember ,clusterNonMember,repclusterMember ,repclusterNonMember )
names(clusterdata ) <- c("clusterMember" , "clusterNonMember","repclusterMember ","repclusterNonMember ")

return(clusterdata )
}

