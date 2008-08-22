`dataReduction` <-
function(arrayData, startData){

rpsel <- logicalMask(Cdata=arrayData,data=startData)

clusterNonMember <- startData[rpsel == F,]
clusterMember <- startData[rpsel==T,]
clusterdata <- list(clusterMember ,clusterNonMember)
names(clusterdata ) <- c("clusterMember" , "clusterNonMember")

return(clusterdata )
}

