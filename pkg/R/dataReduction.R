`dataReduction` <-
function(arrayData, startData){

rpsel <- logicalMask(Cdata=arrayData,data=startData)

clusterNonMember <- startData[!rpsel, ]
clusterMember <- startData[rpsel, ]
clusterdata <- list(clusterMember ,clusterNonMember)
names(clusterdata ) <- c("clusterMember" , "clusterNonMember")

return(clusterdata )
}

