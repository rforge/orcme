`clusteredData` <- function(parentData,clusterRowNames,clusterNumber,monoDir){
    newData <-  parentData[clusterRowNames[[1]],]
    ntmp<- list()
    if(length(unlist(clusterNumber))==1) {
       ntmp[[1]] <- rep(1,clusterNumber[[1]])
    }else{
      ntmp[[1]] <- rep(1,clusterNumber[[1]])
      for(k in 2:length(unlist(clusterNumber)) ){
        newData <-  rbind(newData,parentData[clusterRowNames[[k]],])
        ntmp[[k]]<- rep(k,clusterNumber[[k]])           
      }
    }
    gnames <- paste(monoDir,c(1:length(unique(unlist(ntmp)))),sep="")
    cluster1 <- factor(unlist(ntmp),labels=gnames)
    newData2 <- data.frame(grp=cluster1,newData)
    return(newData2)
}

