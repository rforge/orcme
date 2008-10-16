`OCDMR` <- function(DRdata,lambda,alpha,repData,nDose,phi){
  iter <- 1
  clusterRowNames <- list()
  clusterNrow <- list()
  DDdata <- DRdata
  repDDdata <- repData
  
  delta<- deltaInitR(arrayData=DRdata,arrayrepData=repData,nDose=nDose,lambda=lambda)
  
  while (nrow(DDdata) > phi){
  	repCdata <- repDDdata 
  	Cdata <- DDdata
  	## call multiple node deletions
  	clsData <- multipleNodeDeletionR(arrayData=Cdata,arrayrepData=repCdata, nDose=nDose,alpha=alpha,delta=delta,phi=phi)
  
    ## remove identified clusters
  	cluster <- dataReductionR(arrayData=clsData , startData=Cdata,repstartData=repCdata )
  
    ## add rows that do not increase mean squared residue of the known cluster
  	clusterAddRows <- addRowsR(startData=Cdata,repstartData=repCdata,Rdata1=cluster$clusterMember,repRdata1=cluster$repclusterMember,nDose=nDose)
  
  	## remove  identified clusters after row additions
  	cluster2 <- dataReductionR(arrayData=clusterAddRows , startData=Cdata,repstartData=repCdata)
  
  	clusterRowNames[[iter]]<- row.names(cluster2$clusterMember)
  	clusterNrow[[iter]]<-nrow(cluster2$clusterMember)
    
    ## remaining data to be clustered
  	DDdata <- cluster2$clusterNonMember
  	repDDdata <- cluster2$repclusterNonMember
  
  	cat(iter)
  	iter <- iter+1
  }
  clusterRowNames[[iter]] <- as.numeric(row.names(DDdata))
  clusterNrow[[iter]] <- nrow(DDdata) 
  tmpdata <- list(clusterRowNames,clusterNrow)
  
  names(tmpdata ) <- c("clusterRowNames","clusterNrow")
  
  return(tmpdata)
}

