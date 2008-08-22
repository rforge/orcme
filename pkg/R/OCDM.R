`OCDM` <-
function(DRdata,alpha,lambda,phi){
iter <- 1
clusterRowNames <- list()
clusterNrow <- list()
DDdata <- DRdata 
delta<- deltaInit(arrayData=DDdata ,lambda)
alpha <- alpha
phi <- phi
while(nrow(DDdata) > phi)
{
Cdata<-DDdata

## call multiple node deletions

clsData <- multipleNodeDeletion (arrayData = Cdata,alpha = alpha,delta = delta,phi = phi)

if(sum(dim(clsData) == dim(DDdata))== 2 ){break}

    ##remove identified clusters

cluster <- dataReduction(arrayData = clsData , startData = Cdata)

##add rows that do not increase mean squared residue of the known cluster

clusterAddRows <- addRows(startData = Cdata,Rdata1 = cluster$clusterMember)

##remove  identified clusters after row additions

cluster2 <- dataReduction(arrayData = clusterAddRows , startData = Cdata)

clusterRowNames[[iter]] <- as.numeric(row.names(cluster2$clusterMember))
clusterNrow[[iter]] <- nrow(cluster2$clusterMember)
    
    ##remaining data to be clustered
DDdata <- cluster2$clusterNonMember

    iter<-iter+1
}
clusterRowNames[[iter]] <- as.numeric(row.names(DDdata))
clusterNrow[[iter]] <-nrow(DDdata) 
tmpdata <- list(clusterRowNames,clusterNrow)

names(tmpdata ) <- c("clusterRowNames","clusterNrow")

return(tmpdata)
}

