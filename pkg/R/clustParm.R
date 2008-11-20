
clusteringParameters <- function(directionDoseSpecificMeans, directionObservedData,
      includeObserved, nsamp){
      alpha <- 1.2 # TODO make into argument
      directionDoseSpecificMeans <- data.frame(directionDoseSpecificMeans)
      nDose <- repDose(doseData = Trend)
      parkeep <- list()
      for(q in 1:nsamp){

        if(dim(directionDoseSpecificMeans)[1] > 100){
          sampleval <- sample(c(1:dim(directionDoseSpecificMeans)[1]), 100, replace = TRUE)
          DRdata <-  directionDoseSpecificMeans[sampleval,]
          repData <- directionObservedData[sampleval,]
          row.names(DRdata) <- c(1:nrow(DRdata))
          row.names(repData) <- c(1:nrow(repData))
        }else{
          sampleval <- sample(c(1:dim(directionDoseSpecificMeans)[1]), dim(directionDoseSpecificMeans)[1], replace = TRUE)
          DRdata <-  directionDoseSpecificMeans[sampleval,]
          repData <- directionObservedData[sampleval,]
          row.names(DRdata) <- c(1:nrow(DRdata))
          row.names(repData) <- c(1:nrow(repData))
        }

        modelparm <- seq(0.05,0.9,0.05)
        modelphi <- seq(5, 20, 5)
        modelcrit <- matrix(NA, nrow=length(modelparm),ncol=length(modelphi))
        parlist <- matrix(NA, nrow=length(modelparm),ncol=length(modelphi))
        for( i in 1:length(modelparm)){
          for(j in 1:length(modelphi)){
            if (includeObserved){
              incClusterData <- OCDMR(repData=repData ,DRdata = DRdata,alpha = alpha,lambda = modelparm[i],phi = modelphi[j],nDose=nDose)
              incClusterRowNames <- incClusterData$clusterRowNames
              incClusterNumber <-  incClusterData$clusterNrow
              incCluster <- clusteredData(parentData=DRdata, clusterRowNames = incClusterRowNames,clusterNumber = incClusterNumber,monoDir = 'g')
              incCluster2 <- incCluster[order(as.numeric(rownames(incCluster))),]
              xroot <- log10(nrow(DRdata))
              lgrp <- length(unique(incCluster$grp))
              modelcrit[i,j] <- MBW3(multData = incCluster) + log10(lgrp^(1/xroot))

              parlist[i,j] <-  paste("lambda =",modelparm[i],", phi=",modelphi[j],sep="")
           
           
            } else {
            incClusterData <- OCDM(DRdata = DRdata,alpha = alpha,lambda = modelparm[i],phi = modelphi[j] )
            incClusterRowNames <- incClusterData$clusterRowNames
            incClusterNumber <-   incClusterData$clusterNrow
            incCluster <- clusteredData(parentData=DRdata,clusterRowNames = incClusterRowNames,clusterNumber = incClusterNumber,monoDir = 'g')
            xroot <- log10(nrow(DRdata))
            lgrp <- length(unique(incCluster$grp))
            modelcrit[i,j] <- MBW3(multData = incCluster) + log10(lgrp^(1/xroot))

            parlist[i,j] <-  paste("lambda =",modelparm[i],", phi=",modelphi[j],sep="")
            
          }
 
          }
      }
      
      vecmodelcrit <- as.vector(modelcrit)
      vecparlist <- as.vector(parlist)
      parkeep[[q]] <- vecparlist[vecmodelcrit==min(vecmodelcrit)]
      cat(q)
  }
  kpparval <- unlist(parkeep)
  optprop <- table(kpparval)/nsamp
  namespar <- names(optprop)
  valuepar <- as.numeric(optprop)
  parameterchoice <- namespar[valuepar == max(valuepar)]
  tmpout <- list(parameterchoice, optprop)
  names(tmpout)<- c('parameterchoice', 'optprop')
  return(tmpout)
}
