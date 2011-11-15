plotCluster <- function(DRdata,doseData,ORCMEoutput,clusterID,zeroMean=FALSE,xlabel,ylabel){
      plotData <- DRdata[ORCMEoutput[,clusterID],]
      if(zeroMean){
          plotData2 <- plotData-rowMeans(plotData)
      } else{plotData2 <- plotData}
      ngene <- dim(plotData2)[1]
      plot(unique(doseData),plotData2[1,],ylim=range(plotData2),xlab=xlabel,ylab=ylabel,cex.axis=1.5,cex=1.5,cex.lab=1.5,type="n",xaxt="n")
      axis(1, at = c(1:length(unique(doseData))),labels = c(1:length(unique(doseData))),cex.axis=1.5,cex=1.5,cex.lab=1.5)
      apply(plotData2,1,function(x)lines(unique(doseData),x,cex.axis=1.5,cex=1.5,cex.lab=1.5))
}
