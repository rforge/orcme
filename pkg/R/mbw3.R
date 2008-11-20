
MBW3<- function(multData){
  grp <- multData$grp
  lgrp <- tapply(grp,grp,length)
  mData1 <- as.matrix(multData[,-1])
  if(length(lgrp)>1){
     mData <- MBW4(mData1[grp==unique(grp)[1],])
     for(k in 2:length(unique(grp))){
        if(dim(t(as.data.frame(mData1[grp==unique(grp)[k],])))[1]> 1){
         mData <- rbind(mData, MBW4(mData1[grp==unique(grp)[k],]))
        }else{
          mData <- rbind(mData, rep(0,ncol(mData1)))
       }
     }
    }else{
        mData <- MBW4(mData1[grp==unique(grp)[1],])
    }
  mData <- as.matrix(mData)
  mData1 <- as.matrix(mData1)
  Rmes <-( mean(mData^2) / mean((MBW4(mData1))^2))

  return(Rmes)
}