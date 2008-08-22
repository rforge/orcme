`isoMean` <-
function(geneData,doseData){ 
  x.res <- doseData
  dat.mat <- as.data.frame(geneData )
  unx <- unique(x.res)
  ydf <- as.data.frame(t(dat.mat))
  y.m <- do.call("cbind", unclass(by(ydf, x.res, mean)))
  y.is.u <- t(apply(y.m, 1, function(x) isoreg(unx, x)$yf))
  y.is.d <- t(apply(y.m, 1, function(x) rev(isoreg(rev(unx),
            x)$yf)))
  ud <- cbind(y.is.u,y.is.d)
  iso.dir <- IsoGenem(x.res,as.matrix(dat.mat))[[11]]
  iso.mean <- matrix(0,nrow(geneData) ,length(unique(doseData)))
  for (i in 1: nrow(geneData) ) {
    if (iso.dir[i]=="u"){
       iso.mean[i,] <- y.is.u[i,]
    } else {
  iso.mean[i,] <- y.is.d[i,]

  }

  }
    return(iso.mean)
}

