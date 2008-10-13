getDirection <- function(x, y){

  y <- as.matrix(y)
  ordx <- order(x)
  x <- x[ordx]
  unx <- unique(x)
  
  if (nrow(y) == 1){
    
    y <- y[ordx]  
    y.m <- tapply(y, as.factor(x), mean)
    y.m.tot <- rep(mean(y), length(unx))
    y.is.u <- isoreg(unx,  y.m)$yf
    y.is.d <- rev(isoreg(rev(unx), y.m)$yf)
    n.p <- table(x)
    n.g <- length(n.p)
    
    iso.u <- rep.iso.u <- rep.iso.d <- y.m.all <- NULL
    rep.iso.d <- rep(y.is.d, n.p)
    rep.iso.u <- rep(y.is.u, n.p)
    y.m.all <- rep(y.m, n.p)
    
    SST0 <- sum((y - mean(y))^2)
    
    SSIS.u1 <- sum((rep.iso.u-y)^2)
    SSIS.d1 <- sum((rep.iso.d-y)^2)
    
    SST <- sum((y-y.m.all)^2)
    
    direction <- if (SSIS.u1 <= SSIS.d1) "u" else "d"
   
  } else {
    
    y <- y[, ordx]
    
    ydf <- as.data.frame(t(y))
    y.m <- do.call("cbind", unclass(by(ydf, x, mean)))
    
    y.m.tot <- matrix(rep(rowMeans(y), length(x)), ncol = length(x))
    
    y.is.u <- t(apply(y.m, 1, function(x) isoreg(unx, x)$yf))
    y.is.d <- t(apply(y.m, 1, function(x) rev(isoreg(rev(unx), x)$yf)))
  
    n.p <- table(x)
    n.g <- length(n.p)
     
    rep.iso.d <- y.is.d[, rep(1:length(n.p),n.p)]
    rep.iso.u <- y.is.u[, rep(1:length(n.p),n.p)]
    
    y.m.all <- y.m[, rep(1:length(n.p), n.p)]
  
    SST0 <- rowSums((y - rowMeans(y))^2)
  
    SSIS.u1 <- rowSums((rep.iso.u - y)^2)
    SSIS.d1 <- rowSums((rep.iso.d - y)^2)
  
    SST <- rowSums((y - y.m.all)^2)
    direction <- ifelse(SSIS.u1 <= SSIS.d1, "u", "d")
  }
  return(direction)
}
