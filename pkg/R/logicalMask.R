`logicalMask` <- function(Cdata,data){
  rprow <- as.numeric(row.names(data))
  rpcdata <- as.numeric(row.names(Cdata))
  rpsel <- rprow %in% rpcdata  
  return(rpsel)
}

