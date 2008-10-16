`logicalMask` <- function(Cdata, data){
  rprow <- row.names(data)
  rpcdata <- row.names(Cdata)
  rpsel <- rprow %in% rpcdata  
  return(rpsel)
}

