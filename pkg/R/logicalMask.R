<<<<<<< .mine
`logicalMask` <- function(Cdata,data){
  rprow <- as.vector(row.names(data))
  rpcdata <- as.vector(row.names(Cdata))
=======
`logicalMask` <- function(Cdata, data){
  rprow <- row.names(data)
  rpcdata <- row.names(Cdata)
>>>>>>> .r16
  rpsel <- rprow %in% rpcdata  
  return(rpsel)
}

