\name{ORCME}
\alias{ORCME}
\title{Order restricted clustering for dose-response trends in microarray experiments}


\description{
  The function performs delta-clustering of a microarray data. It can be used for
  clustering of both the time-course or dose-response microarray data.
}

\usage{
ORCME(DRdata, lambda, phi, robust=FALSE)
}

\arguments{
  \item{DRdata}{matrix of a microarray data with rows corresponding to genes and columns 
    corresponding to time points or different doses}
  \item{lambda}{assumed proportion of coherence relative to the observed data, 
    it ranges between 0 and 1. A lambda value of 1 considers the observed data as a 
    cluster and lambda value of 0 finds every possible pattern within the data.} 
  \item{phi}{minimum number of genes in a cluster} 
  \item{robust}{logical variable that determines, if algorithm uses robust version based on median polish and 
  absolute values, instead of mean square error. Default is \code{FALSE}.} 
}

\value{
{The matrix of classification into clusters: each row represents one gene and columns found clusters. The matrix consist of the Booleans values, in each row there is only one of them \code{TRUE} which means that the gene was classified into the respective cluster.}
}


\references{
Lin D., Shkedy Z., Yekutieli D., Amaratunga D., and Bijnens, 
L. (editors). (2012) Modeling Dose-response Microarray Data in EarlyDrug Development Experiments Using 
R. Springer.

Cheng, Y. and Church, G. M. (2000). Biclustering of expression data. In:
     Proceedings of the Eighth International Conference on Intelligent Systems for
     Molecular Biology, 1, 93-103.
}

\author{Adetayo Kasim, Martin Otava and Tobias Verbeke}

\seealso{\code{\link{monotoneDirection}}, \code{\link{plotIsomeans}}}

\examples{
  data(doseData)
  data(geneData)

  dirData <- monotoneDirection(geneData = geneData,doseData = doseData)
  incData <- as.data.frame(dirData$incData)
  
  \donttest{
  print(orcme <- ORCME(DRdata=incData,lambda=0.15,phi=2))
  orcmeRobust <- ORCME(DRdata=incData,lambda=0.15,phi=2, robust=TRUE)
  
  # number of genes within cluster
  colSums(orcme)
  colSums(orcmeRobust)
  }
}

\keyword{cluster}
