\name{ORCME}
\alias{ORCME}
\title{Order restricted clustering for dose-response trends in microarray experiments}


\description{
  The function performs delta-clustering of a microarray data. It can be used for
  clustering of both the time-course or dose-response microarray data.
}

\usage{
ORCME(DRdata,lambda,phi)
}

\arguments{
  \item{DRdata}{the microarray data with rows corresponding to genes and columns 
    corresponding to time points or different doses}
  \item{lambda}{assumed proportion of coherence relative to the observed data, 
    it ranges between 0 and 1. A lambda value of 1 considers the observed data as a 
    cluster and lambda value of 0 finds every possible pattern within the data.} 
  \item{phi}{minimum number of genes in a cluster} 
}

\value{
{the matrix of classification into clusters: each row represents one gene and columns found clusters. The matrix consist of the Booleans values, in each row there is only one of them TRUE which means that the gene was classified into the respective gene.}
}


\references{
  Cheng, Y. and Church, G. M. (2000). Biclustering of expression data. In:
     Proceedings of the Eighth International Conference on Intelligent Systems for
     Molecular Biology, 1, 93-103.
   
  Kasim, A. \emph{et al.} (2008). Order-restricted clustering for dose-response microarray
    experiments, submitted.
}

\author{Adetayo Kasim, Tobias Verbeke}

\seealso{\code{\link{monotoneDirection}}, \code{\link{plotIsomeans}}}

\examples{
  data(doseData)
  data(geneData)

  dirData <- monotoneDirection(geneData = geneData,doseData = doseData)
  incData <- as.data.frame(dirData$incData)
 
  ORCME(DRdata=incData,lambda=0.15,phi=2)
}

\keyword{cluster}
