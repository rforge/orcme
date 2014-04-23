\name{plotCluster}
\alias{plotCluster}
\title{Plotting the gene specific profiles for one given cluster of genes}


\description{
  The function is plotting the profiles of the genes that belongs to the same cluster. It is not providing the clustering itself, just plotting the results of clustering from input. Optionally, the function can center the profiles around the gene-specific means. 
}

\usage{
plotCluster(DRdata, doseData, ORCMEoutput, clusterID, 
zeroMean=FALSE, xlabel, ylabel, main="")
}

\arguments{
  \item{DRdata}{the microarray data with rows corresponding to genes and columns 
    corresponding to time points or different doses}
  \item{doseData}{indicates the dose levels} 
  \item{ORCMEoutput}{the matrix of classification into clusters: each row represents one gene and columns found clusters. The matrix consist of the Booleans values, in each row there is only one of them TRUE which means that the gene was classified into the respective gene}
  \item{clusterID}{id of the cluster to be plotted}
  \item{zeroMean}{if TRUE, it centers the gene profiles around the gene-specific means, default is FALSE}
  \item{xlabel}{a title for the x axis}
  \item{ylabel}{a title for the y axis}
  \item{main}{an overall title for the plot}
}


\value{
  Plot of the gene specific profiles dependent one the dose level (or time point) that are classified into the given cluster. 
}


\references{
Lin D., Shkedy Z., Yekutieli D., Amaratunga D., and Bijnens, 
L. (editors). (2012) Modeling Dose-response Microarray
Data in Early   Drug Development Experiments Using 
R. Springer.

Cheng, Y. and Church, G. M. (2000). Biclustering of expression data. In:
     Proceedings of the Eighth International Conference on Intelligent Systems for
     Molecular Biology, 1, 93-103.
}

\author{Adetayo Kasim, Martin Otava and Tobias Verbeke}

\seealso{\code{\link{ORCME}}, \code{\link{plotIsomeans}}}

\examples{
  data(doseData)
  data(geneData)

  \donttest{
  dirData <- monotoneDirection(geneData = geneData,doseData = doseData)
  incData <- as.data.frame(dirData$incData)
  ORCMEoutput <- ORCME(DRdata=incData,lambda=0.15,phi=2)
  
  plotCluster(DRdata=incData,doseData=doseData, ORCMEoutput=ORCMEoutput,
  clusterID=4,zeroMean=FALSE, xlabel="Dose",ylabel="Gene Expression")
  }
}

\keyword{cluster}
