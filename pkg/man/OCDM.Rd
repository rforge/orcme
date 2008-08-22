\name{OCDM}
\alias{OCDM}
\title{Delta-Clustering of Microarray Data}
\description{
  This function performs delta-clustering of a microarray data. it can be used for
  clustering of both the time-course or dose-response microarray data.
}
\usage{
  OCDM(DRdata, alpha, lambda, phi)
}
\arguments{
  \item{DRdata}{the microarray data with rows corresponding to genes and columns 
    corresponding to time points or different doses}
  \item{alpha}{delta inflating factor for multiple row deletion}
  \item{lambda}{assumed proportion of coherence relative to the observed data, 
    it ranges between 0 and 1. A lambda value of 1 considers the observed data as a 
    cluster and lambda value of 0 finds every possible pattern within the data.}
  \item{phi}{minimum number of genes in a cluster}
}
% TODO \details{}
\value{
  A list with the following two components:
  \item{clusterRowNames}{a list with as many components as clusters; for each cluster
    the row names of genes are given.}
  \item{clusterNumber}{a list with as many components as clusters; for each cluster
    the number of genes are given.}
}
\references{
  Cheng, Y. and Church, G. M. (2000). Biclustering of expression data. In:
     Proceedings of the Eighth International Conference on Intelligent Systems for
     Molecular Biology, 1, 93-103.
   
  Kasim, A. \emph{et al.} (2008). Order-restricted clustering for dose-response microarray
    experiments, submitted. % TODO complete reference
}
\author{Adetayo Kasim, Tobias Verbeke}
\seealso{\code{\link{isoClust}}, \code{\link{OCDMR}}}
% TODO \examples{}
\keyword{cluster}
