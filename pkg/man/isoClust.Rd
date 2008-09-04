\name{isoClust}
\alias{isoClust}
\title{Clustering Algorithm for Dose-Response Data under Monotone Constraints}
\description{
  Clustering algorithm for dose response data under monotone constraints. 
  The algorithm evaluates clusters based on coherence within the cluster 
  relative to observed heterogeneity within the microarray data. 
  It is a node-deletion-based algorithm, starting with all the genes and 
  deleting a gene or group of genes that improves the mean squared residue
  score for the cluster until the delta or phi criterium is satisfied.
}
\usage{
isoClust(geneData, doseData, alpha, lambda, phi, isoDir, includeObserved)
}
\arguments{
  \item{geneData}{A data frame of dose response microarray data, with rows corresponding 
    to genes and coloumns corresponding to the different doses under which the array were obtained. 
    The order of the array should be the same as the order of the doses in the 
    dose data}
  \item{doseData}{contains information about the doses under which the array was obtained}
  \item{alpha}{delta inflating factor for multiple row deletion}
  \item{lambda}{assumed proportion of coherence relative to the observed data, it ranges between 
    0 and 1. A lambda value of 1 considers the observed data as a cluster and lambda value of 0 
    finds every possible pattern within the data.}
  \item{phi}{minimum number of genes in a cluster}
  \item{isoDir}{the monotone direction of genes to be clustered. \code{"u"} corresponds to the 
  clustering of monotone increasing genes, \code{"d"} to monotone decreasing genes, and \code{"b"}
  to clustering of both monotone increasing and decreasing genes}
  \item{includeObserved}{boolean indicating whether the clustering should be based on both isotonic means
    and observed data (\code{TRUE}) or only on isotonic means (\code{FALSE})}
}
\details{
  Order-restricted clustering is a node deletion based algorithm. It draws inspiration on the delta 
  biclustering algorithm of Cheng and Church, 2000)
}
\value{
Data frame with input data and one additional column \code{grp}
indicating the group membership. For monotone increasing profiles,
the groups are named \code{u1}, \code{u2}, etc.; for monotone
decreasing profiles, the groups are named \code{d1}, \code{d2}, etc.
}
\references{
  Cheng, Y. and Church, G. M. (2000). Biclustering of expression data. In:
     Proceedings of the Eighth International Conference on Intelligent Systems for
     Molecular Biology, 1, 93-103.
   
  Kasim, A. \emph{et al.} (2008). Order-restricted clustering for dose-response microarray
    experiments, submitted. % TODO complete reference
}
\author{Adetayo Kasim, Tobias Verbeke}
\seealso{\code{\link{OCDM}}, \code{\link{OCDMR}}}
\examples{
  data(doseData)
  data(geneData)
  incIso <- isoClust(geneData = geneData, doseData = doseData, 
    alpha = 1.2, lambda = 0.05, phi = 10, isoDir = "up", 
    includeObserved = FALSE)
  head(incIso) # cluster membership added as a column
}
\keyword{cluster}
