\name{monotoneDirection}
\alias{monotoneDirection}
\title{The monotone means under increasing/decreasing trend}


\description{
  The function calculates the likelihood for the increasing and decreasing trend in the dose response for all the given genes separately gene-by-gene. The trend with the higher likelihood is chosen and the isotonic regression is applied on the means. 
}

\usage{
monotoneDirection(geneData, doseData)
}

\arguments{
  \item{geneData}{gene expression matrix for all genes}
  \item{doseData}{indicates the dose levels} 
}


\value{
 A list with components 
  \item{direction}{the direction with the higher likelihood of increasing
    (indicated by "up") or decreasing (indicated by "dn") trend.}
  \item{incData}{isotonic means with respect to dose for those genes that were classified as following the increasing trend.}
  \item{decData}{isotonic means with respect to dose for those genes that were classified as following the decreasing trend.}
  \item{obsincData}{observed gene expression matrix for those genes that were classified as following the increasing trend.}
  \item{obsdecData}{observed gene expression matrix for those genes that were classified as following the decreasing trend.}
  \item{arrayMean}{isotonic means with respect to dose for all genes.}
}


\references{
Lin D., Shkedy Z., Yekutieli D., Amaratunga D., and Bijnens, 
L. (editors). (2012) Modeling Dose-response Microarray
Data in Early Drug Development Experiments Using 
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

  dirData <- monotoneDirection(geneData = geneData,doseData = doseData)

  ## direction of monotone trend
  Direction <- dirData$direction
  ## Isotonic means for upward genes
  incData <- as.data.frame(dirData$incData)
  ##Isotonic means for downward genes
  decData <- as.data.frame(dirData$decData)
  ## observd data upward genes
  obsIncData <- as.data.frame(dirData$obsincData)
  ## observed data for downward genes
  obsDecData <- as.data.frame(dirData$obsdecData)
  ## isotonic means for all genes
  isoMeans <- as.data.frame(dirData$arrayMean)

}

\keyword{cluster}
