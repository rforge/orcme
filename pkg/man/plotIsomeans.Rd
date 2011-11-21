\name{plotIsomeans}
\alias{plotIsomeans}
\title{Plot of the observed gene expression and the isotonic means with respect to dose}

\description{
 The function is plotting the observed data points of the gene expression and isotonic means with respect to dose for one particular gene. 
}

\usage{
plotIsomeans(monoData, obsData, doseData, geneIndex)
}

\arguments{
  \item{monoData}{isotonic means with respect to dose for all genes}
  \item{obsData}{observed gene expression for all genes}
  \item{doseData}{indicates the dose levels}
  \item{geneIndex}{index of the gene to be plotted}
}

\value{
  Plot of the data points and the isotonic means for each dose with the isotonic regression curve.

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

\seealso{\code{\link{ORCME}}, \code{\link{monotoneDirection}}}

\examples{
  data(doseData)
  data(geneData)

  dirData <- monotoneDirection(geneData = geneData,doseData = doseData)
  incData <- as.data.frame(dirData$incData)
  obsIncData <- as.data.frame(dirData$obsincData)
  
  ## gene-specific profile plot
  plotIsomeans(monoData=incData,obsData=obsIncData,doseData=
  doseData,geneIndex=10)

}

\keyword{cluster}

