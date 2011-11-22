\name{resampleORCME}
\alias{resampleORCME}
\title{Estimation of the proportion of the heterogeneity in the observed data for clustering}


\description{
 The function is computing within cluster sum of squares for given proportion of heterogeneity. Minimal number of genes per  cluster is fixed as 2. The sum of squares is computed through resampling the 100 data sets with 100 genes randomly sampled with  replacement from the reduced expression data.
}

\usage{
resampleORCME(clusteringData,lambdaVector)
}

\arguments{
  \item{clusteringData}{the microarray data with rows corresponding to genes and columns 
    corresponding to time points or different doses}
  \item{lambdaVector}{vector of assumed proportions of of heterogeneity of the observed data, 
    it ranges between 0 and 1. A lambda value of 1 considers the observed data as a 
    cluster and lambda value of 0 finds every possible pattern within the data} 
}

\value{
 A list of matrices that represent one of the 100 iterations. Every matrix consist of the columns 
  \item{lambda}{vector of the proportions of heterogeneity given as input}
  \item{WSS}{within clusters sum of squares for given proportion of heterogeneity} 
  \item{TSS}{total clusters sum of squares for given proportions of heterogeneity} 
  \item{nc}{number of clusters as a function for given proportions of heterogeneity} 
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

\seealso{\code{\link{ORCME}}, \code{\link{plotLambda}}}

\examples{
  data(doseData)
  data(geneData)

  dirData <- monotoneDirection(geneData = geneData,doseData = doseData)
  incData <- as.data.frame(dirData$incData)
 
  lambdaVector <- c(0.05,0.50,0.95)
  
  \donttest{
  resampleORCME(clusteringData=incData, lambdaVector=lambdaVector)
  }
}

\keyword{cluster}
