\name{plotLambda}
\alias{plotLambda}
\title{Plot the variaty of the properties dependent on the proportion of heterogeneity in observed data set}

\description{
  This function provides the plots of the dependency of the variety of properties on the proportion of heterogeneity in observed data set. It is not using the clustering as simple input, but it is also computing additional properties. The function can plot within cluster sum of squares, number of cluster, penalized within cluster sum of squares, Calsanzik and Harabasx index and Hartigan index.
}

\usage{
  plotLambda(lambdaChoiceOutput,output)
}

\arguments{
  \item{lambdaChoiceOutput}{the output of the function \link{resampleORCME}}
  \item{output}{the variable that determines which output would be plotted, the values are "wss" for the cluster sum of squares, "ncluster" for the number of cluster, "pwss" for the penalized within cluster sum of squares, "ch" for the Calsanzik and Harabasx index and "h" for the Hartigan index}
}


\value{
  A plot of one of the properties mentioned above dependent on the proportion of heterogeneity. The confidence intervals are plotted instead of the point estimates.
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

\author{Adetayo Kasim, Martin Otava and  Tobias Verbeke}

\seealso{\code{\link{ORCME}}, \code{\link{resampleORCME}}}

\examples{
  data(doseData)
  data(geneData)

  dirData <- monotoneDirection(geneData = geneData,doseData = doseData)
  incData <- as.data.frame(dirData$incData)
 
  lambdaVector <- c(0.05,0.50,0.95)
  
  \donttest{
  lambdaChoiceOutput <- resampleORCME(clusteringData=incData, lambdaVector=lambdaVector)
  plotLambda(lambdaChoiceOutput,output="wss")
  plotLambda(lambdaChoiceOutput,output="ncluster")
  plotLambda(lambdaChoiceOutput,output="pwss")
  plotLambda(lambdaChoiceOutput,output="ch")
  plotLambda(lambdaChoiceOutput,output="h")
  }
}

\keyword{cluster}
