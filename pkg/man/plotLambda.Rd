\name{plotLambda}
\alias{plotLambda}
\title{Plot the variaty of the properties dependent on the proportion of heterogeneity in observed data set.}

\description{
  This function provides the plots of the dependency of the variety of properties on the proportion of heterogeneity in observed data set. It is not using the clustering as simple input, but it is also computing additional properties. The function can plot within cluster sum of squares, number of cluster, penalized within cluster sum of squares, Calsanzik and Harabasx index and Hartigan index
}

\usage{
  plotLambda(lambdaChoiceOutput,output)
}

\arguments{
  \item{lambdaChoiceOutput}{the output of the function \link{resampleORCME}}
  \item{output}{the variable that determines which output would be plotted, the values are "wss" for the cluster sum of squares, "ncluster" for the number of cluster, "pwss" for the penalized within cluster sum of squares, "ch" for the Calsanzik and Harabasx index and "h" for the Hartigan index}
}


\value{
  A plot of one of the properties mentioned above dependent on the proportion of heterogeneity. The confidence intervals are plotted instead of the point estimate.
}

\references{
  Cheng, Y. and Church, G. M. (2000). Biclustering of expression data. In:
     Proceedings of the Eighth International Conference on Intelligent Systems for
     Molecular Biology, 1, 93-103.
   
  Kasim, A. \emph{et al.} (2008). Order-restricted clustering for dose-response microarray
    experiments, submitted. % TODO complete reference
}

\author{Adetayo Kasim, Tobias Verbeke}

\seealso{\code{\link{ORCME}}, \code{\link{resampleORCME}}}

\examples{
  data(doseData)
  data(geneData)

  dirData <- monotoneDirection(geneData = geneData,doseData = doseData)
  incData <- as.data.frame(dirData$incData)
 
  lambdaVector <- seq(0.05,0.95,0.05)
  
  lambdaChoiceOutput <- resampleORCME(clusteringData=incData, lambdaVector=lambdaVector)

  plotLambda(lambdaChoiceOutput,output="wss")
  plotLambda(lambdaChoiceOutput,output="ncluster")
  plotLambda(lambdaChoiceOutput,output="pwss")
  plotLambda(lambdaChoiceOutput,output="ch")
  plotLambda(lambdaChoiceOutput,output="h")
}

\keyword{cluster}
