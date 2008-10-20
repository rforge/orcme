\name{timeGeneData}
\alias{timeGeneData}
\docType{data}
\title{Time Trend Gene Expression Data Example}
\description{
  Simulated time trend gene expression data as a numerical matrix; 
  the expression values reflect 5 different linear profiles with random erros from N(0,1).
 There are 10 observations    for each profile.
}
\usage{data(timeGeneData)}
\format{
  The format is:
 num [1:50, 1:15]  8.492629   8.152580   7.333815  20.370423  18.895721  ...
}
\examples{
  data(timeGeneData)
  head(timeGeneData)
}
\keyword{datasets}
