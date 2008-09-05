\name{esData}
\alias{esData}
\docType{data}
\title{Example Data as ExpressionSet}
\description{
  Example Data under the form of an ExpressionSet.
  The actual data are the same as \code{geneData} and
  \code{doseData}.
}
\usage{data(esData)}
\format{
  The format is:
Formal class 'ExpressionSet' [package "Biobase"] with 6 slots
  ..@ assayData        :length 1 <environment> 
  ..@ phenoData        :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	1 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr "dose of the compound administered"
  .. .. ..@ data             :'data.frame':	20 obs. of  1 variable:
  .. .. .. ..$ dose: num [1:20] 1 1 1 1 2 2 2 2 3 3 ...
  .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ featureData      :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	0 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr(0) 
  .. .. ..@ data             :'data.frame':	100 obs. of  0 variables
  .. .. ..@ dimLabels        : chr [1:2] "featureNames" "featureColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ experimentData   :Formal class 'MIAME' [package "Biobase"] with 13 slots
                        (empty for the moment)
  ..@ annotation       : chr(0) (empty) 
  ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. ..@ .Data:List of 4
  .. .. .. ..$ : int [1:3] 2 7 0
  .. .. .. ..$ : int [1:3] 2 0 0
  .. .. .. ..$ : int [1:3] 1 1 0
  .. .. .. ..$ : int [1:3] 1 0 0
}
\seealso{\code{\link{geneData}}, \code{\link{doseData}}}
\examples{
  data(esData)
  esData
}
\keyword{datasets}
