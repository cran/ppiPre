\name{SVMTrain}
\alias{SVMTrain}
\title{ Using Golden Standard Data Sets to Train an SVM Classifier }
\description{ Golden standard data sets is given by user, the trained SVM classifier will be returned}
\usage{
SVMTrain(input, replaceNA = 0)
}
\arguments{
  \item{input}{ Gold-standard data sets, including positive and negative samples }
  \item{replaceNA}{ The value to replace NA in training set }
}
\value{ The trained SVM classifier }

\author{
	Yue Deng <anfdeng@163.com>
}
\seealso{
  \code{\link{SVMPredict}}
}
\keyword{ manip }