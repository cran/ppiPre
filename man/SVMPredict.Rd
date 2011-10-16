\name{SVMPredict}
\alias{SVMPredict}
\title{ Predict false interactions using a training set }
\description{ Using a training set to train the SVM classifier, and predict false interactions in a user given PPI network.}
\usage{
SVMPredict(training_set,predict_set,output="falsePPIs-ppiPre.csv",organism="yeast",drop ="IEA", replaceNA=0)
}
\arguments{
  \item{training_set}{ CSV format golden standard training set }
  \item{predict_set}{ PPI network to do the prediction }
  \item{output}{ CSV format file to ave the result }
  \item{organism}{ One of "anopheles", "arabidopsis", "bovine", "canine", "chicken", "chimp", "ecolik12", "ecsakai", "fly", "human", "malaria", "mouse", "pig", "rat", "rhesus", "worm", "xenopus", "yeast" and "zebrafish." }
  \item{drop}{A set of evidence codes based on which certain annotations are dropped. Use NULL to keep all GO annotations.}
  \item{replaceNA}{ The value to replace NA in training and predict set. }
}
\author{
	Yue Deng <anfdeng@163.com>
}
\seealso{
  \code{\link{TopologicSims}}
  \code{\link{GOKEGGSims}}
  \code{\link{ComputeAllEvidences}}
}
\keyword{ manip }