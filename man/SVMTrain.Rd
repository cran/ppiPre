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
\examples{
  #edges <- data.frame(node1=c("1132", "1133", "1134", "1134", "1145", "1147"),
  #                    node2=c("1134", "1134", "1145", "1147", "1147", "1149"),
  #                    label=c(1, 1, 1, 0, 0, 0))
  #graph<-graph.data.frame(edges,directed=FALSE)
  #samplefile <- "ppiPre-ComputeAllEvidences-sample.csv"
  #write.csv(edges,file=samplefile,row.names=FALSE) 
  #ComputeAllEvidences(input=samplefile,output="AllEvidences-svmtrainsample.csv",organism="human", 
  #                    drop ="IEA", header=TRUE, sep=",") 
  #SVMTrain("AllEvidences-svmtrainsample.csv", replaceNA=0)
}
\keyword{ manip }