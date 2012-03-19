\name{SVMPredict}
\alias{SVMPredict}
\title{ Predict false interactions using a training set }
\description{ Using a training set to train the SVM classifier, and predict false interactions in a user given PPI network.}
\usage{
SVMPredict(training_set,predict_set,output="falsePPIs-ppiPre.csv",organism="yeast",
           drop ="IEA", replaceNA=0)
}
\arguments{
  \item{training_set}{ CSV format golden standard training set }
  \item{predict_set}{ PPI network to do the prediction }
  \item{output}{ CSV format file to ave the result }
  \item{organism}{ One of "anopheles", "arabidopsis", "bovine", "canine", "chicken", "chimp", "ecolik12", "ecsakai", "fly", "human", "malaria", "mouse", "pig", "rat", "rhesus", "worm", "xenopus", "yeast", "coelicolor" and "zebrafish."}
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
\examples{
  #edges <- data.frame(node1=c("1132", "1133", "1134", "1134", "1145", "1147"),
  #                    node2=c("1134", "1134", "1145", "1147", "1147", "1149"),
  #                    label=c(1, 1, 1, 0, 0, 0))
  #graph<-graph.data.frame(edges,directed=FALSE)
  #trainsample <- "ppiPre-SVMPredict-trainsample.csv"
  #write.csv(edges,file=trainsample,row.names=FALSE) 
  #edges <- data.frame(node1=c("1132", "1133", "1134", "1134", "1146", "1147"),
  #                    node2=c("1133", "1134", "1142", "1147", "1147", "1149"),
  #                    label=c(1, 0, 1, 0, 1, 0))
  #graph<-graph.data.frame(edges,directed=FALSE)
  #predictsample <- "ppiPre-SVMPredict-predictsample.csv"
  #write.csv(edges,file=predictsample,row.names=FALSE) 
  #SVMPredict(trainsample, predictsample, organism="human", replaceNA=0)
}

\keyword{ manip }