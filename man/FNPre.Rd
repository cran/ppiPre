\name{FNPre}
\alias{FNPre}
\title{ Predict false negative interactions based on topological similarities}
\description{ User can predict false negative interactions from an given PPI network, based on one to three topological similarities.
}
\usage{
FNPre(file,indicator=c("RA","AA","Jaccard"),threshold=0.1, 
      output="FalseNegativePreResult-ppiPre.csv") 
}
\arguments{
  \item{file}{ CSV format PPI network given by user, each line of which includes two interacting proteins. }
  \item{indicator}{ Can be any combination of "RA", "AA", and "Jaccard", indicate the similarities used.}
  \item{threshold}{ The ratio of false negative interactions to positive interactions in the network.}
  \item{output}{ Result will be saved in the output file.}
}
\references{
[1] T. Zhou, L. Lv, and Y.-C. Zhang, "Predicting missing links via local information", The European Physical Journal B - Condensed Matter and Complex Systems, vol. 71, no. 4, pp. 623-630, Oct. 2009

[2] L. A. Adamic and E. Adar, "Friends and neighbors on the Web", Social Networks, vol. 25, no. 3, pp. 211-230, 2003.

[3] P. Jaccard, "Etude comparative de la distribution florale dans une portion des Alpes et des Jura", Bull. Soc. Vaud. Sci. Nat, vol. 37, p. 541, 1901.
}
\author{
	Yue Deng <anfdeng@163.com>
}
\seealso{
  \code{\link{AASim}}
  \code{\link{JaccardSim}}
  \code{\link{RASim}}
}
\examples{
  edges <- data.frame(node1=c("1132", "1133", "1134"),node2=c("1134", "1134", "1145"))
  graph<-graph.data.frame(edges,directed=FALSE)
  samplefile <- "ppiPre-FNPre-sample.csv"
  write.csv(edges,file=samplefile,row.names=FALSE) 
  FNPre(file=samplefile, indicator = c("RA", "AA"), threshold = 0.1)
  result<-read.csv(file="FalseNegativePreResult-ppiPre.csv")
  print(result)
}
\keyword{ manip }
