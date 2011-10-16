\name{TopologicSims}
\alias{TopologicSims}
\title{ Compute topological similarities from user input file }
\description{ Given an PPI network, this function will calculate three topological similarities of each protein pair in the network. }
\usage{
TopologicSims(inputfile, outputfile = "TopologicSims-ppiPre.csv", header = TRUE, sep = ",")
}
\arguments{
  \item{inputfile}{ CSV format PPI network, each line of which is an two interacting proteins. }
  \item{outputfile}{ Result will be saved in output file. }
  \item{header}{ If the input file contains header. }
  \item{sep}{ The separator of the input file. }
}
\references{
[1] T. Zhou, L. Lv, and Y.-C. Zhang, "Predicting missing links via local information", The European Physical Journal B - Condensed Matter and Complex Systems, vol. 71, no. 4, pp. 623-630, Oct. 2009
[2] L. A. Adamic and E. Adar, "Friends and neighbors on the Web", Social Networks, vol. 25, no. 3, pp. 211-230, 2003.
[3] P. Jaccard, "?tude comparative de la distribution florale dans une portion des Alpes et des Jura", Bull. Soc. Vaud. Sci. Nat, vol. 37, p. 541, 1901.
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
  edges <- data.frame(node1=c("1132", "1133", "1134", "1134", "1145", "1147", "1147", "1147"),node2=c("1134", "1134", "1145", "1147", "1147", "1148", "1149", "1150"))
  graph<-graph.data.frame(edges,directed=FALSE)
  samplefile <- "ppiPre-TopologicSims-sample.csv"
  write.csv(edges,file=samplefile,row.names=FALSE) 
  TopologicSims(samplefile,outputfile="TopologicSims-ppiPre.csv", header=TRUE, sep=",") 
  result<-read.csv(file="TopologicSims-ppiPre.csv")
}
\keyword{ manip }
