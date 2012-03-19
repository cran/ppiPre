\name{JaccardSim}
\alias{JaccardSim}
\title{ Compute Jaccard Index Between Two Nodes in PPI Network }
\description{ Given a network and two nodes, this function will calculate the Jaccard Index  between the two nodes
}
\usage{
JaccardSim(node1, node2, graph)
}
\arguments{
  \item{node1}{ Name of one node in the PPI network.}
  \item{node2}{ Name of another node in the PPI network.}
  \item{graph}{ A PPI network}
}
\value{
	The Jaccard Index  between the two nodes
}
\references{
[1] P. Jaccard, "Etude comparative de la distribution florale dans une portion des Alpes et des Jura", Bull. Soc. Vaud. Sci. Nat, vol. 37, p. 541, 1901.
}
\author{
	Yue Deng <anfdeng@163.com>
}
\seealso{
  \code{\link{AASim}}
  \code{\link{RASim}}
}
\examples{
  edges <- data.frame(node1=c("1132", "1133", "1134", "1134", "1145", "1147", "1147", "1147"),
                      node2=c("1134", "1134", "1145", "1147", "1147", "1148", "1149", "1150"))
  graph<-graph.data.frame(edges,directed=FALSE)
  JaccardSim("1134","1147",graph)
}
\keyword{ manip }
