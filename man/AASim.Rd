\name{AASim}
\alias{AASim}
\title{ Compute Adamic-Adar Index Between Two Nodes in PPI Network }
\description{ Given a network and two nodes, this function will calculate the Adamic-Adar Index  between the two nodes
}
\usage{
AASim(node1, node2, graph)
}
\arguments{
  \item{node1}{ Name of one node in the PPI network.}
  \item{node2}{ Name of another node in the PPI network.}
  \item{graph}{ A PPI network}
}
\value{
	The Adamic-Adar Index  between the two nodes
}
\references{
[1] L. A. Adamic and E. Adar, "Friends and neighbors on the Web", Social Networks, vol. 25, no. 3, pp. 211-230, 2003.
}
\author{
	Yue Deng <anfdeng@163.com>
}
\seealso{
  \code{\link{RASim}}
  \code{\link{JaccardSim}}
}
\examples{
  edges <- data.frame(node1=c("1132", "1133", "1134", "1134", "1145", "1147", "1147", "1147"),
                      node2=c("1134", "1134", "1145", "1147", "1147", "1148", "1149", "1150"))
  graph<-graph.data.frame(edges,directed=FALSE)
  AASim("1134","1147",graph)
}
\keyword{ manip }

