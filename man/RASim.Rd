\name{RASim}
\alias{RASim}
\title{ Compute Resource Allocation Index Between Two Nodes in PPI Network }
\description{ Given a network and two nodes, this function will calculate the Resource Allocation Index  between the two nodes
}
\usage{
RASim(node1, node2, graph)
}
\arguments{
  \item{node1}{ Name of one node in the PPI network.}
  \item{node2}{ Name of another node in the PPI network.}
  \item{graph}{ A PPI network}
}
\value{
	The Resource Allocation Index  between the two nodes
}
\references{
[1] T. Zhou, L. Lv, and Y.-C. Zhang, "Predicting missing links via local information", The European Physical Journal B - Condensed Matter and Complex Systems, vol. 71, no. 4, pp. 623-630, Oct. 2009.
}
\author{
	Yue Deng <anfdeng@163.com>
}
\seealso{
  \code{\link{AASim}}
  \code{\link{JaccardSim}}
}
\examples{
  edges <- data.frame(node1=c("1132", "1133", "1134", "1134", "1145", "1147", "1147", "1147"),
                      node2=c("1134", "1134", "1145", "1147", "1147", "1148", "1149", "1133"))
  graph<-igraph::graph.data.frame(edges,directed=FALSE)
  RASim("1134","1147",graph)
}
\keyword{ manip }
