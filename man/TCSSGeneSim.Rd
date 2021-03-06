\name{TCSSGeneSim}
\alias{TCSSGeneSim}
\title{Topological Clustering Semantic Similarity(TCSS) Between two Genes}
\description{Given two genes, this function will calculate the TCSS between them}
\usage{
TCSSGeneSim(gene1, gene2, ont = "MF", organism = "yeast", drop = "IEA")
}
\arguments{
  \item{gene1}{ Entrez gene id. }
  \item{gene2}{ Another entrez gene id.}
  \item{ont}{One of "MF", "BP", and "CC" subontologies.}
  \item{organism}{ One of "anopheles", "arabidopsis", "bovine", "canine", "chicken", "chimp", "ecolik12", "ecsakai", "fly", "human", "malaria", "mouse", "pig", "rat", "rhesus", "worm", "xenopus", "yeast" and "zebrafish."}
  \item{drop}{A set of evidence codes based on which certain annotations are dropped. Use NULL to keep all GO annotations.}
}
\value{
  \item{geneSim }{Semantic similarity.}
  \item{GO1 }{Corresponding GO terms for gene1.}
  \item{GO2 }{Corresponding GO terms for gene2.}
}
\references{
[1] S. Jain and G. Bader, "An improved method for scoring protein-protein interactions using semantic similarity within the gene ontology", BMC Bioinformatics, vol. 11, no. 1, p. 562, 2010.

[2] G. Yu, F. Li, Y. Qin, X. Bo, Y. Wu, and S. Wang, "GOSemSim: an R package for measuring semantic similarity among GO terms and gene products", Bioinformatics, vol. 26, no. 7, pp. 976-978, Apr. 2010.
}
\author{
	Yue Deng <anfdeng@163.com>
}
\seealso{
  \code{\link{IntelliGOGeneSim}}
}
\examples{
  #TCSSGeneSim("1134", "1145", ont="MF", organism="human")
  #TCSSGeneSim("YOR065W", "YEL024W", ont="CC", organism="yeast", drop="IEA")
}
\keyword{ manip }
