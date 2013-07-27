\name{WangGeneSim}
\alias{WangGeneSim}
\title{Semantic Similarity Between Two Genes by Wang's Method}
\description{Given two genes, this function will calculate the Wang's Semantic Similarity between them}
\usage{
  WangGeneSim(gene1, gene2, ont = "MF", organism = "yeast", drop = "IEA")
}
\arguments{
  \item{gene1}{ Entrez gene id. }
  \item{gene2}{ Another entrez gene id. }
  \item{ont}{One of "MF", "BP", and "CC" subontologies. }
  \item{organism}{ One of "anopheles", "arabidopsis", "bovine", "canine", "chicken", "chimp", "ecolik12", "ecsakai", "fly", "human", "malaria", "mouse", "pig", "rat", "rhesus", "worm", "xenopus", "yeast", "coelicolor" and "zebrafish."}
  \item{drop}{A set of evidence codes based on which certain annotations are dropped. Use NULL to keep all GO annotations. }
}
\value{
  \item{geneSim}{Semantic similarity.}
  \item{GO1}{Corresponding GO terms for gene1.}
  \item{GO2}{Corresponding GO terms for gene2.}
}
\references{
[1] J. Z. Wang, Z. Du, R. Payattakool, P. S. Yu, and C.-F. Chen, "A new method to measure the semantic similarity of GO terms", Bioinformatics, vol. 23, no. 10, pp. 1274-1281, May. 2007.}
\author{
  Yue Deng <anfdeng@163.com>
}
\seealso{
  \code{\link{TCSSGeneSim}}
  \code{\link{IntelliGOGeneSim}}
}
\examples{
	#WangGeneSim("1134", "1145", ont="MF", organism="human")
	#WangGeneSim("YOR065W", "YEL024W", ont="MF", organism="yeast")
}
\keyword{ manip }
