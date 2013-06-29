\name{IntelliGOGeneSim}
\alias{IntelliGOGeneSim}
\title{IntelliGO Semantic Similarity Between two Genes}
\description{Given two genes, this function will calculate the IntelliGO Semantic Similarity between them}
\usage{
IntelliGOGeneSim(gene1, gene2, w1 = 1, w2 = 1, ont = "MF", 
                 organism = "yeast", drop = "IEA")
}
\arguments{
  \item{gene1}{ Entrez gene id. }
  \item{gene2}{ Another entrez gene id.}
  \item{w1}{EC weight }
  \item{w2}{EC weight }
  \item{ont}{One of "MF", "BP", and "CC" subontologies.}
  \item{organism}{ One of "anopheles", "arabidopsis", "bovine", "canine", "chicken", "chimp", "ecolik12", "ecsakai", "fly", "human", "malaria", "mouse", "pig", "rat", "rhesus", "worm", "xenopus", "yeast", "coelicolor" and "zebrafish."}
  \item{drop}{A set of evidence codes based on which certain annotations are dropped. Use NULL to keep all GO annotations.}
}
\value{
  \item{geneSim }{Semantic similarity.}
  \item{GO1 }{Corresponding GO terms for gene1.}
  \item{GO2 }{Corresponding GO terms for gene2.}
}
\references{
[1] S. Benabderrahmane, M. Smail-Tabbone, O. Poch, A. Napoli, and M.-D. Devignes, "IntelliGO: a new vector-based semantic similarity measure including annotation origin", BMC Bioinformatics, vol. 11, no. 1, p. 588, 2010.
}
\author{
	Yue Deng <anfdeng@163.com>
}
\seealso{
  \code{\link{TCSSGeneSim}}
  \code{\link{WangGeneSim}}
}
\examples{
  #IntelliGOGeneSim("1134", "1145", ont="MF", organism="human")
  #IntelliGOGeneSim("YOL001W", "YPL031C", ont="MF", organism="yeast")
}
\keyword{ manip }
