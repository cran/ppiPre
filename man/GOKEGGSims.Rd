\name{GOKEGGSims}
\alias{GOKEGGSims}
\title{ GO- and KEGG- based Similarities Between two Genes }
\description{Given two genes, this function will calculate KEGG-based similarity and three GO-based similarities between them}
\usage{
GOKEGGSims(gene1, gene2, organism = "yeast", drop = "IEA")
}

\arguments{
  \item{gene1}{ Entrez gene id. }
  \item{gene2}{ Another entrez gene id.}
  \item{organism}{ One of "anopheles", "arabidopsis", "bovine", "canine", "chicken", "chimp", "ecolik12", "ecsakai", "fly", "human", "malaria", "mouse", "pig", "rat", "rhesus", "worm", "xenopus", "yeast" and "zebrafish."}
  \item{drop}{A set of evidence codes based on which certain annotations are dropped. Use NULL to keep all GO annotations.}
}
\references{
[1] S. Benabderrahmane, M. Smail-Tabbone, O. Poch, A. Napoli, and M.-D. Devignes, "IntelliGO: a new vector-based semantic similarity measure including annotation origin", BMC Bioinformatics, vol. 11, no. 1, p. 588, 2010.

[2] J. Z. Wang, Z. Du, R. Payattakool, P. S. Yu, and C.-F. Chen, "A new method to measure the semantic similarity of GO terms", Bioinformatics, vol. 23, no. 10, pp. 1274-1281, May. 2007.

[3] S. Jain and G. Bader, "An improved method for scoring protein-protein interactions using semantic similarity within the gene ontology", BMC Bioinformatics, vol. 11, no. 1, p. 562, 2010.

[4 ]G. Yu, F. Li, Y. Qin, X. Bo, Y. Wu, and S. Wang, "GOSemSim: an R package for measuring semantic similarity among GO terms and gene products", Bioinformatics, vol. 26, no. 7, pp. 976-978, Apr. 2010.
}
\author{
	Yue Deng <anfdeng@163.com>
}
\seealso{
  \code{\link{TCSSGeneSim}}
  \code{\link{IntelliGOGeneSim}}
  \code{\link{GOKEGGSimsFromFile}}
}
\examples{
  #GOKEGGSims("1134", "1145", organism = "human", drop = "IEA")
  #GOKEGGSims("YER133W", "YBR045C", organism = "yeast")
}

\keyword{ manip }
