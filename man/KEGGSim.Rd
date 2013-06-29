\name{KEGGSim}
\alias{KEGGSim}
\title{ KEGG Semantic Similarity Between two Genes }
\description{ Given two genes, this function will calculate the KEGG semantic similarity between them, based on the shared pathways between the genes.}
\usage{
KEGGSim(protein1, protein2)
}
\arguments{
  \item{protein1}{ Entrez gene id. }
  \item{protein2}{ Another entrez gene id.}
}

\value{ The KEGG semantic similarity, using the Jaccard index.}
\references{
[1] M. Kanehisa and S. Goto, "KEGG: Kyoto Encyclopedia of Genes and Genomes", Nucleic Acids Research, vol. 28, no. 1, pp. 27 -30, 2000.

[2] P. Jaccard, "Etude comparative de la distribution florale dans une portion des Alpes et des Jura", Bull. Soc. Vaud. Sci. Nat, vol. 37, p. 541, 1901.
}
\author{
	Yue Deng <anfdeng@163.com>
}
\seealso{
  \code{\link{GOKEGGSims}}
  \code{\link{GOKEGGSimsFromFile}}
}
\examples{
  KEGGSim("YJL026W","YGR180C")  
  KEGGSim("1134","1145")
}
\keyword{ manip }
