\name{ComputeAllEvidences}
\alias{ComputeAllEvidences}
\title{ Compute the Biological and Topological Similarities Between Protein Pairs }
\description{ Given an input file, this function will calculate 4 biological and 3 topological similarities between each protein pair in the file.}
\usage{
ComputeAllEvidences(input,output="AllEvidences-ppiPre.csv",organism="yeast", drop ="IEA", 
                    header=TRUE, sep=",")
}
\arguments{
   \item{input}{ CSV format PPI network, each line of which includes two interacting proteins and another column which is 1 or 0 indicating that the interaction is positive or negative. }
   \item{output}{ Result will be saved in the output file.}
   \item{organism}{ One of "anopheles", "arabidopsis", "bovine", "canine", "chicken", "chimp", "ecolik12", "ecsakai", "fly", "human", "malaria", "mouse", "pig", "rat", "rhesus", "worm", "xenopus", "yeast" and "zebrafish".}
   \item{drop}{ A set of evidence codes based on which certain annotations are dropped. Use NULL to keep all GO annotations.}
   \item{header}{ If the input file contains header. }
   \item{sep}{ The separator of the input file. }
}
\references{ 
[1] S. Benabderrahmane, M. Smail-Tabbone, O. Poch, A. Napoli, and M.-D. Devignes, "IntelliGO: a new vector-based semantic similarity measure including annotation origin", BMC Bioinformatics, vol. 11, no. 1, p. 588, 2010.

[2] J. Z. Wang, Z. Du, R. Payattakool, P. S. Yu, and C.-F. Chen, "A new method to measure the semantic similarity of GO terms", Bioinformatics, vol. 23, no. 10, pp. 1274-1281, May. 2007.

[3] S. Jain and G. Bader, "An improved method for scoring protein-protein interactions using semantic similarity within the gene ontology", BMC Bioinformatics, vol. 11, no. 1, p. 562, 2010.

[4] T. Zhou, L. Lv, and Y.-C. Zhang, "Predicting missing links via local information", The European Physical Journal B - Condensed Matter and Complex Systems, vol. 71, no. 4, pp. 623-630, Oct. 2009

[5] L. A. Adamic and E. Adar, "Friends and neighbors on the Web", Social Networks, vol. 25, no. 3, pp. 211-230, 2003.

[6] P. Jaccard, "Etude comparative de la distribution florale dans une portion des Alpes et des Jura", Bull. Soc. Vaud. Sci. Nat, vol. 37, p. 541, 1901.

[7] M. Kanehisa and S. Goto, "KEGG: Kyoto Encyclopedia of Genes and Genomes", Nucleic Acids Research, vol. 28, no. 1, pp. 27 -30, 2000.

[8] G. Yu, F. Li, Y. Qin, X. Bo, Y. Wu, and S. Wang, "GOSemSim: an R package for measuring semantic similarity among GO terms and gene products", Bioinformatics, vol. 26, no. 7, pp. 976-978, Apr. 2010.
}
\author{
	Yue Deng <anfdeng@163.com>
}
\seealso{
  \code{\link{TCSSGeneSim}}
  \code{\link{IntelliGOGeneSim}}
  \code{\link{GOKEGGSims}}
}
\examples{
  #edges <- data.frame(node1=c("1132", "1133", "1134", "1134", "1145", "1147"),
  #                    node2=c("1134", "1134", "1145", "1147", "1147", "1149"),
  #                    label=c(1, 1, 1, 0, 0, 0))
  #graph<-igraph::graph.data.frame(edges,directed=FALSE)
  #samplefile <- "ppiPre-ComputeAllEvidences-sample.csv"
  #write.csv(edges,file=samplefile,row.names=FALSE) 
  #ComputeAllEvidences(input=samplefile,output="AllEvidences-ppiPre.csv",organism="human") 
  #result<-read.csv(file="AllEvidences-ppiPre.csv")
  #print(result)
}
\keyword{ manip }
