\name{GOKEGGSimsFromFile}
\alias{GOKEGGSimsFromFile}
\title{ GO- and KEGG- based Similarities Between two Genes }
\description{ Given an input file, this function will calculate KEGG-based similarity and three GO-based similarities between each protein pair in the file}
\usage{
GOKEGGSimsFromFile(input, output = "GOKEGGSims-ppiPre.csv", header = TRUE, sep = ",", 
                   organism = "yeast", drop = "IEA")
}
\arguments{
  \item{input}{ CSV format PPI network, each line of which is two interacting proteins. }
  \item{output}{ Result will be saved in output file. }
  \item{header}{ If the input file contains header. }
  \item{sep}{ The separator of the input file. }
  \item{organism}{ One of "anopheles", "arabidopsis", "bovine", "canine", "chicken", "chimp", "ecolik12", "ecsakai", "fly", "human", "malaria", "mouse", "pig", "rat", "rhesus", "worm", "xenopus", "yeast", "coelicolor" and "zebrafish."}
  \item{drop}{A set of evidence codes based on which certain annotations are dropped. Use NULL to keep all GO annotations.}
}
\references{
[1] S. Benabderrahmane, M. Smail-Tabbone, O. Poch, A. Napoli, and M.-D. Devignes, "IntelliGO: a new vector-based semantic similarity measure including annotation origin", BMC Bioinformatics, vol. 11, no. 1, p. 588, 2010.

[2] J. Z. Wang, Z. Du, R. Payattakool, P. S. Yu, and C.-F. Chen, "A new method to measure the semantic similarity of GO terms", Bioinformatics, vol. 23, no. 10, pp. 1274-1281, May. 2007.

[3] S. Jain and G. Bader, "An improved method for scoring protein-protein interactions using semantic similarity within the gene ontology", BMC Bioinformatics, vol. 11, no. 1, p. 562, 2010.
}
\author{
	Yue Deng <anfdeng@163.com>
}
\seealso{
  \code{\link{TCSSGeneSim}}
  \code{\link{WangGeneSim}}
  \code{\link{IntelliGOGeneSim}}
  \code{\link{GOKEGGSims}}
}
\examples{
  #edges <- data.frame(node1=c("1132", "1133", "1134", "1134", "1145", "1147"),
  #                    node2=c("1134", "1134", "1145", "1147", "1147", "1149"))
  #graph<-graph.data.frame(edges,directed=FALSE)
  #samplefile <- "ppiPre-GOKEGGSimsFromFile-sample.csv"
  #write.csv(edges,file=samplefile,row.names=FALSE) 
  #GOKEGGSimsFromFile(input=samplefile,output="GOKEGGSims-ppiPre.csv", header=TRUE, sep=",",
  #                   organism="human") 
  #result<-read.csv(file="GOKEGGSims-ppiPre.csv")
  #print(result)
}
\keyword{ manip }
