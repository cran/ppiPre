\name{ppiPre-package}
\alias{ppiPre-package}
\alias{ppiPre}
\docType{package}
\title{
	  Predicting protein-protein interactions
}
\description{
	Predicting protein-protein interactions based on functional and topologic similarities
}
\details{
	This package integrates functions for computing similarities between proteins based on their GO annotation, KEGG annotation and PPI network topology. It integrates seven features(similarities) to predict PPIs using an SVM classifier. 
	Some internal functions to manipulate GO.db are derived from R package GOSemSim authored by Guangchuang Yu <guangchuangyu@gmail.com>. 
\tabular{ll}{
Package: \tab ppiPre\cr
Type: \tab Package\cr
Version: \tab 1.8\cr
Date: \tab 2015-06-24\cr
Depends:\tab R (>= 3.1.0)\cr
Imports:\tab AnnotationDbi, igraph, e1071, GOSemSim, GO.db\cr
Enhances:\tab org.Sc.sgd.db, org.Ag.eg.db, org.Xl.eg.db, org.Ss.eg.db,\cr
         \tab org.At.tair.db, org.Bt.eg.db, org.Ce.eg.db, org.Cf.eg.db,\cr
         \tab org.Dm.eg.db, org.Dr.eg.db, org.EcK12.eg.db, org.EcSakai.eg.db,\cr
         \tab org.Gg.eg.db, org.Mm.eg.db, org.Mmu.eg.db, org.Pf.plasmo.db,\cr
         \tab org.Pt.eg.db, org.Rn.eg.db, KEGG.db, org.Hs.eg.db\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}
\author{
Yue Deng, Rongjie Shao, Gang Wang, Yuanjun Sun
Maintainer: Yue Deng <anfdeng@163.com>
}
\seealso{
	\code{\link{GOKEGGSims}}
	\code{\link{IntelliGOGeneSim}}
	\code{\link{TCSSGeneSim}}
	\code{\link{SVMPredict}}
	\code{\link{TopologicSims}}
}

\keyword{ package }