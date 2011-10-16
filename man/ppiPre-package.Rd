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
	This package implements several functions useful for computing similarities between proteins based on their GO annotation, KEGG annotation and network topology. It integrates thees similarities to predict PPIs using an SVM classifier. Four species are supported, including Fly, Human, Worm, and Yeast.
\tabular{ll}{
Package: \tab ppiPre\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-07-06\cr
Depends:\tab AnnotationDbi, GO.db, KEGG.db, org.Sc.sgd.db, igraph, e1071\cr
Suggests:\tab org.Ce.eg.db, org.Dm.eg.db, org.Hs.eg.db\cr
License: \tab GPL Version 2\cr
LazyLoad: \tab yes\cr
}
}
\author{
Yue Deng

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