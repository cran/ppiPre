.initial<-function(pos = 1,envir = as.environment(pos)){
  if(!exists("ppiPreEnv") || length(ppiPreEnv)<1) {
    print("initializing ppiPre package ...")		
    assign("ppiPreEnv",new.env(),envir=envir)  
    assign("ppiPreCache", new.env(),envir=envir)
    print("finished.")
  }
}
################
KEGGSim <- function(protein1, protein2)    # KEGG-based similarity of two proteins
{

    if(!require("KEGG.db")){
    	stop("package KEGG.db is needed.")
    }
	path1 <- KEGGEXTID2PATHID[[protein1]]
	path2 <- KEGGEXTID2PATHID[[protein2]]
	intersec <- length(na.omit(match(path1, path2)))
	if(intersec==0)
		sim<-0
	else
		sim<-intersec/(length(path1)+length(path2)-intersec)
	return(sim)
}
GOKEGGSims <- function(gene1, gene2, organism="yeast", drop ="IEA")  #KEGG- and GO-based similarity of two proteins
{
	wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine","canine","anopheles","ecsakai","chicken","chimp","malaria","rhesus","pig","xenopus", "coelicolor"))

	dropcodes <- drop
	Sims <- data.frame(protein1=gene1,protein2=gene2,BPWang=0,MFWang=0,CCWang=0,BPTCSS=0,MFTCSS=0,CCTCSS=0,BPIG=0,MFIG=0,CCIG=0,KEGGSim=0)

	Sims[[3]][1]<-WangGeneSim(gene1,gene2,ont="BP",organism=wh_organism,drop = dropcodes)$geneSim
	Sims[[4]][1]<-WangGeneSim(gene1,gene2,ont="MF",organism=wh_organism,drop = dropcodes)$geneSim
	Sims[[5]][1]<-WangGeneSim(gene1,gene2,ont="CC",organism=wh_organism,drop = dropcodes)$geneSim

	Sims[[6]][1]<-TCSSGeneSim(gene1,gene2,ont="BP",organism=wh_organism,drop = dropcodes)$geneSim
	Sims[[7]][1]<-TCSSGeneSim(gene1,gene2,ont="MF",organism=wh_organism,drop = dropcodes)$geneSim
	Sims[[8]][1]<-TCSSGeneSim(gene1,gene2,ont="CC",organism=wh_organism,drop = dropcodes)$geneSim

	Sims[[9]][1]<-IntelliGOGeneSim(gene1,gene2,ont="BP",organism=wh_organism,drop = dropcodes)$geneSim
	Sims[[10]][1]<-IntelliGOGeneSim(gene1,gene2,ont="MF",organism=wh_organism,drop = dropcodes)$geneSim
	Sims[[11]][1]<-IntelliGOGeneSim(gene1,gene2,ont="CC",organism=wh_organism,drop = dropcodes)$geneSim

	Sims[[12]][1]<-KEGGSim(gene1,gene2)
	return(Sims)
}

GOKEGGSimsFromFile <- function(input,output="GOKEGGSims-ppiPre.csv",header=TRUE,sep=",", organism="yeast", drop ="IEA") ##KEGG- and GO-based similarity of protein pairs in an input file
{	
	cache<-read.csv(file=input,header=header,sep=sep)
	wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine","canine","anopheles","ecsakai","chicken","chimp","malaria","rhesus","pig","xenopus", "coelicolor"))

	dropcodes <- drop
	SimsFromFile<-data.frame(protein1=cache[1],protein2=cache[2],BPWang=0,MFWang=0,CCWang=0,BPTCSS=0,MFTCSS=0,CCTCSS=0,BPIG=0,MFIG=0,CCIG=0,KEGGSim=0)
	i<-1
	for(n in 1:length(cache[[1]]))
	{
		print(paste("Computing GO- & KEGG-based similarities of",as.character(cache[[1]][i]),"and", as.character(cache[[2]][i])))
		SimsFromFile[[3]][i]<-WangGeneSim(as.character(cache[[1]][i]),as.character(cache[[2]][i]),ont="BP",organism=wh_organism,drop = dropcodes )$geneSim
     	 	SimsFromFile[[4]][i]<-WangGeneSim(as.character(cache[[1]][i]),as.character(cache[[2]][i]),ont="MF",organism=wh_organism,drop = dropcodes )$geneSim
      		SimsFromFile[[5]][i]<-WangGeneSim(as.character(cache[[1]][i]),as.character(cache[[2]][i]),ont="CC",organism=wh_organism,drop = dropcodes )$geneSim
		SimsFromFile[[6]][i]<-TCSSGeneSim(as.character(cache[[1]][i]),as.character(cache[[2]][i]),ont="BP",organism=wh_organism,drop = dropcodes )$geneSim
		SimsFromFile[[7]][i]<-TCSSGeneSim(as.character(cache[[1]][i]),as.character(cache[[2]][i]),ont="MF",organism=wh_organism,drop = dropcodes )$geneSim
		SimsFromFile[[8]][i]<-TCSSGeneSim(as.character(cache[[1]][i]),as.character(cache[[2]][i]),ont="CC",organism=wh_organism,drop = dropcodes )$geneSim
		SimsFromFile[[9]][i]<-IntelliGOGeneSim(as.character(cache[[1]][i]),as.character(cache[[2]][i]),ont="BP",organism=wh_organism,drop = dropcodes )$geneSim
		SimsFromFile[[10]][i]<-IntelliGOGeneSim(as.character(cache[[1]][i]),as.character(cache[[2]][i]),ont="MF",organism=wh_organism,drop = dropcodes )$geneSim
		SimsFromFile[[11]][i]<-IntelliGOGeneSim(as.character(cache[[1]][i]),as.character(cache[[2]][i]),ont="CC",organism=wh_organism,drop = dropcodes )$geneSim
		SimsFromFile[[12]][i]<-KEGGSim(as.character(cache[[1]][i]),as.character(cache[[2]][i]))	
		i <- i+1
	}
	write.csv(SimsFromFile,file=output,row.names=FALSE)
}

`WangGeneSim` <-
function(gene1, gene2, ont="MF", organism="yeast", drop="IEA"){
	wh_ont <- match.arg(ont, c("MF", "BP", "CC")) 
	wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine","canine", "anopheles", "ecsakai", "chicken", "chimp", "malaria", "rhesus", "pig", "xenopus", "coelicolor"))


	go1 <- GetOntology(gene1, organism= wh_organism, ontology= wh_ont, dropCodes=drop) 
	go2 <- GetOntology(gene2, organism= wh_organism, ontology= wh_ont, dropCodes=drop)
	if (sum(!is.na(go1)) == 0 || sum(!is.na(go2)) == 0) {
		return (list(geneSim=NA, GO1=go1, GO2=go2)) 
	}

	go1 <- unlist(go1)
	go2 <- unlist(go2)
	m <- length(go1)
	n <- length(go2)
	 
	scores <- matrix(nrow=m, ncol=n)
	rownames(scores) <- go1
	colnames(scores) <- go2

	for( i in 1:m) {
		for (j in 1:n) {
			scores[i,j] <- WangGoSim(go1[i], go2[j], wh_ont, wh_organism)
		}
	}
#原始	if (!sum(!is.na(scores))) return (NA)	
       if (!sum(!is.na(scores))) return (list(geneSim=NA, GO1=go1, GO2=go2)) 

	if (n ==1 || m == 1) {
#原始		return (max(scores))
		return (list(geneSim=max(scores), GO1=go1, GO2=go2)) 
	}	
	sim <- (sum(sapply(1:m, function(x) {max(scores[x,], na.rm=TRUE)})) + sum(sapply(1:n, function(x) {max(scores[,x], na.rm=TRUE)})))/(m+n)	

	sim <- round(sim, digits=3)
	return (list(geneSim=sim, GO1=go1, GO2=go2)) 
}

`WangGoSim` <-
function(GOID1, GOID2, ont="MF", organism="yeast"){
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine", "canine", "anopheles", "ecsakai", "chicken", "chimp", "malaria", "rhesus", "pig", "xenopus", "coelicolor"))

	sim <- WangMethod(GOID1, GOID2, ont=wh_ont, wh_organism) 
	sim <- unname(sim, force=TRUE)
	return(round(sim, digits=3))
}

WangMethod <- function(GOID1, GOID2, ont="MF", organism="yeast") {
	if(!exists("ppiPreEnv")) .initial()
	weight.isa = 0.8
	weight.partof = 0.6

	if (GOID1 == GOID2)
		return (1)		

	Parents.name <- switch(ont,
		MF = "MFParents",
		BP = "BPParents",
		CC = "CCParents"
	)
	if (!exists(Parents.name, envir=ppiPreEnv)) { 
		GetGOParents(ont)
	}
	Parents <- get(Parents.name, envir=ppiPreEnv) 
	
	sv.a <- 1
	sv.b <- 1
	sw <- 1
	names(sv.a) <- GOID1
	names(sv.b) <- GOID2 
	
	sv.a <- WangSemVal(GOID1, ont, Parents, sv.a, sw, weight.isa, weight.partof) 
	sv.b <- WangSemVal(GOID2, ont, Parents, sv.b, sw, weight.isa, weight.partof)
	
	sv.a <- uniqsv(sv.a)
	sv.b <- uniqsv(sv.b)
	
	idx <- intersect(names(sv.a), names(sv.b))
	inter.sva <- unlist(sv.a[idx])
	inter.svb <- unlist(sv.b[idx])
	sim <- sum(inter.sva,inter.svb) / sum(sv.a, sv.b)
	return(sim)
}
WangSemVal <- function(goid, ont, Parents, sv, w, weight.isa, weight.partof) {
	if(!exists("ppiPreCache")) 
		return(WangSemVal_internal(goid, ont, Parents, sv, w, weight.isa, weight.partof))
	goid.ont <- paste(goid, ont, sep=".")
	if (!exists(goid.ont, envir=ppiPreCache)) {
	  	value <- WangSemVal_internal(goid, ont, Parents, sv, w, weight.isa, weight.partof) 
	  	assign(eval(goid.ont), value, envir=ppiPreCache)
	}
	return(get(goid.ont, envir=ppiPreCache))
}

WangSemVal_internal <- function(goid, ont, Parents, sv, w, weight.isa, weight.partof) {
	p <- Parents[goid] 
	p <- unlist(p[[1]]) 
	if (length(p) == 0) {
		return(0)
	}
	relations <- names(p)   
	old.w <- w
	for (i in 1:length(p)) {
		if (relations[i] == "is_a") {
			w <- old.w * weight.isa
		} else {
			w <- old.w * weight.partof
		}
		names(w) <- p[i]   
		sv <- c(sv,w)
		if (p[i] != "all") {
			sv <- WangSemVal_internal(p[i], ont, Parents, sv, w, weight.isa, weight.partof)
		}
	}
	return (sv)
}

uniqsv <- function(sv) { 
	sv <- unlist(sv)
	una <- unique(names(sv))
	sv <- unlist(sapply(una, function(x) {max(sv[names(sv)==x])}))
	return (sv)
}


################
TCSSGetChildren <- function(ont="MF") {
	if(!exists("ppiPreEnv")) .initial()

	wh_Children <- switch(ont,
		MF = "MFChildren",
		BP = "BPChildren",
		CC = "CCChildren"
	)	
		
	Children <- switch(ont,
		MF = AnnotationDbi::as.list(GOMFCHILDREN) ,
		BP = AnnotationDbi::as.list(GOBPCHILDREN) , 
		CC = AnnotationDbi::as.list(GOCCCHILDREN)	
	)
	assign(eval(wh_Children), Children, envir=ppiPreEnv)
}

TCSSGetAncestors <- function(ont="MF") {
	if(!exists("ppiPreEnv")) .initial()
	wh_Ancestors <- switch(ont,
		MF = "MFAncestors",
		BP = "BPAncestors",
		CC = "CCAncestors"	
	)		
	Ancestors <- switch(ont,
		MF = AnnotationDbi::as.list(GOMFANCESTOR) ,
		BP = AnnotationDbi::as.list(GOBPANCESTOR) , 
		CC = AnnotationDbi::as.list(GOCCANCESTOR)	
	)
	assign(eval(wh_Ancestors), Ancestors, envir=ppiPreEnv)
}

TCSSGetOffsprings <- function(ont="MF") {
	if(!exists("ppiPreEnv")) .initial()
	wh_Offsprings <- switch(ont,
		MF = "MFOffsprings",
		BP = "BPOffsprings",
		CC = "CCOffsprings"	
	)		
	Offsprings <- switch(ont,
		MF = AnnotationDbi::as.list(GOMFOFFSPRING) ,
		BP = AnnotationDbi::as.list(GOBPOFFSPRING) , 
		CC = AnnotationDbi::as.list(GOCCOFFSPRING)	
	)
	assign(eval(wh_Offsprings), Offsprings, envir=ppiPreEnv)
}

CheckAnnotationPackage <- function(species){
# 	pkgname <- switch (species,
# 		human = "org.Hs.eg.db",
# 		fly = "org.Dm.eg.db",
# 		mouse = "org.Mm.eg.db",
# 		rat = "org.Rn.eg.db",
# 		yeast = "org.Sc.sgd.db",
# 		zebrafish = "org.Dr.eg.db",
# 		worm = "org.Ce.eg.db",
# 		arabidopsis = "org.At.tair.db",
# 		ecolik12 = "org.EcK12.eg.db", 
# 		bovine	= "org.Bt.eg.db",
# 		canine	= "org.Cf.eg.db", 
# 		anopheles	=	"org.Ag.eg.db", 
# 		ecsakai	=	"org.EcSakai.eg.db", 
# 		chicken	=	"org.Gg.eg.db", 
# 		chimp	=	"org.Pt.eg.db", 
# 		malaria	=	"org.Pf.plasmo.db", 
# 		rhesus	=	"org.Mmu.eg.db", 
# 		pig	= 	"org.Ss.eg.db", 
# 		xenopus	=	"org.Xl.eg.db",
# 		coelicolor	=	"org.Sco.eg.db"
# 	)
	if (species == "human")
		if(!require(org.Hs.eg.db))
			stop("The package org.Hs.eg.db is needed.")
 	if (species == "yeast")
		if(!require(org.Sc.sgd.db))
			stop("The package org.Sc.sgd.db is needed.")
 	if (species == "fly")
		if(!require(org.Dm.eg.db))
			stop("The package org.Dm.eg.db is needed.")
 	if (species == "mouse")
		if(!require(org.Mm.eg.db))
			stop("The package org.Mm.eg.db is needed.")
 	if (species == "rat")
		if(!require(org.Rn.eg.db))
			stop("The package org.Rn.eg.db is needed.")
 	if (species == "zebrafish")
		if(!require(org.Dr.eg.db))
			stop("The package org.Dr.eg.db is needed.")
 	if (species == "worm")
		if(!require(org.Ce.eg.db))
			stop("The package org.Ce.eg.db is needed.")
 	if (species == "arabidopsis")
		if(!require(org.At.tair.db))
			stop("The package org.At.tair.db is needed.")
 	if (species == "ecolik12")
		if(!require(org.EcK12.eg.db))
			stop("The package org.EcK12.eg.db is needed.")
 	if (species == "bovine")
		if(!require(org.Bt.eg.db))
			stop("The package org.Bt.eg.db is needed.")
 	if (species == "canine")
		if(!require(org.Cf.eg.db))
			stop("The package org.Cf.eg.db is needed.")
 	if (species == "anopheles")
		if(!require(org.Ag.eg.db))
			stop("The package org.Ag.eg.db is needed.")
 	if (species == "ecsakai")
		if(!require(org.EcSakai.eg.db))
			stop("The package org.EcSakai.eg.db is needed.")
 	if (species == "chicken")
		if(!require(org.Gg.eg.db))
			stop("The package org.Gg.eg.db is needed.")
 	if (species == "chimp")
		if(!require(org.Pt.eg.db))
			stop("The package org.Pt.eg.db is needed.")
 	if (species == "malaria")
		if(!require(org.Pf.plasmo.db))
			stop("The package org.Pf.plasmo.db is needed.")
 	if (species == "rhesus")
		if(!require(org.Mmu.eg.db))
			stop("The package org.Mmu.eg.db is needed.")
 	if (species == "pig")
		if(!require(org.Ss.eg.db))
			stop("The package org.Ss.eg.db is needed.")
 	if (species == "xenopus")
		if(!require(org.Xl.eg.db))
			stop("The package org.Xl.eg.db is needed.")
 	if (species == "coelicolor")
		if(!require(org.Sco.eg.db))
			stop("The package org.Sco.eg.db is needed.")
}

GetGOMap <- function(organism="yeast") {
	if(!exists("ppiPreEnv")) .initial()
	CheckAnnotationPackage(organism) #download and install the packages
	species <- switch(organism,
		human = "Hs",
		fly = "Dm",
		mouse = "Mm",
		rat = "Rn",
		yeast = "Sc",
		zebrafish = "Dr",
		worm = "Ce",
		arabidopsis = "At",
		ecolik12 = "EcK12",
		bovine	= "Bt",
		canine	= "Cf", 
		anopheles	=	"Ag", 
		ecsakai	=	"EcSakai", 
		chicken	=	"Gg", 
		chimp	=	"Pt", 
		malaria	=	"Pf", 
		rhesus	=	"Mmu", 
		pig	= "Ss", 
		xenopus	=	"Xl",
		coelicolor	=	"Sco"
	)

	gomap <- switch(organism,
		human = org.Hs.egGO,
		fly = org.Dm.egGO,
		mouse = org.Mm.egGO,
		rat = org.Rn.egGO,
		yeast = org.Sc.sgdGO,
		zebrafish = org.Dr.egGO,
		worm = org.Ce.egGO,
		arabidopsis = org.At.tairGO,
		ecoli = org.EcK12.egGO,
		bovine	= org.Bt.egGO,
		canine	= org.Cf.egGO, 
		anopheles	=	org.Ag.egGO, 
		ecsakai	=	org.EcSakai.egGO, 
		chicken	=	org.Gg.egGO, 
		chimp	=	org.Pt.egGO, 
		malaria	=	org.Pf.plasmoGO, 
		rhesus	=	org.Mmu.egGO, 
		pig	= org.Ss.egGO, 
		xenopus	=	org.Xl.egGO,	
		coelicolor	=	org.Sco.egGO
	)

	assign(eval(species), gomap, envir=ppiPreEnv) 
}

`GetOntology` <-  function(gene, organism, ontology, dropCodes) {
	.initial() 
	species <- switch(organism,
		human = "Hs",
		fly = "Dm",
		mouse = "Mm",
		rat = "Rn",
		yeast = "Sc",
		zebrafish = "Dr",
		worm = "Ce",
		arabidopsis = "At",
		ecolik12 = "EcK12",
		bovine	= "Bt",
		canine	= "Cf", 
		anopheles	=	"Ag", 
		ecsakai	=	"EcSakai", 
		chicken	=	"Gg", 
		chimp	=	"Pt", 
		malaria	=	"Pf", 
		rhesus	=	"Mmu", 
		pig	= "Ss", 
		xenopus	=	"Xl",
		coelicolor	=	"Sco"
	)

	if (!exists(species, envir=ppiPreEnv)) {
		GetGOMap(organism)
	}
	gomap <- get(species, envir=ppiPreEnv) 

    	allGO <- gomap[[gene]] 
	if (is.null(allGO)) {
    		return (NA)
   	}
    	if (sum(!is.na(allGO)) == 0) {
    		return (NA)
   	}
    	if(!is.null(dropCodes)) { 
     	 	evidence<-sapply(allGO, function(x) x$Evidence) 
      		drop<-evidence %in% dropCodes 
      		allGO<-allGO[!drop] 
    	}

    	category<-sapply(allGO, function(x) x$Ontology) 
   	allGO<-allGO[category %in% ontology] 

   	if(length(allGO)==0) return (NA)
  	return (unlist(unique(names(allGO)))) #return the GOIDs
}

TCSSComputeIC <- function(dropCodes="IEA", ont, organism) {
print("Calulating IC...")
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine","canine","anopheles","ecsakai","chicken","chimp","malaria","rhesus","pig","xenopus", "coelicolor"))
	CheckAnnotationPackage(wh_organism)
	gomap <- switch(organism,
		human = org.Hs.egGO,
		fly = org.Dm.egGO,
		mouse = org.Mm.egGO,
		rat = org.Rn.egGO,
		yeast = org.Sc.sgdGO,
		zebrafish = org.Dr.egGO,
		worm = org.Ce.egGO,
		arabidopsis = org.At.tairGO,
		ecoli = org.EcK12.egGO,
		bovine	= org.Bt.egGO,
		canine	= org.Cf.egGO, 
		anopheles	=	org.Ag.egGO, 
		ecsakai	=	org.EcSakai.egGO, 
		chicken	=	org.Gg.egGO, 
		chimp	=	org.Pt.egGO, 
		malaria	=	org.Pf.plasmoGO, 
		rhesus	=	org.Mmu.egGO, 
		pig	= org.Ss.egGO, 
		xenopus	=	org.Xl.egGO,		
		coelicolor	=	org.Sco.egGO
	)

	mapped_genes <- mappedkeys(gomap)
	gomap = AnnotationDbi::as.list(gomap[mapped_genes])
	if (!is.null(dropCodes)){
		gomap<-sapply(gomap,function(x) sapply(x,function(y) c(y$Evidence %in% dropCodes, y$Ontology %in% wh_ont)))
		gomap<-sapply(gomap, function(x) x[2,x[1,]=="FALSE"])
		gomap<-gomap[sapply(gomap,length) >0]		
	}else {
		gomap <- sapply(gomap,function(x) sapply(x,function(y) y$Ontology %in% wh_ont))
	}
	
	goterms<-unlist(sapply(gomap, function(x) names(x)), use.names=FALSE) # all GO terms appearing in an annotation	
	goids <- toTable(GOTERM)

	goids <- unique(goids[goids[,"Ontology"] == wh_ont, "go_id"])  	
	gocount <- table(goterms)
	goname <- names(gocount) #goid of specific organism and selected category.

	go.diff <- setdiff(goids, goname)
	m <- double(length(go.diff)) 
	names(m) <- go.diff
	gocount <- as.vector(gocount)
	names(gocount) <- goname
	gocount <- c(gocount, m)

	Offsprings.name <- switch(wh_ont,
		MF = "MFOffsprings",
		BP = "BPOffsprings",
		CC = "CCOffsprings"	
	)	
	if (!exists(Offsprings.name, envir=ppiPreEnv)) {
		TCSSGetOffsprings(wh_ont)
	}
	Offsprings <- get(Offsprings.name, envir=ppiPreEnv)	
	cnt <- sapply(goids,function(x){ c=gocount[unlist(Offsprings[x])]; gocount[x]+sum(c[!is.na(c)])})		
	names(cnt) <- goids	
	IC<- -log(cnt/sum(gocount))
print("done...")		
	return (IC)
}
rebuildICdata <- function(){
	ont <- c("MF","CC", "BP")
	species <- c("human", "rat", "mouse", "fly", "yeast", "zebrafish", "arabidopsis","worm", "ecolik12", "bovine","canine","anopheles","ecsakai","chicken","chimp","malaria","rhesus","pig","xenopus","coelicolor") 
	cat("------------------------------------\n")
	cat("calulating Information Content...\nSpecies:\t\tOntology\n")
	for (i in ont) {
		for (j in species) {
			cat(j)
			cat("\t\t\t")
			cat(i)
			cat("\n")
			TCSSComputeIC(ont=i, organism=j)
		}
	}
	cat("------------------------------------\n")
	print("done...")
}

GetLatestCommonAncestor<-function(GOID1, GOID2, ont, organism){
#print("Calulating Latest Common Ancestor...")
	if(!exists("ppiPreEnv")) .initial()
	
	fname <- paste("Info_Contents", ont, organism, sep="_")
	tryCatch(utils::data(list=fname, package="ppiPre", envir=ppiPreEnv))	
	InfoContents <- get("IC", envir=ppiPreEnv)

	rootCount <- max(InfoContents[InfoContents != Inf])
	InfoContents["all"] = 0
	p1 <- InfoContents[GOID1]/rootCount
	p2 <- InfoContents[GOID2]/rootCount    
	if(is.na(p1) || is.na(p2)) return (NA)
	if (p1 == 0 || p2 == 0) return (NA)
	Ancestor.name <- switch(ont,MF = "MFAncestors",BP = "BPAncestors",CC = "CCAncestors")	
	if (!exists(Ancestor.name, envir=ppiPreEnv)) {
		TCSSGetAncestors(ont)
	}

	Ancestor <- get(Ancestor.name, envir=ppiPreEnv)					
	ancestor1 <- unlist(Ancestor[GOID1])
	ancestor2 <- unlist(Ancestor[GOID2])
	if (GOID1 == GOID2) { 
		commonAncestor <- GOID1
	} else if (GOID1 %in% ancestor2) { 
		commonAncestor <- GOID1
	} else if (GOID2 %in% ancestor1) {
		commonAncestor <- GOID2
	} else { 
		commonAncestor <- intersect(ancestor1, ancestor2)
	}
	if (length(commonAncestor) == 0)
		LCA<-NULL
	max<- -100
	LCA<-NULL
	for(a in commonAncestor){
		if(!is.na(InfoContents[a])) {
	  		if(InfoContents[a]>max){
				max<-InfoContents[a]
				LCA<-a
			}
		}
	}
#print("done...")
	return (LCA)

}

TCSSCompute_ICA<- function(dropCodes="IEA", ont, organism) {

  	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine","canine","anopheles","ecsakai","chicken","chimp","malaria","rhesus","pig","xenopus","coelicolor"))

	ICA.name <- paste(wh_organism,wh_ont,"ICA",sep="_")
	if (exists(ICA.name, envir=ppiPreEnv)) {
		return(get(ICA.name, envir=ppiPreEnv))
	}
	CheckAnnotationPackage(wh_organism)
#print("Calulating ICA...")
	gomap <- switch(wh_organism,
		human = org.Hs.egGO,
		fly = org.Dm.egGO,
		mouse = org.Mm.egGO,
		rat = org.Rn.egGO,
		yeast = org.Sc.sgdGO,
		zebrafish = org.Dr.egGO,
		worm = org.Ce.egGO,
		arabidopsis = org.At.tairGO,
		ecoli = org.EcK12.egGO,
		bovine	= org.Bt.egGO,
		canine	= org.Cf.egGO, 
		anopheles	=	org.Ag.egGO, 
		ecsakai	=	org.EcSakai.egGO, 
		chicken	=	org.Gg.egGO, 
		chimp	=	org.Pt.egGO, 
		malaria	=	org.Pf.plasmoGO, 
		rhesus	=	org.Mmu.egGO, 
		pig	= org.Ss.egGO, 
		xenopus	=	org.Xl.egGO,	
		coelicolor	=	org.Sco.egGO	
	)
	mapped_genes <- mappedkeys(gomap)
	gomap = AnnotationDbi::as.list(gomap[mapped_genes])
	if (!is.null(dropCodes)){
		gomap<-sapply(gomap,function(x) sapply(x,function(y) c(y$Evidence %in% dropCodes, y$Ontology %in% wh_ont)))
		gomap<-sapply(gomap, function(x) x[2,x[1,]=="FALSE"])
		gomap<-gomap[sapply(gomap,length) >0]		
	}else {
		gomap <- sapply(gomap,function(x) sapply(x,function(y) y$Ontology %in% wh_ont))
	}
	
	goterms<-unlist(sapply(gomap, function(x) names(x)), use.names=FALSE)	
	goids <- toTable(GOTERM)
	
	goids <- unique(goids[goids[,"Ontology"] == wh_ont, "go_id"])  	
	gocount <- table(goterms)
	goname <- names(gocount) 
	
	go.diff <- setdiff(goids, goname)
	m <- double(length(go.diff)) 
	names(m) <- go.diff
	gocount <- as.vector(gocount)
	names(gocount) <- goname
	gocount <- c(gocount, m)
	Children.name <- switch(wh_ont,
		MF = "MFChildren",
		BP = "BPChildren",
		CC = "CCChildren"	
	)	
	if (!exists(Children.name, envir=ppiPreEnv)) {
		TCSSGetChildren(wh_ont)
	}
	Children<- get(Children.name, envir=ppiPreEnv)	
	cnt <- sapply(goids,function(x){ c=gocount[unlist(Children[x])]; gocount[x]+sum(c[!is.na(c)])})
	cnt<-cnt+1;
	names(cnt)<-goids;
	ICA<- -log(cnt/sum(gocount))
#print("done...")
	assign(eval(ICA.name), ICA, envir=ppiPreEnv)
	return (ICA)
}	

`TCSSGoSim` <- function(GOID1, GOID2, ont, organism, ICA) {
	if(!exists("ppiPreEnv")) .initial()
	AnnotationIC <- ICA
	
	AnnotationIC["all"] = 0
	lca<-GetLatestCommonAncestor(GOID1, GOID2, ont, organism)
	if(length(lca)==0)
		sim<-0 
	ica_max<-max(AnnotationIC)
	ica_lca<-max(AnnotationIC[lca])
	sim<-ica_lca/ica_max
	return(round(sim, digits=3))
}
#m<-GetLatestCommonAncestor("GO:0030563", "GO:0004131", "MF", "yeast")
#n<-TCSSCompute_ICA(NULL,"MF","human")
#gosim<-TCSSGoSim("GO:0004022", "GO:0005515", "MF", "human")
#gosim1<-TCSSGoSim("GO:0043121", "GO:0019838", "MF", "human")
#go1 <- GetOntology("1019", organism= "human", ontology= "CC", dropCodes="IEA")
#go2 <- GetOntology("4831", organism= "human", ontology= "CC", dropCodes="IEA")

`TCSSGeneSim` <-
function(gene1, gene2, ont="MF", organism="yeast", drop="IEA"){

	wh_ont <- match.arg(ont, c("MF", "BP", "CC")) 
	wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine","canine","anopheles","ecsakai","chicken","chimp","malaria","rhesus","pig","xenopus","coelicolor"))

	go1 <- GetOntology(gene1, organism= wh_organism, ontology= wh_ont, dropCodes=drop) 
	go2 <- GetOntology(gene2, organism= wh_organism, ontology= wh_ont, dropCodes=drop)
	#if (sum(!is.na(go1)) == 0 || sum(!is.na(go2)) == 0) {
	if (is.na(go1)|| is.na(go2)) {
		return (list(geneSim=NA, GO1=go1, GO2=go2)) 
	}
	go1<-unlist(go1)
	go2<-unlist(go2)
	m<-length(go1)
	n<-length(go2)
	scores <- matrix(nrow=m, ncol=n)
	rownames(scores) <- go1
	colnames(scores) <- go2
	ICA<-TCSSCompute_ICA(drop,wh_ont,wh_organism) 
	
	for( i in 1:m) {
		for (j in 1:n) {
			scores[i,j] <-TCSSGoSim(go1[i], go2[j], wh_ont, wh_organism, ICA)
		}
	}
	if (!sum(!is.na(scores))) return (list(geneSim=NA, GO1=go1, GO2=go2)) 
	sim<-max(scores)
	sim <- round(sim, digits=3)

	return (list(geneSim=sim, GO1=go1, GO2=go2)) 
}	

#go1<-GetOntology("YPL094C", "yeast", "CC", "IEA")

#TCSSgene1<-TCSSGeneSim("1019","4831",ont="CC", organism="human")
#TCSSgene2<-TCSSGeneSim("241", "2561", ont="MF", organism="human")
#TCSSgene3<-TCSSGeneSim("snR18", "YPR062W", ont="MF", organism="yeast")
#TCSSgene2<-TCSSGeneSim("YPL070W", "YPR193C", ont="BP", organism="yeast")
#TCSSgene2<-TCSSGeneSim("YDR388W", "YNL094W", ont="MF", organism="yeast")
#TCSSgene2<-TCSSGeneSim("YDR388W", "YCR009C", ont="MF", organism="yeast")	
#geneSim("241", "2561", ont = "MF", organism = "human", measure = "Jiang")


GetGOParents <- function(ont="MF") {
	if(!exists("ppiPreEnv")) .initial()

	wh_Parents <- switch(ont,
		MF = "MFParents",
		BP = "BPParents",
		CC = "CCParents"	
	)#wh_Parents值为CCParents
		
	Parents <- switch(ont,
		MF = AnnotationDbi::as.list(GOMFPARENTS) ,
		BP = AnnotationDbi::as.list(GOBPPARENTS) , 
		CC = AnnotationDbi::as.list(GOCCPARENTS)	
	)
	assign(eval(wh_Parents), Parents, envir=ppiPreEnv)
}

IntelliGOInverseAnnotationFrequency<-function(dropCodes="IEA", goid,ont,organism){ 
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine","canine","anopheles","ecsakai","chicken","chimp","malaria","rhesus","pig","xenopus","coelicolor"))

	CheckAnnotationPackage(wh_organism)
	gomap <- switch(wh_organism,
		human = org.Hs.egGO,
		fly = org.Dm.egGO,
		mouse = org.Mm.egGO,
		rat = org.Rn.egGO,
		yeast = org.Sc.sgdGO,
		zebrafish = org.Dr.egGO,
		worm = org.Ce.egGO,
		arabidopsis = org.At.tairGO,
		ecoli = org.EcK12.egGO,
		bovine	= org.Bt.egGO,
		canine	= org.Cf.egGO, 
		anopheles	=	org.Ag.egGO, 
		ecsakai	=	org.EcSakai.egGO, 
		chicken	=	org.Gg.egGO, 
		chimp	=	org.Pt.egGO, 
		malaria	=	org.Pf.plasmoGO, 
		rhesus	=	org.Mmu.egGO, 
		pig	= org.Ss.egGO, 
		xenopus	=	org.Xl.egGO,
		coelicolor	=	org.Sco.egGO		
	)

	if (!is.null(dropCodes)){
		gomap<-sapply(gomap,function(x) sapply(x,function(y) c(y$Evidence %in% dropCodes, y$Ontology %in% wh_ont)))
		gomap<-sapply(gomap, function(x) x[2,x[1,]=="FALSE"])
		gomap<-gomap[sapply(gomap,length) >0]		
	}else {
		gomap <- sapply(gomap,function(x) sapply(x,function(y) y$Ontology %in% wh_ont))
	}

	Gti <- switch(wh_organism,
		human = length(org.Hs.egGO2EG[[goid]])	,
		fly = length(org.Dm.egGO2EG[[goid]]),
		yeast = length(org.Sc.sgdGO2ORF[[goid]]),
		worm = length(org.Ce.egGO2EG[[goid]]),	
		mouse = length(org.Mm.egGO2EG[[goid]]),
		rat = length(org.Rn.egGO2EG[[goid]]),
		zebrafish = length(org.Dr.egGO2EG[[goid]]),
		arabidopsis = length(org.At.tairGO2TAIR[[goid]]),
		ecoli = length(org.EcK12.egGO2EG[[goid]]),
		bovine	= length(org.Bt.egGO2EG[[goid]]),
		canine	= length(org.Cf.egGO2EG[[goid]]), 
		anopheles	=	length(org.Ag.egGO2EG[[goid]]), 
		ecsakai	=	length(org.EcSakai.egGO2EG[[goid]]), 
		chicken	=	length(org.Gg.egGO2EG[[goid]]), 
		chimp	=	length(org.Pt.egGO2EG[[goid]]), 
		malaria	=	length(org.Pf.plasmoGO2ORF[[goid]]), 
		rhesus	=	length(org.Mmu.egGO2EG[[goid]]), 
		pig	= length(org.Ss.egGO2EG[[goid]]), 
		xenopus	=	length(org.Xl.egGO2EG[[goid]]),
		coelicolor	=	length(org.Sco.egGO2EG[[goid]])
	)
	Gtot <- length(mappedkeys(gomap))
 	IAF <- log(Gtot/Gti)	
	return (IAF)
}

IntelliGOGetParents<-function(goid,ont,organsim){
	if(!exists("ppiPreEnv")) .initial()
	Parents.name <- switch(ont,
		MF = "MFParents",
		BP = "BPParents",
		CC = "CCParents"	
	)
	if (!exists(Parents.name, envir=ppiPreEnv)){ 
		GetGOParents(ont)
	}
	Parents <- get(Parents.name, envir=ppiPreEnv)
	p<- Parents[goid]
	par<-p[[1]][[1]]
	return (par)
}
IntelliGOGetDepth<-function(goid,ont,organsim){
	wh_ont<-ont
	depth<- 0
	parent<- IntelliGOGetParents(goid,wh_ont,organsim)
	if(length(parent)==0)
		return (depth+1)
	while(length(parent)!=0){
		depth<- depth+1
		parent<- IntelliGOGetParents(parent,wh_ont,organsim)
	}
	return (depth+1)
}
IntelliGOComputeEiEj<-function(GOID1,GOID2,ont,organsim){
	EiEj <- 0
	lca <- GetLatestCommonAncestor(GOID1,GOID2,ont,organsim)
  	depth_ei <- IntelliGOGetDepth(GOID1,ont,organsim)
	depth_ej <- IntelliGOGetDepth(GOID2,ont,organsim)
	depth_lca <- IntelliGOGetDepth(lca,ont,organsim)
	EiEj <- 2*depth_lca/(depth_ei+depth_ej)
	return (EiEj)
}
IntelliGOGOSim<-function(GOID1,GOID2,w1,w2,ont,organsim){
	sim<-0
	go1<-unlist(GOID1)
	go2<-unlist(GOID2)
	if(length(go1)==0||length(go2)==0)
		return (0)
	a<-IntelliGOComputeEiEj(go1,go2,ont,organsim)
	x<-sqrt(IntelliGOComputeEiEj(go1,go1,ont,organsim))
	y<-sqrt(IntelliGOComputeEiEj(go2,go2,ont,organsim))
	sim<- a/(x*y)
	return(round(sim, digits=3))
}

#pa<-IntelliGOGetParents("GO:0004022","MF","Human")
#depth_goid<-IntelliGOGetDepth("GO:0004022","MF","Human")

`IntelliGOGeneSim` <-
function(gene1, gene2, w1=1,w2=1,ont="MF", organism="yeast", drop="IEA"){
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine","canine","anopheles","ecsakai","chicken","chimp","malaria","rhesus","pig","xenopus","coelicolor"))

	go1 <- GetOntology(gene1, organism= wh_organism, ontology= wh_ont, dropCodes=drop) 
	go2 <- GetOntology(gene2, organism= wh_organism, ontology= wh_ont, dropCodes=drop)
	if (sum(!is.na(go1)) == 0 || sum(!is.na(go2)) == 0) {
		return (list(geneSim=NA, GO1=go1, GO2=go2)) 
	}
	weight_go1<-w1
	weight_go2<-w2
	go1<-unlist(go1)
	go2<-unlist(go2)
	m<-length(go1)
	n<-length(go2)
	scores <- matrix(nrow=m, ncol=n)
	rownames(scores) <- go1
	colnames(scores) <- go2

	for( i in 1:m ) {
		for (j in 1:n) {
			scores[i,j] <-IntelliGOGOSim(go1[i], go2[j],weight_go1,weight_go2,wh_ont, wh_organism)
		}
	}
	if (!sum(!is.na(scores))) return(list(geneSim=NA, GO1=go1, GO2=go2)) 
	if (n ==1 || m == 1) {
		sim<-max(scores)
	}
	
	sim <- (sum(sapply(1:m, function(x) {max(scores[x,], na.rm=TRUE)})) + sum(sapply(1:n, function(x) {max(scores[,x], na.rm=TRUE)})))/(m+n)
	sim <- round(sim, digits=3)
	return (list(geneSim=sim, GO1=go1, GO2=go2)) #返回值包括三部分：相似性、两个蛋白质对应的GO条目编号
}
	
#print(KEGGSim("YJL026W","YGR180C"))  #0.75
#iaf1<-IntelliGOInverseAnnotationFrequency(dropCodes="IEA","GO:0005515",ont = "MF", organism = "yeast")

#print(IntelliGOGeneSim("YOL001W", "YPL031C", w1=1,w2=1,ont="MF", organism="yeast")) #官网0.23，程序0.267（YOL001W（854161），YPL031C（856076））
#print(IntelliGOGeneSim("YNL037C", "YOR136W", w1=1,w2=1,ont="MF", organism="yeast")) #程序与官网均为1（YNL037C(855691)	YOR136W（854303））
#TCSSGeneSim("YOR065W", "YEL024W", ont="MF", organism="yeast", drop="IEA")
#print(WangGeneSim("YOR065W", "YEL024W", ont="CC", organism="yeast"))  #0.947
#print(PPSim("YOR065W", "YEL024W"))
#GOKEGGSims("YBL045C", "YPR191W", organism = "yeast", drop = "IEA")
#GOKEGGSims("YLR417W", "YPL002C", organism = "yeast", drop = "IEA")
#GOKEGGSims("YOR132W", "YOR069W", organism = "yeast", drop = "IEA")
#GOKEGGSims("YER133W", "YBR045C", organism = "yeast", drop = "IEA")