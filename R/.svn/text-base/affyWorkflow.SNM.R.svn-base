affyWorkflow.SNM <- function(rawDataDir, gene2row=NULL, bio.var=NULL, adj.var=NULL, int.var=NULL, rm.adj=TRUE, verbose=TRUE, logFile=TRUE, ...) {
	if(is.logical(logFile)) {
		logFile <- setUpLogFile(logFile)
	}	
	env <- new.env()
	# Load available platforms to catch situations where we don't have workflows devised for an input platform.
	platforms <- .loadAvailablePlatforms()
	cel.files <- .list.celfiles(rawDataDir)
	if(length(cel.files) == 0) {
		msg <- "No CEL files detected\n"
		sendLogMessage(msg,logFile)	
		cat(msg)
		class(msg) <- "try-error"
		return(msg)
	}
	cdfs <- sapply(paste(rawDataDir, cel.files, sep="/"), whatcdf)
# For each platform we implement the unsupervised affy workflow, summarize with EPSA, then create the mGenomics object,
# and finally add the output to the environment.
	sapply(unique(cdfs), function(x){
				if(x %in% platforms & sum(cdfs==x) > 1) {
					if(verbose) {
						cat("Platform:", x, "\n")
					}
					sendLogMessage(paste("Platform:", x),logFile)	
					
					if(verbose) {cat("Number of Samples:", sum(cdfs==x), "\n")}
					sendLogMessage(paste("Number of Samples:", sum(cdfs==x)),logFile)	
					if(x == "GenomeWideSNP_6"){
						tmp <- try(.runSNPCNV(rawDataDir, env=env, gene2row=gene2row, bio.var=bio.var, adj.var=adj.var, int.var=int.var, rm.adj=rm.adj, verbose=verbose, logFile=logFile, x=x, cdfs=cdfs, cel.files=cel.files, ...))
					}else{
						tmp <- try(.run3PrimeUtr(rawDataDir, env=env, gene2row=gene2row, bio.var=bio.var, adj.var=adj.var, int.var=int.var, rm.adj=rm.adj, verbose=verbose, logFile=logFile, x=x, cdfs=cdfs, cel.files=cel.files, ...))
					}
					if(class(tmp) == "try-error"){
						return(tmp)
					}
					string <- tmp[[1]]
					retval <- tmp[[2]]
					# Adds it to the environment
					assign(string,retval,envir=env)
					if(length(cdfs) > 1) { cat("---------------------------------------------------\n") }  
				}else if(sum(cdfs==x) < 2){
					retval <- paste("Too few cel files for", x,"platform.")
					cat("Too few cel files for platform", x,".\n")
					sendLogMessage(paste("Too few cel files for platform", x), logFile)
					if(length(cdfs) > 1) { cat("---------------------------------------------------\n") }  
					assign(x, retval, envir=env)
				}else{
					retval <- paste("Platform ", x, " not yet supported.")
					cat("Platform ", x, " not yet supported.\n")
					sendLogMessage(paste("Platform ", x, " not yet supported."), logFile)
					if(length(cdfs) > 1) { cat("---------------------------------------------------\n") }  
					assign(x, retval, envir=env)
				}
			}) -> eme
	as.list(env)
}

.loadAvailablePlatforms <- function(){ 
	res <- c("HG-Focus", "HG-U133A_2", "HG-U133A", "HG-U133B", "HG-U133_Plus_2", 
			"HG_U95Av2", "HG_U95B", "HG_U95C", "HG_U95D", "HG_U95E", "HT_HG-U133A", 
			"HT_HG-U133B", "HT_MG-430A", "HT_MG-430B", "HuEx-1_0-st-v2", "HuGene-1_0-st-v1", 
			"HuGene-1_1-st-v1", "MG_U74Av2", "MG_U74Bv2", "MG_U74Cv2", "MOE430A", "MOE430B", 
			"MoEx-1_0-st-v1", "MoGene-1_0-st-v1", "MoGene-1_1-st-v1", "Mouse430_2", 
			"Mouse430A_2", "RAE230A", "RAE230B", "RaEx-1_0-st-v1", "RaGene-1_0-st-v1", 
			"Rat230_2", "RG_U34A", "RG_U34B", "RG_U34C","GenomeWideSNP_6")
	res
}

.list.celfiles <- function (rawDataDir) 
{
	files <- list.files(rawDataDir)
	return(files[grep("\\.[cC][eE][lL]\\.gz$|\\.[cC][eE][lL]$", files)])
}

.runSNPCNV <- function(rawDataDir, gene2row, bio.var, adj.var, int.var, rm.adj, verbose, logFile, x, cdfs, cel.files, ...) {
	cel.files <- .list.celfiles(rawDataDir)
	platforms <- sapply(paste(rawDataDir,cel.files,sep="/"), affy:::whatcdf)
	cat("Loading CEL Files\n")
	sendLogMessage("Loading CEL Files",logFile)	
	scan.date <- cnames <- vector(mode = "character", length = length(cel.files))
	
	if(length(cel.files) > 100){
		bins <- seq(0,length(cel.files),by=100)
		if(bins[length(bins)] != length(cel.files)){
			bins <- c(bins, length(cel.files))
		}
		tDir <- tempdir()
		cat(1+bins[1], "\t", bins[2], "\n")
		abatch <- read.celfiles(filenames=cel.files[(1+bins[1]):bins[2]])
		scan.date[(1+bins[1]):bins[2]] <- runDate(abatch)
		cnames[(1+bins[1]):bins[2]] <- sampleNames(abatch)
		int <- oligo::intensity(abatch)
		data(list=c(annotation(abatch)))
		tmp.raw.int <- log2(int[cnv2row,]+1)
		raw.int <- matrix(0, nr=nrow(tmp.raw.int), nc=length(cel.files))
		rownames(raw.int) <- rownames(tmp.raw.int)
		raw.int[,(1+bins[1]):bins[2]] <- tmp.raw.int		
		for(i in 2:(length(bins)-1)){
			cat(1+bins[i], "\t", bins[i+1], "\n")
			abatch <- read.celfiles(filenames=cel.files[(1+bins[i]):bins[i+1]])
			scan.date[(1+bins[i]):bins[i+1]] <- runDate(abatch)
			cnames[(1+bins[i]):bins[i+1]] <- sampleNames(abatch)
			int <- oligo::intensity(abatch)
			raw.int[,(1+bins[i]):bins[i+1]] <- log2(int[cnv2row,]+1)		
			rm(int)
		}
		colnames(raw.int) <- cnames
	}else{
		abatch <- try(read.celfiles(filenames=paste(rawDataDir,cel.files,sep="/")))
		data(list=c(annotation(abatch)))
		scan.date <- runDate(abatch)
		int <- oligo::intensity(abatch)
		raw.int <- log2(int[cnv2row,]+1)
		cnames <- sampleNames(abatch)
	}
	if(is.null(int.var)) {
		int.var <- data.frame(array = factor(1:ncol(raw.int)))
	}
	sendLogMessage("Normalizing Data",logFile)	
	snm.fit <- tryCatch.W.E(snm(raw.int, bio.var, adj.var, int.var, diagnose=FALSE, rm.adj=rm.adj, verbose=FALSE),
			logFile,
			"snm",
			"Trying to normalize data with snm:")
	if(class(snm.fit) != "snm") {
		sendLogMessage("Could not fit model.  Please consult manual and try again.\n",logFile)
	}
	colnames(snm.fit$norm.dat) <- cnames
	dat <- snm.fit$norm.dat
	
	mGenomicsObject <- buildmGenomicsObject(dat, 
			annotation=annotation(abatch))
	mGenomicsObject[[1]]@protocolData@data$ScanDate <- as.character(unlist(scan.date))
	
	string <- annotation(abatch)
	list(string = string, 
			retval=mGenomicsObject)
}

.run3PrimeUtr <- function(rawDataDir, gene2row, bio.var, adj.var, int.var, rm.adj, verbose, logFile, x, cdfs, cel.files, ...) {
	# Load the data
	cat("Loading CEL Files\n")
	sendLogMessage("Loading CEL Files",logFile)	
	abatch <- tryCatch.W.E(ReadAffy(filenames=names(which(cdfs==x))), 
			logFile, "readAffy",
			"Trying to read data with ReadAffy")
	if(class(abatch) == "try-error"){
		sendLogMessage("Could not read data with ReadAffy.\n",logFile)
		return(abatch)
	}
	scan.date <- abatch@protocolData$ScanDate
	int <- exprs(abatch)
	# Gets the gene2row object that helps us define probe sets.
	gene2row.cdf <- getGene2Row(gene2row, annotation(abatch))
	# Loads the pm data and calls snm
	pms <- int[unlist(gene2row.cdf), ]
	pms[pms <= 1] <- 1
	data <- log2(pms)
	if(is.null(int.var)) {
		int.var <- data.frame(array = factor(1:ncol(data)))
	}
	# Calls snm
	if(verbose) {cat("Normalizing Data\n")}
	sendLogMessage("Normalizing Data",logFile)
	snm.fit <- tryCatch.W.E(snm(data, bio.var, adj.var, int.var, diagnose=FALSE, rm.adj=rm.adj,verbose=FALSE),
			logFile,
			"snm",
			"Trying to normalize data with snm:")
	if(class(snm.fit) != "snm") {
		sendLogMessage("Could not fit model.  Please consult manual and try again.\n",logFile)
	}
	colnames(snm.fit$norm.dat) <- sampleNames(abatch)
	gene2row.tmp <- split(1:nrow(pms), rep(names(gene2row.cdf), sapply(gene2row.cdf,length)))
	
	# Calls EPSA
	if(verbose) {cat("Summarizing Data\n")}
	sendLogMessage("Summarizing Data",logFile)	
	fits <- fit.pset(snm.fit$norm.dat, gene2row.tmp)
	
	# Builds mGenomics object
	if(verbose) {cat("Building mGenomics Object\n")}
	sendLogMessage("Building mGenomics object", logFile)
	stats <- list(singular.values=fits$singular.values,
			probe.weights=fits$probe.weights)
	retval <- buildmGenomicsObject(fits$estimated.rna.concentration,
			annotation=annotation(abatch),
			statistics=stats)
	retval$eset@protocolData$ScanDate <- abatch@protocolData$ScanDate
	string <- annotation(abatch)
	return(list(string=string,retval=retval))
}