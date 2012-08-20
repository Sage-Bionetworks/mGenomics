affyWorkflow <- function(rawDataDir,workflow,verbose=TRUE, logFile=TRUE, ...){ 
	if( ! any(workflow %in% c("rma","dchip","mas5","gcrma","frma","snm")) ){
		stop("workflow", workflow,"is currently NOT available.\n")
	}
	if(workflow=="snm"){
		return(affyWorkflow.SNM(rawDataDir, logFile=logFile, ...))
	}
	if(is.logical(logFile)) {
		logFile <- setUpLogFile(logFile)
	}	
	env <- new.env()
# Load available platforms to catch situations where we don't have workflows devised for an input platform.
	cel.files <- list.celfiles(rawDataDir, full.names=T, recursive=T)
	if(length(cel.files) == 0) {
		msg <- "No CEL files detected\n"
		sendLogMessage(msg,logFile)	
		cat(msg)
		class(msg) <- "try-error"
		return(msg)
	}
	cdfs <- sapply(cel.files, whatcdf)
# For each platform we implement the unsupervised affy workflow, summarize with EPSA, then create the mGenomics object,
# and finally add the output to the environment.
	sapply(unique(cdfs), function(x){
				if(sum(cdfs==x) >1){ 
					if(verbose) {
						cat("Platform:", x, "\n")
					}
					sendLogMessage(paste("Platform:", x),logFile)	
					if(verbose) {cat("Number of Samples:", sum(cdfs==x), "\n")}
					sendLogMessage(paste("Number of Samples:", sum(cdfs==x)),logFile)
					retval <- switch(workflow, 
							"rma"=.callRMAWorkflow(names(which(cdfs==x)), logFile),
							"dchip"=.callDchipWorkflow(names(which(cdfs==x)), logFile),
							"mas5"=.callMas5Workflow(names(which(cdfs==x)), logFile),
							"gcrma"=.callGcrmaWorkflow(names(which(cdfs==x)), logFile),
							"frozenRMA"=.callFrmaWorkflow(names(which(cdfs==x)), logFile))
					string <- annotation(retval)
					assign(string,retval,envir=env)
					if(length(cdfs) > 1) { cat("---------------------------------------------------\n") }  
				}else if(sum(cdfs==x) < 2){
							retval <- paste("Too few cel files for", x,"platform.")
							cat("Too few cel files for platform", x,".\n")
							sendLogMessage(paste("Too few cel files for platform", x), logFile)
							if(length(cdfs) > 1) { cat("---------------------------------------------------\n") }  
							assign(x, retval, envir=env)
			}
		})
  as.list(env)
}

.callRMAWorkflow <- function(filenames, logFile) {
	cat("Running justRMA\n")
	sendLogMessage("Running justRMA",logFile)	
	rmaResult <- justRMA(filenames=filenames)
	return(rmaResult)
}

.callDchipWorkflow <- function(filenames, logFile) {
	cat("Running dCHIP\n")
	sendLogMessage("Running dCHIP",logFile)
	abatch <- ReadAffy(filenames=filenames)
	dchipResult <- expresso (abatch, 
			normalize.method="invariantset", 
			bg.correct =FALSE, 
			pmcorrect.method="pmonly", 
			summary.method="liwong")
	return(dchipResult)
}

.callMas5Workflow <- function(filenames, logFile) {
	cat("Running MAS5\n")
	sendLogMessage("Running MAS5",logFile)	
	abatch <- ReadAffy(filenames=filenames)
	mas5Result<- mas5(abatch)
	return(mas5Result)
}

.callGcrmaWorkflow <- function(filenames, logFile) {
	cat("Running gcRMA\n")
	sendLogMessage("Running gcRMA",logFile)	
	abatch <- ReadAffy(filenames=filenames)
	gcRmaResult <- gcrma(abatch)
	return(gcRmaResult )
}

	