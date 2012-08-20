
callWorkflows <- function(rawDataDir, workflow=NULL, logFile=TRUE, ...) {

	# Call affy workflow
	res.affy <- affyWorkflow(rawDataDir, workflow=workflow, logFile=logFile, ...)
	# Call agilent workflow
	res.agilent <- agilentTxtWorkflow(rawDataDir, logFile=logFile, ...)
	
	if(class(res.affy) == "try-error" & class(res.agilent) == "try-error") {
		res <- "Could not detect Affymetrix CEL files or Agilent TXT files.  Currently these are the only two technologies that are supported\n"
		class(res) <- "try-error"
	}else if(class(res.affy) == "try-error" & class(res.agilent) != "try-error") {
		res <- res.agilent	
	}else if(class(res.affy) != "try-error" & class(res.agilent) == "try-error") {
		res <- res.affy	
	}else{
		res <- c(res.affy, res.agilent)
	}
	
	return(res)
}


