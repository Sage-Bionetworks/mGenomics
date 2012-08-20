
detectTechnologies <- function(rawDataDir=NULL) {
	if(is.null(rawDataDir)){ 
		stop("Please provide a raw data directory\n")
	}
	all.txt.files <- list.files(path=rawDataDir)[grep(".txt",list.files())]
	platforms <- sapply(paste(rawDataDir, all.txt.files, sep=""), getAgilentPlatform)
	platforms <- platforms[which(platforms != "No Agilent Platforms Detected\n")]
	if(all(platforms=="No Agilent Platforms Detected\n"))  
	{ 
		agilent=FALSE
	}else{
		agilent=list(files=names(platforms))
	}
	cel.files <- list.celfiles(rawDataDir, full.names=T)
	if(length(cel.files) > 0){
		affy=list(files=cel.files)
	}else{
		affy=FALSE
	}
	list(agilent=agilent,affy=affy)
}

