agilentTxtWorkflow <- function(rawDataDir, bio.var=NULL, adj.var=NULL, int.var=NULL, rm.adj=TRUE, verbose=TRUE, logFile=TRUE) {
	#
	# Function implements the agilent workflows.
	# First we check to see if the logfile has been passed in,
	# if not we set one up.
	#
	if(is.logical(logFile)) {
		logFile <- setUpLogFile(logFile)
	}
	
	#
	# Next we get all the txt files in the supplied raw data directory.
	# If there are no txt files then we return a message.  Otherwise,
	# we get the agilent platforms for each file.  
	# If we don't detect any agilent files amongst the txt files we 
	# return a message indicating this to the user.  Otherwise, we
	# attempt to process them
	#
	all.txt.files <-list.txtfiles(path=rawDataDir, full.names=T, recursive=T)
	if(length(all.txt.files) == 0){
		cat("No Agilent Platforms Detected\n")
		sendLogMessage("No Agilent platforms detected.",logFile)
		fits <- "No Agilent Platforms Detected\n"
		class(fits) <- "try-error"
		return(fits)
	}
	platforms <- sapply(all.txt.files, getAgilentPlatform)
	if(all(is.na(platforms))) {
		cat("No Agilent Platforms Detected\n")
		sendLogMessage("No Agilent platforms detected.",logFile)
		fits <- "No Agilent Platforms Detected\n"
		class(fits) <- "try-error"
		return(fits)
	}else{		
		agilent.files <- as.numeric(which(!is.na(platforms)))
		all.agilent.files <- all.txt.files[agilent.files]
		platforms <- platforms[agilent.files]
		platforms <- standardizePlatforms(platforms) # Maps the barcode ID to the Agilent platform name
		platform2file <- split(all.agilent.files, platforms)
		fits <- list()
		
		# For each platform we load the files, then normalize the data, then build the mGenomics object.
		# If there is only a single file for a given platform then we return an error message and do
		# not process it.  Otherwise we push ahead.
		for(i in 1:length(platform2file)) {
			if(length(platform2file[[i]]) > 1) { 
				cat("Platform:", names(platform2file)[i],"\n")			
				sendLogMessage(paste("Processing platform.", names(platform2file)[i]), logFile)
				files <- platform2file[[i]]
				cat("Number of Samples:", length(files),"\n")	
				sendLogMessage(paste("Number of samples.", length(files)), logFile)
				cat("Loading Samples\n")
				sendLogMessage("Loading samples.", logFile)
				# Load the files
				obj <- loadAgilentFiles(files, logFile)
				if(length(files) > 1){	
					cat("Normalizing Data\n")
					sendLogMessage("Normalizing data.", logFile)
					# Normalize the data
					platformObj <- normalizeAgilent(obj, logFile, names(platform2file)[i])
				}else{
					if(sum(!is.na(obj$red)) > 0) { is.red = TRUE }
					if(sum(!is.na(obj$green)) > 0) { is.green = TRUE }
					if(is.red & is.green){ data <- cbind(obj$green, obj$red)}
					if(is.red & !is.green){ data <- cbind(obj$red)}
					if(!is.red & is.green){ data <- cbind(obj$green)}
					platformObj <- buildmGenomicsObject(data, 
							annotation=names(platform2file)[i])
				}
				if(length(platform2file) > 1){ 
					cat("---------------------------------------------------\n")
				}
				if(class(platformObj) == "try-error"){
					fits[[i]] <- platformObj
					nms <- names(platform2file)[i]
					nms <- gsub("-",".",nms)
					names(fits)[i] <- nms
				}else{
					fits[[i]] <- platformObj
					nms <- names(platform2file)[i]
					nms <- gsub("-",".",nms)
					names(fits)[i] <- nms
				}
			}else{
				retval <- paste("Too few files for", names(platform2file)[i],"platform.")
				cat("Too few files for platform", names(platform2file)[i],".\n")
				sendLogMessage(paste("Too few files for platform", names(platform2file)[i]), logFile)
				if(length(platform2file) > 1) { cat("---------------------------------------------------\n") }  
				fits[[i]] <- retval
				nms <- names(platform2file)[i]
				nms <- gsub("-",".",nms)
				names(fits)[i] <- nms
			}
		}
	}
	fits
}

standardizePlatforms <- function(platforms){ 
	platformMap <- list()
	platformMap["1601197"] <- "Mouse_G4121A"
	platformMap["2511978"] <- "Mouse_G4121A"
	platformMap["1601239"] <- "Human_G4112A"
	platformMap["2512391"] <- "Human_G4112A"
	platformMap["1601485"] <- "Human_G4112F"
	platformMap["2514850"] <- "Human_G4112F"	
	platformMap["1601209"] <- "Human_G4110B"
	platformMap["2512097"] <- "Human_G4110B"
	platformMap["1601486"] <- "Mouse_G4112F"
	platformMap["2514868"] <- "Mouse_G4112F"
	platformMap["1601487"] <- "Rat_G4131F"
	platformMap["2514879"] <- "Rat_G4131F"
	platformMap["1601186"] <- "Rat_4130A"
	platformMap["2511868"] <- "Rat_4130A"
	platformMap["1601269"] <- "2512694"
	platformMap["1601332"] <- "2513326"
	platformMap["1601469"] <- "2514695"
	platformMap["113256"] <- "Human_G450"
	platformMap["113260"] <- "Human_G450"
	platformMap["113144"] <- "Human_G450"
	platformMap["113157"] <- "Human_G450"
	platformMap["105343"] <- "Human_G450"
	platformMap["2516436"] <- "miRNA_8x15K"
	platformMap["2521529"] <- "cgh-1x1m_g4447a"
	platformMap["2523364"] <- "hg-cgh-415k_g4124a"
	platformMap["2514693"] <- "Human_G4411B"
	platformMap["1601469"] <- "Human_G4411B"
	sapply(platforms, function(x){ 
				if(x %in% names(platformMap)) {
					platformMap[x]
				}else{
					x
				}
			}) -> ret
	names(ret) <- names(platforms)
	as.character(unlist(ret))
}

normalizeAgilent <- function(obj, logFile, annotation){ 
	is.green <- is.red <- FALSE;
	if(sum(!is.na(obj$red)) > 0) {
		red <- obj$red[,obj$valid.red]
		red[red < 1] <- 1
		red <- log2(red)
		is.red = TRUE
	}
	if(sum(!is.na(obj$green)) > 0) {
		green <- obj$green[,obj$valid.green]
		green[green < 1] <- 1
		green <- log2(green)
		is.green = TRUE
	}
	if(is.red & is.green) {
		if(ncol(green) != ncol(red)){ 
			sendLogMessage("Cannot fit unsupervised model due to unequal numbers of Cy5 and Cy3 channels across samples.  Most likely a result of NA values for one channel inside of a single file.\n  Please fix input and try again.", logFile)		
			ret <- "Cannot fit unsupervised model due to unequal numbers of Cy5 and Cy3 channels across samples.  Most likely a result of NA values for one channel inside of a single file.\n  Please fix input and try again.\n";
			class(ret) <- "try-error"
      return(ret)
		}
	}
	return.mat <- NULL
	if(is.green & is.red){ 
		array <- factor(rep(1:ncol(green), times=2))
		dye <- factor(rep(c("Cy3","Cy5"), each=ncol(green)))
		int.var <- data.frame(array=array, dye=dye)
		snm.fit <- tryCatch.W.E(snm(cbind(green,red), bio.var=NULL, adj.var=NULL, int.var=int.var),
				logFile,"snm",msg="Trying to normalize data with SNM.")
		if(class(snm.fit) != "snm"){
			sendLogMessage("Could not fit model.  Please consult manual and try again.\n",logFile)
		}else{
			return.mat <- matrix(NA, nr=nrow(obj$green), nc=ncol(obj$green)+ncol(obj$red))
			return.mat[,c(obj$valid.green,(obj$valid.red + ncol(obj$green)))] <- snm.fit$norm.dat
			colnames(return.mat) <- c(paste("Cy3_",basename(colnames(obj$green)),sep=""),
					paste("Cy5_",basename(colnames(obj$red)),sep=""))
			rownames(return.mat) <- rownames(obj$green)
		}
	}else if(is.green & !is.red){
		array <- factor(1:ncol(green))
		int.var <- data.frame(array=array)
		snm.fit <- tryCatch.W.E(snm(green, bio.var=NULL, adj.var=NULL, int.var=int.var),logFile,"snm",msg="Trying to normalize data with SNM.")
		if(class(snm.fit) != "snm"){
			sendLogMessage("Could not fit model.  Please consult manual and try again.\n",logFile)
		}else{
			return.mat <- matrix(NA, nr=nrow(obj$green), nc=ncol(obj$green))
			return.mat[,obj$valid.green] <- snm.fit$norm.dat
			colnames(return.mat) <- c(paste("Cy3_",basename(colnames(obj$green)),sep=""))
			rownames(return.mat) <- rownames(obj$green)
		}
	}else if(is.red){
		array <- factor(1:ncol(red))
		int.var <- data.frame(array=array)
		snm.fit <- tryCatch.W.E(snm(red, bio.var=NULL, adj.var=NULL, int.var=int.var),logFile,"snm",msg="Trying to normalize data with SNM.")
		if(class(snm.fit) != "snm"){
			sendLogMessage("Could not fit model.  Please consult manual and try again.\n",logFile)
		}else{
			return.mat <- matrix(NA, nr=nrow(red), nc=ncol(obj$red))
			return.mat[,obj$valid.red] <- snm.fit$norm.dat
			colnames(return.mat) <- c(paste("Cy5_",basename(colnames(obj$red)),sep=""))
			rownames(return.mat) <- rownames(obj$red)		
		}
	}
	gene2row <- split(1:nrow(return.mat), rownames(return.mat))
	t(sapply(gene2row, function(x){ 
			if(length(x) == 1){
				return.mat[x,]
			}else{
				colMeans(return.mat[x,])
			}
		})) -> ret.mat

	ret <- buildmGenomicsObject(ret.mat,
			annotation=annotation)
	return(ret)
}

removeNAs <- function(mat){ 
	# Removes columns with more than 0 NA values.
	ids <- which(colSums(is.na(mat)) == 0)
	return(mat[,ids])
}

loadAgilentFiles <- function(files, logFile, ...) {
	cat("Loading file 1 of",length(files))
	sendLogMessage(paste("Loading file 1 of",length(files)),logFile)

	#################################################################################
	# First we have to find a file to define the dimensions of the
	# data.  In some situations the first few files provided might 
	# contain issues that restrict our ability to load it into memory.
	# To account for this we look at each file until we find one that
	# can be read in.
	################################################################################
	for( i in 1:length(files)) {
		cat(i,"\n")
		data <- try(readAgilentFile(files[i], logFile), 
				silent=TRUE) # Read in the file
		# If we successfully read it in :
		if(class(data) != "try-error" & length(files) > 1) { 
			greenPValue <- redPValue <- 
					red <- green <- matrix(NA, 
							nr=nrow(data), 
							nc=length(files),
							dimnames=list(data$ProbeName, files))  # BOTS REQUEST
			is.greenPValues <- is.redPValues <- is.green <- is.red <- FALSE;
			if("gprocessedsignal" %in% tolower(colnames(data))){ 
				green[,i] <- data[,match("gprocessedsignal", tolower(colnames(data)))]
				is.green <- TRUE
			}
			if("rprocessedsignal" %in% tolower(colnames(data))){ 
				red[,i] <- data[,match("rprocessedsignal", tolower(colnames(data)))]
				is.red <- TRUE
			}
			if("gpvalfeateqbg" %in% tolower(colnames(data))) {
				greenPValue[,i] <- data[,match("gpvalfeateqbg", tolower(colnames(data)))]
				is.greenPValues <- TRUE
			}
			if("rpvalfeateqbg" %in% tolower(colnames(data))) {
				redPValue[,i] <- data[,match("rpvalfeateqbg", tolower(colnames(data)))]
				is.redPValues <- TRUE
			}
			non.controls <- try(which(data$ControlType == 0))
			if(class(non.controls) == "try-error") {
				cat("Warning: cannot find column ControlType, all probes treated as non-controls\n")
			}
			break			
		}else if(class(data) == "try-error" & length(files) > 1) {
			# if we get an error and there is more than file being input
			# then do nothing
		}else if(class(data) == "try-error" & length(files) == 1) {
			# if we get an error and there is only one file, then return NA
			return(NA)
		}else{
			return(NA)
		}
	}

	#
	# Now that we have the dimensions we can continue to load the 
  # files into memory.  
  #
	next.i <- i+1
	
	if(length(files) > 1){
#		for( i in next.i:length(files)) {
		for( i in 6:10) {
			cat("\rLoading file",i,"of",length(files))
			sendLogMessage(paste("Loading file",i,"of",length(files)),logFile)
			data <- try(readAgilentFile(files[i], logFile), silent=TRUE)
			if(class(data) == "try-error"){
				#	cat(data,"\n")
				green[,i] <- red[,i] <- greenPValue[,i] <- redPValue[,i] <- rep(NA, nrow(green))
			}
			else if(nrow(data) != nrow(green)) {
				cat("Incorrect # of probe features in", files[i],"\n")
				sendLogMessage(paste("Incorrect # of probe features in", files[i]),logFile)
				green[,i] <- red[,i] <- greenPValue[,i] <- redPValue[,i] <- rep(NA, nrow(green)) 
			}else{
				if("gprocessedsignal" %in% tolower(colnames(data))){ 
					green[,i] <- data[,match("gprocessedsignal", tolower(colnames(data)))]
				}
				if("rprocessedsignal" %in% tolower(colnames(data))){ 
					red[,i] <- data[,match("rprocessedsignal", tolower(colnames(data)))]
				}
				if("gpvalfeateqbg" %in% tolower(colnames(data))) {
					greenPValue[,i] <- data[,match("gpvalfeateqbg", tolower(colnames(data)))]
				}
				if("rpvalfeateqbg" %in% tolower(colnames(data))) {
					redPValue[,i] <- data[,match("rpvalfeateqbg", tolower(colnames(data)))]
				}
			}
		}
	}
	
	#
  # Split the input data into the various matrices 
  # and vectors of interest before returning successfully 
  # to the user.
  #
	if(is.green){ green <- green[non.controls,]; }
	if(is.red){ red <- red[non.controls,]; }
	if(is.greenPValues & is.redPValues) { 
		probePValues = cbind(greenPValue, redPValue)[non.controls,]
	}else if(is.greenPValues & !is.redPValues){ 
		probePValues = greenPValue[non.controls,]
	}else if(!is.greenPValues & is.redPValues){ 
		probePValues = redPValue[non.controls,]
	}else{
		probePValues = NA
	}
	valid.red <- which(colSums(is.na(red)) == 0)
	valid.green <- which(colSums(is.na(green)) == 0)
	obj <- list(red=red,
			green=green,
			valid.red=valid.red,
			valid.green=valid.green,
			probePValues=probePValues)
	cat("\n")
	obj
}

readAgilentFile <- function(file, logFile){
	perl.dir <- file.path(.path.package("mGenomics"),"Perl")
	perl.file <- file.path(perl.dir, "removeTrailing.pl") 
	if(grepl(".gz$",file,perl=TRUE)) {
		# The file is of type .gz
		uncompressedFileName <- tempfile()
		system(  paste("gunzip -c",file,">",uncompressedFileName )  )
		new.file <- tempfile()
		system(paste("perl",perl.file,uncompressedFileName,new.file))
	}else{
		# File is not compressed
		new.file <- tempfile()
		system(paste("perl",perl.file,file,new.file))
	}
	skip=0
	nrows=1
	colTypes <- NA	
	#Get columns of interest
	for(skip in skip:1000) {
		tmp <- read.table(new.file, sep="\t", quote="",nrows=nrows,skip=skip,stringsAsFactors=FALSE)
		if(tmp[which(!is.na(tmp))[1]] == "FEATURES") { 
			tmp <- tolower(tmp)
			sapply(tmp, function(x){
						if(x =="row") {"numeric"}
						else if(x == "col") {"numeric"}
						else if(x=="probename") {"character"}
						else if(x=="gprocessedsignal") {"numeric"}
						else if(x=="rprocessedsignal") {"numeric"}
						else if(x=="gpvalfeateqbg") {"numeric"}
						else if(x=="rpvalfeateqbg") {"numeric"}
						else if(x=="controltype") {"numeric"}
						else if(x=='logratio') {"numeric"}
						else if(x=='logratioerror') {"numeric"}
						else if(x=='pvaluelogratio') {"numeric"}
						else if(x=='gnumpix') {"numeric"}
						else if(x=='rnumpix') {"numeric"}
						else if(x=='gpixsdev') {"numeric"}
						else if(x=='rpixsdev') {"numeric"}
						else{ "NULL" }
					}) -> colTypes
			break
		}
	}
	if(is.na(colTypes[1])){ 
		# File did not contain standard agilent header
		sendLogMessage(paste("File",file,"did not contain standard agilent header"),logFile)
		var <- paste("File",file,"did not contain standard agilent header\n")
		class(var) <- "try-error"
		return(var)
	}
	data <- tryCatch.W.E(read.table(new.file,sep="\t", 
					quote="",
					skip=skip,
					stringsAsFactors=FALSE,
					colClasses=colTypes,
					header=TRUE),logFile,NULL,msg=paste("Trying to load file: ",file))
	system(paste("rm",new.file))
	if(class(data) == "try-error") {
		var <- paste("Error reading ", file, ":",data,"\n")
		sendLogMessage(var,logFile)
		class(var) <- "try-error"
		return(var)
	}else{	
		return(data)
	}
}


agilentColumnsOfInterest <- function(){
	interest <- c("Row","Col","ProbeName","gProcessedSignal", "rProcessedSignal", "gPValFeatEqBG","rPValFeatEqBG","ControlType")
	hmm <- list("numeric","numeric","character","numeric","numeric","numeric","numeric","numeric")
	names(hmm) <- tolower(interest)
}

getAgilentPlatform <- function(file){
	# Identifies the agilent file and removes the trailing blank values from each row using
	# a perl script.
	skip=0
	nrows=1
	platform <- NA
	#Test if file is a raw agilent data file
	header <- scan(file,sep="",n=20,na.strings=c("NA",""), what="character",quiet=TRUE)
	if(length(header) > 0 & header[	which(header!= "")[1]] == "TYPE") {
		for(skip in skip:10) { 	  
			tmp <- read.table(file,sep="\t", quote="",nrows=nrows,skip=skip,stringsAsFactors=FALSE)
			if(tmp[which(!is.na(tmp))[1]] == "FEPARAMS") { 
#				ids <- which(tmp == "FeatureExtractor_PatternName")
				ids <- which(tmp == "FeatureExtractor_Barcode")				
				tmp2 <- read.table(file,sep="\t", quote="",nrows=nrows,skip=skip+1,stringsAsFactors=FALSE)
				platform <- as.character(tmp2[ids])
				platform <- paste(strsplit(platform,"")[[1]][1:7],collapse="")
				break
			}	
		}
	}
	if(length(platform) == 0 | is.na(platform)) {
#		cat("No Agilent Platforms Detected\n")
		platform <- NA
	}
	platform
}


list.txtfiles <- function (...) 
{
	files <- list.files(...)
	return(files[grep("\\.[tT][xX][tT]\\.gz$|\\.[tT][xX][tT]$", 
							files)])
}
