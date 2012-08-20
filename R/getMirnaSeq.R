getMirnaSeq <- function(rawDataDir, platform, metadata){
	if(is.null(platform)){
		var <- "Must provide valid platform.\n"
		class(var) <- 'try-error'
		return(var)
	}
	# Get the miRNA files to load
	tmp <- list.files(rawDataDir, full.names=T, recursive=T)
	mirnaFiles <- tmp[grep("mirna.quantification.txt", basename(tmp), fixed=T)]
	
	# Read in the counts for both the read counts and counts per million.
	# First we load them into a list. Then we merge them into a matrix.
  # We do this in case the number of miRNAs varies across files.
	readCountsPerMillion <- readCounts <- vector("list", length(mirnaFiles))
	for(i in 1:length(mirnaFiles)){
		cat("\rloading file ", i, " ", basename(mirnaFiles[i]))
		tmpMat <- read.delim(file=mirnaFiles[i], as.is=T, header=T, row.names=NULL, quote="")
		readCounts[[i]] <- tmpMat$read_count
		readCountsPerMillion[[i]] <- tmpMat$reads_per_million_miRNA_mapped
		names(readCountsPerMillion[[i]]) <- names(readCounts[[i]]) <- tmpMat$miRNA_ID
		if(i == length(mirnaFiles)){
			cat("\n")
		}
	}
	# Now merge them into a matrix
	allNmsReadCounts  <- unique(names(unlist(readCounts)))
	matReadCounts <- matrix(NA, 
			nr=length(allNmsReadCounts), 
			nc=length(mirnaFiles), 
			dimnames=list(allNmsReadCounts,basename(mirnaFiles)))
	for(i in 1:length(mirnaFiles)){
		matReadCounts[names(readCounts[[i]]),i] <- readCounts[[i]]
	}	
	allNmsReadCountsPerMillion  <- unique(names(unlist(readCountsPerMillion)))
	matReadCountsPerMillion <- matrix(NA, 
			nr=length(allNmsReadCountsPerMillion), 
			nc=length(mirnaFiles), 
			dimnames=list(allNmsReadCountsPerMillion,basename(mirnaFiles)))
	for(i in 1:length(mirnaFiles)){
		matReadCountsPerMillion[names(readCountsPerMillion[[i]]),i] <- readCountsPerMillion[[i]]
	}	
	
	matReadCounts <- .reduceToAliquots(matReadCounts, metadata)
	matReadCountsPerMillion <- .reduceToAliquots(matReadCountsPerMillion, metadata)
	
	mirnaEsetReadCounts <- new("ExpressionSet", exprs=matReadCounts, annotation=platform)
	mirnaEsetReadCountsPerMillion <- new("ExpressionSet", exprs=matReadCountsPerMillion, annotation=platform)
	
	isoformFiles <- tmp[grep("isoform.quantification.txt", basename(tmp), fixed=T)]
	readCountsPerMillion <- readCounts <- vector("list", length(isoformFiles))
	
	for(i in 1:length(isoformFiles)){
		cat("\rloading file ",i," ", basename(isoformFiles[i]))
		tmpMat <- read.delim(file=isoformFiles[i], as.is=T, header=T, row.names=NULL, quote="")
		readCounts[[i]] <- tmpMat$read_count
		readCountsPerMillion[[i]] <- tmpMat$reads_per_million_miRNA_mapped
		names(readCountsPerMillion[[i]]) <- names(readCounts[[i]]) <- tmpMat$isoform_coords
		if(i == length(isoformFiles)){
			cat("\n")
		}
	}
	
	# Now merge them into a matrix
	allNmsReadCounts  <- unique(names(unlist(readCounts)))
	matReadCounts <- matrix(NA, 
			nr=length(allNmsReadCounts), 
			nc=length(isoformFiles), 
			dimnames=list(allNmsReadCounts,basename(isoformFiles)))
	for(i in 1:length(isoformFiles)){
		matReadCounts[names(readCounts[[i]]),i] <- readCounts[[i]]
	}	
	allNmsReadCountsPerMillion  <- unique(names(unlist(readCountsPerMillion)))
	matReadCountsPerMillion <- matrix(NA, 
			nr=length(allNmsReadCountsPerMillion), 
			nc=length(isoformFiles), 
			dimnames=list(allNmsReadCountsPerMillion,basename(isoformFiles)))
	for(i in 1:length(isoformFiles)){
		matReadCountsPerMillion[names(readCountsPerMillion[[i]]),i] <- readCountsPerMillion[[i]]
	}	
	
	matReadCounts <- .reduceToAliquots(matReadCounts, metadata)
	matReadCountsPerMillion <- .reduceToAliquots(matReadCountsPerMillion, metadata)
	
	isoformEsetReadCounts <- new("ExpressionSet", exprs=matReadCounts, annotation=platform)
	isoformReadCountsPerMillion <- new("ExpressionSet", exprs=matReadCountsPerMillion, annotation=platform)
	
	return(list(isoformEsetReadCounts = isoformEsetReadCounts,
					isoformReadCountsPerMillion = isoformReadCountsPerMillion,
					mirnaEsetReadCounts = mirnaEsetReadCounts,
					mirnaEsetReadCountsPerMillion = mirnaEsetReadCountsPerMillion))
}



.reduceToAliquots <- function(dat, metadata){
	id <- which(apply(metadata, 2, function(x){ sum(!is.na(match(colnames(dat),x)))}) > 0)
	tmp <- metadata[match(colnames(dat),metadata[,id]),]
	tmp <- tmp[,which(apply(tmp, 2, function(x){ sum(!is.na(x))}) > 0)]
	barcodeToColumn <- split(1:nrow(tmp), tmp$bcr_aliquot_barcode)
	newDat <- matrix(0, nr=nrow(dat), ncol=length(barcodeToColumn))
	colnames(newDat) <- names(barcodeToColumn)
	rownames(newDat) <- rownames(dat)
	for(j in 1:length(barcodeToColumn)){
		if(length(barcodeToColumn[[j]]) > 1){
			newDat[,j] <- rowMeans(dat[,barcodeToColumn[[j]]])
		}else{
			newDat[,j] <- dat[,barcodeToColumn[[j]]]
		}
	}
	newDat
}
