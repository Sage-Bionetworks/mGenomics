
getRnaSeq <- function(rawDataDir, platform, metadata){
	if(is.null(platform)){
		var <- "Must provide valid platform.\n"
		class(var) <- 'try-error'
		return(var)
	}
	
	tmp <- list.files(rawDataDir, full.names=T, recursive=T)
	geneFiles <- tmp[grep("gene.quantification.txt", basename(tmp), fixed=T)]
	
	if(length(geneFiles) >0 ){
		rawCounts <- medianNorm <- rpkm <- vector("list", length(geneFiles))
  	for(i in 1:length(geneFiles)){
#		for(i in 1:5){
			cat("\rloading file ", i, " ", basename(geneFiles[i]))
			tmpMat <- read.delim(file=geneFiles[i], as.is=T, header=T, row.names=NULL, quote="")
			rawCounts[[i]] <- tmpMat$raw_counts
			medianNorm[[i]] <- tmpMat$median_length_normalized
			rpkm[[i]] <- tmpMat$RPKM
			symbol <- sapply(tmpMat$gene, function(x) { strsplit(x,"\\|")[[1]][1] })
			entrezId <- sapply(tmpMat$gene, function(x) { strsplit(x,"\\|")[[1]][2] })
			names(rpkm[[i]]) <- names(medianNorm[[i]]) <- names(rawCounts[[i]]) <- paste(as.character(entrezId), '_eg',sep="")
			if(i == length(geneFiles)){
				cat("\n")
			}
		}
	}

	# Now merge them into a matrix
	allNmsReadCounts  <- unique(names(unlist(rawCounts)))
	matReadCounts <- matrix(NA, 
			nr=length(allNmsReadCounts), 
			nc=length(geneFiles), 
			dimnames=list(allNmsReadCounts,basename(geneFiles)))
	for(i in 1:length(geneFiles)){
		matReadCounts[names(rawCounts[[i]]),i] <- rawCounts[[i]]
	}	

	allNmsReadCountsPerMillion  <- unique(names(unlist(rpkm)))
	matReadCountsPerMillion <- matrix(NA, 
			nr=length(allNmsReadCountsPerMillion), 
			nc=length(geneFiles), 
			dimnames=list(allNmsReadCountsPerMillion,basename(geneFiles)))
	for(i in 1:length(geneFiles)){
		matReadCountsPerMillion[names(rpkm[[i]]),i] <- rpkm[[i]]
	}	

	allNmsMedianReadCounts <- unique(names(unlist(medianNorm)))
	matMedianReadCounts <- matrix(NA, 
			nr=length(allNmsMedianReadCounts), 
			nc=length(geneFiles), 
			dimnames=list(allNmsMedianReadCounts,basename(geneFiles)))
	for(i in 1:length(geneFiles)){
		matMedianReadCounts[names(medianNorm[[i]]),i] <- medianNorm[[i]]
	}
	
	matReadCounts <- .reduceToAliquots(matReadCounts, metadata)
	matReadCountsPerMillion <- .reduceToAliquots(matReadCountsPerMillion, metadata)
	matMedianReadCounts <- .reduceToAliquots(matMedianReadCounts, metadata)
	
	mrnaEsetReadCounts <- new("ExpressionSet", exprs=matReadCounts, annotation=platform)
	mrnaEsetReadCountsPerMillion <- new("ExpressionSet", exprs=matReadCountsPerMillion, annotation=platform)
	mrnaEsetMedianReadCounts <- new("ExpressionSet", exprs=matMedianReadCounts, annotation=platform)

	#### Exon counts
	exonFiles <- tmp[grep("exon.quantification.txt", basename(tmp), fixed=T)]

	rawCounts <- medianNorm <- rpkm <- vector("list", length(exonFiles))
	for(i in 1:length(exonFiles)){
		cat("\rloading file ", i, " ", basename(exonFiles[i]))
		tmpMat <- read.delim(file=exonFiles[i], as.is=T, header=T, row.names=NULL, quote="")
		rawCounts[[i]] <- tmpMat$raw_counts
		medianNorm[[i]] <- tmpMat$median_length_normalized
		names(medianNorm[[i]]) <- names(rawCounts[[i]]) <- as.character(tmpMat$exon)
		if(i == length(exonFiles)){
			cat("\n")
		}
	}

	# Now merge them into a matrix
	allNmsReadCounts  <- unique(names(unlist(rawCounts)))
	matReadCounts <- matrix(NA, 
			nr=length(allNmsReadCounts), 
			nc=length(geneFiles), 
			dimnames=list(allNmsReadCounts,basename(geneFiles)))

	allNmsReadCounts  <- unique(names(unlist(medianNorm)))
	matMedianReadCounts <- matrix(NA, 
			nr=length(allNmsReadCounts), 
			nc=length(geneFiles), 
			dimnames=list(allNmsReadCounts,basename(geneFiles)))
	for(i in 1:length(geneFiles)){
		matReadCounts[names(rawCounts[[i]]),i] <- rawCounts[[i]]
		matMedianReadCounts[names(medianNorm[[i]]),i] <- medianNorm[[i]]
	}	

	matReadCounts <- .reduceToAliquots(matReadCounts, metadata)
	matMedianReadCounts <- .reduceToAliquots(matMedianReadCounts, metadata)
	exonEsetReadCounts <- new("ExpressionSet", exprs=matReadCounts)
	exonEsetMedianReadCounts <- new("ExpressionSet", exprs=matMedianReadCounts)

	#### Junction Files	
	junctionFiles <- tmp[grep("spljxn.quantification.txt", basename(tmp), fixed=T)]
	rawCounts <- medianNorm <- rpkm <- vector("list", length(junctionFiles))
	for(i in 1:length(junctionFiles)){
		cat("\rloading file ", i, " ", basename(junctionFiles[i]))
		tmpMat <- read.delim(file=junctionFiles[i], as.is=T, header=T, row.names=NULL, quote="")
		rawCounts[[i]] <- tmpMat$raw_counts
		names(rawCounts[[i]]) <- as.character(tmpMat$junction)
		if(i == length(junctionFiles)){
			cat("\n")
		}
	}

# Now merge them into a matrix
	allNmsReadCounts  <- unique(names(unlist(rawCounts)))
	matReadCounts <- matrix(NA, 
			nr=length(allNmsReadCounts), 
			nc=length(geneFiles), 
			dimnames=list(allNmsReadCounts,basename(geneFiles)))
	for(i in 1:length(geneFiles)){
		matReadCounts[names(rawCounts[[i]]),i] <- rawCounts[[i]]
	}	
	matReadCounts <- .reduceToAliquots(matReadCounts, metadata)
	
	junctionEset <- new("ExpressionSet", exprs=matReadCounts)
	return(list(junctionEset = junctionEset,
					exonEsetMedianReadCounts = exonEsetMedianReadCounts, 
					exonEsetReadCounts = exonEsetReadCounts,
					mrnaEsetReadCounts = mrnaEsetReadCounts, 
					mrnaEsetReadCountsPerMillion = mrnaEsetReadCountsPerMillion,
					mrnaEsetMedianReadCounts  = mrnaEsetMedianReadCounts ))
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
