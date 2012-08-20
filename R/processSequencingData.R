#platforms <- c('illuminaga_rnaseq', 'illuminaga_mirnaseq', 'illuminahiseq_rnaseq', 'illuminahiseq_mirnaseq')
#tcgaQry <- synapseQuery('select id, name from study where study.repository == "TCGA" and study.parentId == "syn150935"')
#mgTcgaQry <- synapseQuery('select id, name, acronym from study where study.repository == "TCGA" and study.parentId == "syn275039"')
#
## loop across all platforms
#for(i in 1:length(platforms)){
#	entitiesToProcess <- sapply(tcgaQry$study.id, function(x){.getDataToProcess(x, platforms[i]) })
#	entitiesToProcess  <- entitiesToProcess [which(!sapply(entitiesToProcess, is.null))]
## loop across all studies
#	sapply(names(entitiesToProcess)[7:9], function(studyId){
#				cat("\n\n\n Starting study ", studyId)
#				entities <- processSequencingData(studyId, platforms[i])
#			}) -> res	
#}

processSequencingData <- function(studyId, platform){
	cat("Step 1: Querying Synapse to identify most recent batches.\n")
	myQ <- synapseQuery(paste('SELECT id, name, platform FROM entity WHERE entity.parentId == "', studyId, '" AND entity.platform == "', platform, '"', sep=""))
	myQ <- myQ[!grepl('mage',myQ$entity.name),]	
	plat2row <- split(1:nrow(myQ), as.character(unlist(myQ$entity.platform)))
	sapply(plat2row, function(p){
				nms2 <- myQ[p,]
				batchInfo <- sapply(nms2$entity.name, .getBatchInfo)
				serial2col <- split(1:ncol(batchInfo), batchInfo[1,])
				serial2rev <- split(as.numeric(batchInfo[2,]), batchInfo[1,])
				rev2col <- split(1:ncol(batchInfo), batchInfo[2,])		
				ids <- c(9)
				for(j in 1:length(serial2col)){
					ids[j] <- serial2col[[j]][which.max(serial2rev[[j]])]
				}
				p[ids]
			}) -> ids
	myQ <- myQ[unlist(ids),]
	
	if( nrow(myQ) == 0 ){
		return(NA)
	}
	
	cat("Step 2: Loading entities.\n")
	theseFiles <- lapply(as.list(myQ$entity.id), function(x){
				tmp <- downloadEntity(x)
				list.files(tmp$cacheDir, full.names=T, recursive=T)
			})
	
	theseFiles <- unlist(theseFiles)
	theseDirs <- unique(dirname(theseFiles))
	myDir <- tempfile()
	dir.create(myDir)
	for( i in theseDirs){
		system(paste("ln -s ", i, '/* ', myDir, sep=""), ignore.stdout=T, ignore.stderr=T)
	}
		
	cat("Step 3: Loading metadata\n")
	studyEntity <- getEntity(studyId)
	acronym <- annotValue(studyEntity, 'acronym')
	mgStudy <- mgTcgaQry$study.id[match(acronym, mgTcgaQry$study.acronym)]
	qry <- synapseQuery(paste('select id from entity where entity.parentId == "', mgStudy, '" and entity.name=="', acronym, '_mergedClinical"',sep=""))
	clinEntity <- loadEntity(qry$entity.id)
	metadata <- clinEntity$objects$metadata
	
	if(platform == 'illuminaga_mirnaseq' | platform == 'illuminahiseq_mirnaseq'){
		cat("Step 4: Buildng esets\n")
		esets <- getMirnaSeq(myDir, platform, metadata)
		cat("Step 5: Contributing to Synapse\n")
		newEntity_1 <- .createCoherentSequencingEntity(platform, studyId, 'isoformReadCounts', "isoform",esets[[1]],myQ)
		newEntity_2 <- .createCoherentSequencingEntity(platform, studyId, 'isoformReadCountsPerMillion', "isoform",esets[[2]],myQ)
		newEntity_3 <- .createCoherentSequencingEntity(platform, studyId, 'mirnaReadCounts', "miRNA",esets[[3]],myQ)
		newEntity_4 <- .createCoherentSequencingEntity(platform, studyId, 'mirnaReadCountsPerMillion', "miRNA",esets[[4]],myQ)
		return(list(newEntity_1, newEntity_2, newEntity_3, newEntity_4))		
	}else{
		cat("Step 4: Buildng esets\n")
		esets <- getRnaSeq(myDir, platform, metadata)
		cat("Step 5: Contributing to Synapse\n")
		newEntity_1 <- try(.createCoherentSequencingEntity(platform, studyId, 'junction', "junction",esets[[1]],myQ))
		newEntity_2 <- try(.createCoherentSequencingEntity(platform, studyId, 'medianReadCounts', "exon",esets[[2]],myQ))
		newEntity_3 <- try(.createCoherentSequencingEntity(platform, studyId, 'exonReadCounts', "exon",esets[[3]],myQ))
		newEntity_4 <- try(.createCoherentSequencingEntity(platform, studyId, 'mrnaReadCounts', "mRNA",esets[[4]],myQ))
		newEntity_5 <- try(.createCoherentSequencingEntity(platform, studyId, 'mrnaReadCountsPerMillion', "mRNA",esets[[5]],myQ))
		newEntity_6 <- try(.createCoherentSequencingEntity(platform, studyId, 'mrnaMedianReadCounts', "mRNA",esets[[6]],myQ))
	}
}

.createCoherentSequencingEntity <- function(platform, studyId, nameSuffix=NULL, summaryType=NULL, eset, batch){
	molecFeatureType <- .loadPlatformMolecularFeatures(platform)
	dataTypes <- .loadPlatformDataTypes(platform)
	if(!is.null(nameSuffix)){ 
		if(platform == 'pd.genomewidesnp.6' & grepl('geneRatios', nameSuffix)){
			summaryType <- 'gene'
		}
	}
	thisName <- paste("Normalized ", dataTypes," - ", platform, sep="") 
	if(!is.null(nameSuffix)){
		thisName <- paste(thisName, nameSuffix, sep=" - ")
	}
	## SET UP THE LAYER - AND CLEAR OUT OLD OBJECTS AND FILES IF ALREADY EXIST
	# metaGenomics project Id
	tmpDS <- getEntity(studyId)
	mgId <- synapseQuery(paste('select id from study where study.name=="',paste(propertyValue(tmpDS,"name"),' - SNM Normalized',sep=""),'" and study.parentId=="syn275039"',sep=""))$study.id
	thisLay <- synapseQuery(paste("SELECT name, id FROM entity WHERE entity.parentId == '", mgId, "'", sep=""))
	if( any(thisLay$entity.name == thisName, na.rm=T) ){
		myLayer <- loadEntity(thisLay$entity.id[thisLay$entity.name == thisName])
		myLayer <- deleteObject(myLayer, names(myLayer$objects))
		myLayer <- deleteFile(myLayer, myLayer$files)
		myLayer <- updateEntity(myLayer)
	} else{
		myLayer <- createEntity(Data(list(name=thisName,
								parentId=mgId)))
	}
	cat(propertyValue(myLayer,'id'),"\n")
	annotValue(myLayer,'dataType') <- dataTypes
	annotValue(myLayer,'molecFeatureType') <- molecFeatureType
	myLayer <- addObject(myLayer, eset)
	
	## GRAB THE DATASET FOR TRANSFER OF ANNOTATIONS AND PROPERTIES
	desc <- paste("This entity consists of the level 3 data from TCGA.  The data were merged into a single matrix.",
			"  Included in this entity is the R binary object eset of class 'ExpressionSet'.", sep="")
	propertyValue(myLayer, "description") <- desc
	propertyValue(myLayer, "disease") <- propertyValue(tmpDS, "disease")
	propertyValue(myLayer, "species") <- propertyValue(tmpDS, "species")
	propertyValue(myLayer, "tissueType") <- propertyValue(tmpDS, "tissueType")	
	propertyValue(myLayer, "platform") <- platform
	propertyValue(myLayer, "numSamples") <- ncol(exprs(eset))
	annotValue(myLayer, "status") <- "processed"
	annotValue(myLayer, "processingAlgorithm") <- summaryType
	annotValue(myLayer, "repository") <- annotValue(tmpDS, "repository")
	annotValue(myLayer, "acronym") <- annotValue(tmpDS, "acronym")
	annotValue(myLayer, "numBatches") <- nrow(batch)	
	myLayer <- updateEntity(myLayer)
	myLayer <- storeEntity(myLayer)
	myLayer
}


.loadPlatformDataTypes <- function(x){ 	
	types <- list('hthgu133a' = "Gene Expression", #
			'agilentg4502a_07_2' = "Gene Expression",
			'hg-cgh-244a' = "Copy Number",
			'hg-cgh-415k_g4124a' = "Copy Number", #
			'h-mirna_8x15k' = "miRNA Expression",
			'agilentg4502a' = "Gene Expression", #
			'agilentg4502a_07_1' = "Gene Expression", #
			'agilentg4502a_07_3' = "Gene Expression", #
			'pd.genomewidesnp.6' = "Copy Number", #
			'hgu133plus2' = "Gene Expression", #
			'huex10stv2' = "Gene Expression",			
			'cgh-1x1m_g4447a' = "Copy Number", #
			'humanmethylation450' = 'Methylation',
			'h-mirna_8x15kv2' = "miRNA Expression",
			'illuminaga_rnaseq' = "Gene Expression",
			'illuminaga_mirnaseq' = 'miRNA Expression',
			'illuminahiseq_rnaseq' = 'Gene Expression',
			"mda_rppa_core" = "Protein Expression",
			'illuminahiseq_mirnaseq' = 'miRNA Expression',
			'illuminaga_dnaseq' = 'somatic Mutations') 
	if(x %in% names(types)){
		return(as.character(unlist(types[x])))
	}else{
		stop("Cannot find platform", x,"\n\n\n")
	}	
}

.loadPlatformSummaryTypes <- function(x){ 	
	types <- list('hthgu133a' = "gene", #
			'agilentg4502a_07_2' = "probe",
			'hg-cgh-244a' = "probe",
			'hg-cgh-415k_g4124a' = "probe", #
			'h-mirna_8x15k' = "probe",
			'agilentg4502a' = "probe", #
			'agilentg4502a_07_1' = "probe", #
			'agilentg4502a_07_3' = "probe", #
			'pd.genomewidesnp.6' = "probe", #
			'hgu133plus2' = "gene", #
			'huex10stv2' = "gene",
			'cgh-1x1m_g4447a' = "probe", #
			'h-mirna_8x15kv2' = "probe",			
			"mda_rppa_core" = "anitbody",
			'humanmethylation450' = 'probe',			
	) 
	if(x %in% names(types)){
		return(as.character(unlist(types[x])))
	}else{
		stop("Cannot find platform", x,"\n\n\n")
	}	
}

.loadPlatformMolecularFeatures <- function(x){ 	
	types <- list('hthgu133a' = "RNA", #
			'agilentg4502a_07_2' = "RNA",
			'hg-cgh-244a' = "DNA",
			'hg-cgh-415k_g4124a' = "DNA", #
			'h-mirna_8x15k' = "RNA",
			'agilentg4502a' = "RNA", #
			'agilentg4502a_07_1' = "RNA", #
			'agilentg4502a_07_3' = "RNA", #
			'pd.genomewidesnp.6' = "DNA", #
			'hgu133plus2' = "RNA", #
			'huex10stv2' = "RNA",
			'cgh-1x1m_g4447a' = "DNA", #
			'h-mirna_8x15kv2' = "DNA",
			'humanmethylation450' = 'DNA',
			'illuminaga_rnaseq' = "RNA",
			'illuminaga_mirnaseq' = 'RNA',
			'illuminahiseq_rnaseq' = 'RNA',
			'illuminahiseq_mirnaseq' = 'RNA',
			"mda_rppa_core" = "Protein",
			'illuminaga_dnaseq' = 'DNA') 
	if(x %in% names(types)){
		return(as.character(unlist(types[x])))
	}else{
		stop("Cannot find platform", x,"\n\n\n")
	}	
}

.getBatchInfo <- function(name) {
	if(!grepl('Level',name)){
		return(name)
	}
	domain <- strsplit(name,"_")[[1]][1]
	m <- regexpr("Level_\\d\\.\\d+\\.\\d+", name,perl=TRUE)
	if(m[1] == -1){
		return(c(NA,NA))
	}
	mtch <- regmatches(name, m)
	mtch <-strsplit(mtch, '\\.')[[1]]
	serialIndex <- mtch[2]
	revisionNumber <- mtch[3]
	c(serialIndex, revisionNumber)
}

.getDataToProcess <- function(studyId, platform){
		myQ <- synapseQuery(paste('SELECT id, name, platform FROM entity WHERE entity.parentId == "', studyId, '" AND platform == "', platform, '"', sep=""))
		myQ <- myQ[!grepl('mage',myQ$entity.name),]	
	} -> mis
