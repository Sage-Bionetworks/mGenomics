
setMethod(
		f = "buildmGenomicsObject",
		signature = signature("matrix", "character"),
		definition = function(mat, annotation){
			eset<- new('ExpressionSet', 
					exprs=mat)
			annotation(eset) <- annotation
			list(eset=eset,
					statistics=NULL)
		})

setMethod(
		f = "buildmGenomicsObject",
		signature = signature("matrix", "character","list"),
		definition = function(mat, annotation, statistics){
			eset<- new('ExpressionSet', 
					exprs=mat)
			annotation(eset) <- annotation
			list(eset=eset, 
					statistics=statistics)
		})


#buildmGenomicsObject <- function(mat, statistics=NULL, annotation=NULL, ... ) {	
#	if(class(data) == 'matrix'){
#		if(is.null(statistics)){
#			
#		}
#		list(eset=eset,
#				statistics=statistics)
	#}else{
#		stop("ERROR: data object must be a matrix or object of class snmPSET.\n")
#	}
#}

#setClass("mGenomics", 
#		representation(singular.values='matrix', 
#				probe.weights='numeric', 
#				top10studies='matrix',
#				inference='list',
#				keywords="character"), 
#		contains="ExpressionSet")

#setGeneric("mGenomics", function(x, ...) standardGeneric("mGenomics"))
#setMethod(mGenomics, "missing", function(x, ...) new("mGenomics", ...))
#setMethod(mGenomics, "ExpressionSet", function(x, ...) {
#			new("mGenomics", assayData=assayData(x), phenoData=phenoData(x),
#					featureData=featureData(x), experimentData=experimentData(x),
#					protocolData=protocolData(x), annotation=annotation(x), ...)
#		})
