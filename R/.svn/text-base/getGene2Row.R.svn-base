getGene2Row <- function(gene2row, abatchAnnotation){ 
	if(is.null(gene2row)) {
		other.env <- new.env()
		data(list=c(abatchAnnotation),envir=other.env)
		gene2row.cdf <- get(ls(other.env)[1],1,other.env)
	}else{
		gene2row.cdf <- try(get(abatchAnnotation,gene2row))
		if(class(try(gene2row.cdf))=="try-error"){
			stop("user provided gene2row list does not contain mapping for platform ", abatchAnnotation)
			gene2row.cdf <- NULL
		}
	}
	gene2row.cdf
}
