perform.all.inference <- function(dat, adj.var, meta.data, i=1, sva=FALSE) { 
# Move through each unique covariate
	results <- list();
	for( cov.num in 1:ncol(meta.data)) {
		covariate <- meta.data[[cov.num]]
#		cat("variable",cov.num," class",class(covariate),"\n")
		which.nas <- is.na(covariate)
		if(class(covariate) == "factor" | class(covariate) == "character") {
# Its a factor, remove unused levels			
			if(class(covariate) == "factor") {
				covariate <- as.character(unlist(covariate))
			}
			
			if(class(covariate) == "character") {
				covariate <- factor(covariate)
			}
			
# Handle missing values.  This adjusts the covariate and data matrix to remove samples with missing data
			if(sum(is.na(covariate)) > 0) { 
				valid.samples <- which(!is.na(covariate))
				covariate <- factor(covariate[valid.samples])
				dat.c <- dat[,valid.samples]
			}else{
				dat.c <- dat
			}
			
# Move to the next covariate if there is only one level
			if(length(unique(covariate)) == 1){
				results[[cov.num]] <- NA
				next
			}
			
# If one unique level for each sample then stop			
			if(ncol(dat.c) == nlevels(covariate)){
				results[[cov.num]] <- NA
				next;
			}
			
			X <- model.matrix(~ covariate)
			res <- dat.c - t(X %*% solve(t(X) %*% X) %*% t(X) %*% t(dat.c))
			lf <- lm(dat.c[i,] ~ X)
			results[[cov.num]] <- list(
					var.num = cov.num,
					lf.fit = lf,
					no.samples = ncol(dat.c),
					which.nas = which.nas,
					n.levels = unique(covariate))
		}
	}
	return(results)
}

