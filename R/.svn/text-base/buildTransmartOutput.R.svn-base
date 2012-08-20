buildTransmartOutput <- function(inf.obj, varName, varNames,fname){	
	transmart.matrix <- data.frame(rownames(inf.obj$cfs),
			round(inf.obj$cfs[,1],2),
			round(inf.obj$cfs[,2],2),
			round(inf.obj$pval,3),
  		round(inf.obj$qval,3))
	if(sum(!is.na(inf.obj$pval)) > 0) { 
		transmart.matrix[which(transmart.matrix[,4] < 0.001),4] <- 0.001
		transmart.matrix[which(transmart.matrix[,5] < 0.001),5] <- 0.001
	}
	rownames(transmart.matrix) <- NULL
	colnames(transmart.matrix) <- c("Probe Set ID",  
			paste(varName, "=> Expression in", varNames[1], ".Estimate"),
			paste(varName, "=> Effect of", varNames[2],".FoldChange"),
			paste(varName, "=> P value",".RawPValue"),
			paste(varName, "=> Q value",".RawPValue"))
	write.table(transmart.matrix, 
			quote=FALSE,
			row.names=FALSE,
			file=fname)
}
