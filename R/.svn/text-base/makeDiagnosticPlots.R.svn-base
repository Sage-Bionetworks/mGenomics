makeDiagnosticPlots <- function(dat, root.dir=NULL, mat.rank=1, sampleRows=10000) {
	# Down sample the rows.  Useful if the matrix has many rows.
	if(sampleRows > 0 & sampleRows < nrow(dat)){ 
		dat <- dat[sample(1:nrow(dat), sampleRows),]
	}
	# Make directory to store latent structure plots
	if(is.null(root.dir)){
		system("mkdir -p diagnostics/latent.structure/")
		root.dir <- '.'
	}else{
		system(paste("mkdir -p ",root.dir,"/diagnostics/latent.structure/",sep=""))
	}
	
	# Take decomposition
	non.nas <- which(rowSums(is.na(dat)) == 0)
	dat <- dat[non.nas,]
	u <- fs(dat)
	
	# Make xyplot with x-axis the samples, y-axis the eigengene
	sapply(1:mat.rank, function(x) {
				png(file=paste(root.dir,"/diagnostics/latent.structure/ea_",x,".png",sep=""))
				try(print(xyplot(u$v[,x] ~ 1:ncol(dat), 
										xlab="Samples",ylab=paste("Eigengene",x),  
										main=paste("Eigengene",x))))
				dev.off()
			}) -> tmp
	
	# Make xyplot comparing adjacent eigengenes.
	stopAt <- mat.rank
	if(mat.rank == 1){
		stopAt <- 2
	}
	sapply(1:stopAt, function(x) {
				png(file=paste(root.dir,"/diagnostics/latent.structure/ea_",x,"_ea_",x+1,".png",sep=""))
				print(xyplot(u$v[,x+1] ~ u$v[,1], xlab=paste("Eigengene",x),ylab=paste("Eigengene",x+1)))
				dev.off()
			}) -> tmp
	
	# Make xyplot comparing mean to eigenarray.
	m <- rowMeans(dat,na.rm=TRUE)
	sapply(1:mat.rank, function(x) { 				
				png(file=paste(root.dir,"/diagnostics/latent.structure/mean_eg_",x,".png",sep=""))
				print(xyplot(u$u[,x] ~ m[non.nas], xlab="Average Estimated \nRNA Concentration",ylab=paste("Eigenarray",x)))
				dev.off() 
			}) -> tmp
	
	# Make xyplot showing the eigenweights
	png(file=paste(root.dir,"/diagnostics/latent.structure/Singular_Values.png",sep=""))
	print(xyplot(u$d ~ 1:ncol(dat), 
					xlab="Basis Vectors", 
					ylab="Singular Values"))
	dev.off()
	list.files(paste(root.dir,"/diagnostics",sep=""),
			recursive=TRUE,
			full.names=TRUE)
}
