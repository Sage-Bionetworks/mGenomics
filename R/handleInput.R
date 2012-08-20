setMethod(
		f = "handleInput",
		signature = signature("character", "character", "character"),
		definition = function(rawDataDir, workflow, logFile){
			workflow <- tolower(workflow)
			if( !(any (workflow %in% c("rma","dchip","mas5","gcrma","snm")) & !is.null(workflow))) {
				msg <- paste("User supplied workflow",
						workflow,
						"not valid.  Must be one of snm, rma, dchip, gcrma, mas5")
				write.table(msg,file=logFile,append=TRUE,
						col.names=FALSE,row.names=FALSE,
						quote=FALSE)
				stop(msg)
			}
			if(!file.exists(rawDataDir)){
				msg <- paste("Cannot find raw data directory", rawDataDir,"\n")
				write.table(msg,file=logFile,append=TRUE,
						col.names=FALSE,row.names=FALSE,
						quote=FALSE)
				stop(msg)
			}
			return(list(rawDataDir=rawDataDir, workflow=workflow))		
		})

