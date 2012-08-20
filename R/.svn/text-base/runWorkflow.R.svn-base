#
#  call this on a single directory which should contain raw data
#  for a single layer.
#

runWorkflow <- function(rawDataDir, workflow=NULL, logFile=TRUE, ...) {  
	# Set up log file
	if(logFile == "TRUE") {
		logFile <- setUpLogFile(logFile)
		cat("log file for this workflow: ", logFile, "\n")
	}

	# Handle user input
	userInput <- handleInput(rawDataDir, workflow, logFile)
	rawDataDir <- userInput$rawDataDir
	workflow <- userInput$workflow
	
	# Call the workflow
	res <- callWorkflows(rawDataDir=rawDataDir, 
			workflow=workflow, 
			logFile=logFile, ...)
	return(res)
}




