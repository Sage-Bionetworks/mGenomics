\name{agilentTxtWorkflow}
\alias{agilentTxtWorkflow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Implement workflow on Agilent Gene Expression or Array CGH Data
}
\description{

This function allows the user to implement supervised or unsupervised workflows for
Agilent data.

}
\usage{
agilentTxtWorkflow(rawDataDir, bio.var=NULL, adj.var=NULL, int.var=NULL, 
	rm.adj=TRUE, verbose=TRUE, logFile=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rawDataDir}{
Character string describing the directory in which we can find tab-delimited text files
containing raw data from an Agilent microarray experiment. Other files (including other 
tab-delimiated text files) are allowed in the directory and recursive search is used
to search any directories within directories.
}
  \item{bio.var}{
Model matrix parameterizing the biological variables of interest.  Note this option 
can only be set if it applies to all platforms that will be processed by the workflow.  
If multiple platforms exist we recommend first splitting the data into directories based
on their platform and then running the workflow on each individually.     
}
  \item{adj.var}{
Model matrix parameterizing the adjustment variables.  Note this option can only be set if it 
applies to all platforms that will be processed.  See note above on how to proceed if
data from multiple platforms exists.
}
  \item{int.var}{
A data frame with number of rows equal to the number of samples, with the unique levels of 
intensity-dependent effects. Each column parameterizes a unique source of intensity-dependent 
effect (e.g., array effects for column 1 and dye effects for column 2).   Note this option 
can only be set if it applies to all platforms that will be processed. See note above on how 
to proceed if data from multiple platforms exists.
}
  \item{rm.adj}{
A logical scalar.  If set to 'TRUE' the influence of the adjustment variables will be removed 
from the normalized data. 
 }
  \item{verbose}{
A logical scalar.  If set to 'TRUE' various messages will be passed to the R console providing 
updates for each workflow.
}
  \item{logFile}{
Used to create (if set equal to TRUE) or define a logFile (if set to name of an existing file).
  }
}
\details{

This function implements supervised and unsupervised workflows for Agilent data.  The normalization is 
based on the supervised normalization of microarrays (snm) framework.  The unsupervised workflows 
use a model containing a probe intercept term and intensity-dependent effects for dye and array to normalize the data.  
Supervised workflows can be built by specifying the biological, adjustment, or intensity-dependent variables, but
only if the data in the specified directory comes from a single platform. 

Note that this workflow defines the data for the Cy3 (green) or Cy5 (red) channels as  
the gProcessedSignal and rProcessedSignal values.

More information is available on the wiki for the metaGenomics project:
http://sagebionetworks.jira.com/metaGenomics/WIKI/home

}
\value{

A list of lists, one for each platform processed by the workflow.  The objects of the list that is returned
are named according to the platform type.  A list for a given platform contains the normalized data as a matrix 
in an eset as well as an empty list.  In the future this list will be poulated with statistics useful for donwstream
analysis of the data.  Note that for data obtained from experiments that use both green and red channels we return the data matrices 
so that the green channels are in columns 1 to n, where n is the # of distinct microarrays measuring that 
platform, and the red channels are in columns n+1 to 2n.  In general if the green channel is in column j for the
jth sample, then its corresponding red channel is in column n+j.

Note that we have encountered many strange features of Agilent .TXT files that make them difficult to process.
Given this, we developed this software to be robust to improperly formatted input data.  When we encounter such 
a file its corresponding columns in the data matrices are set equal to NA for every probe in every channel.  
These columns are returned to the user, but do not influence the normalization procedure.

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Brig Mecham <brig.mecham@sagebase.org>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (rawDataDir, bio.var, adj.var, int.var, rm.adj, verbose) 
{
    all.txt.files <- list.files(path = rawDataDir)[grep(".txt", 
        list.files())]
    if (length(all.txt.files) == 0) {
        cat("No Agilent Platforms Detected\n")
        fits <- "No Agilent Platforms Detected\n"
        class(fits) <- "try-error"
        return(fits)
    }
    platforms <- sapply(all.txt.files, getAgilentPlatform)
    if (platforms == "No Agilent Platforms Detected\n") {
        fits <- platforms
    }
    else {
        agilent.files <- as.numeric(which(!is.na(platforms)))
        all.agilent.files <- all.txt.files[agilent.files]
        platforms <- platforms[agilent.files]
        platforms <- standardizePlatforms(platforms)
        platform2file <- split(all.agilent.files, platforms)
        fits <- list()
        for (i in 1:length(platform2file)) {
            cat("Platform:", names(platform2file)[i], "\n")
            files <- platform2file[[i]]
            cat("Number of Samples:", length(files), "\n")
            obj <- loadAgilentFiles(files)
            if (length(files) > 1) {
                cat("Normalizing Data\n")
                platformObj <- normalizeAgilent(obj)
            }
            else {
                if (sum(!is.na(obj$red)) > 0) {
                  is.red = TRUE
                }
                if (sum(!is.na(obj$green)) > 0) {
                  is.green = TRUE
                }
                if (is.red & is.green) {
                  data <- cbind(obj$green, obj$red)
                }
                if (is.red & !is.green) {
                  data <- cbind(obj$red)
                }
                if (!is.red & is.green) {
                  data <- cbind(obj$green)
                }
                platformObj <- list(data = data, probePValues = obj$probePValues)
            }
            if (length(platform2file) > 1) {
                cat("---------------------------------------------------\n")
            }
            fits[[i]] <- platformObj
            nms <- names(platform2file)[i]
            nms <- gsub("-", ".", nms)
            names(fits)[i] <- nms
        }
    }
    fits
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
