\name{affyWorkflow.SNM}
\alias{affyWorkflow.SNM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Implement workflow on Affymetrix gene expression data.
}
\description{
This function allows the user to implement supervised or unsupervised workflows for
affymetrix data.  
}
\usage{
affyWorkflow.SNM(rawDataDir, gene2row = NULL, bio.var = NULL, adj.var = NULL, 
int.var = NULL, rm.adj = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{rawDataDir}{
Character string describing directory that contains raw data files.  These directories should contain raw data
from either Affymetrix (e.g. .CEL files) or Agilent (e.g. .txt files).  They can also contain files of other types.
}
  \item{gene2row}{
List of lists specifying the probe set to probe mapping for each platform that will be processed by the 'affy' workflows.  This is intended
to allow users to easily process data using custom designed probe sets. 
}
  \item{bio.var}{
Model matrix parameterizing the biological variables of interest.  Note this option can only be set if it applies to all platforms that will be processed.
}
  \item{adj.var}{
Model matrix parameterizing the adjustment variables.  Note this option can only be set if it applies to all platforms that will be processed.}
  \item{int.var}{
A data frame with number of rows equal to the number of samples, with the unique levels of intensity-dependent effects. Each column parameterizes a unique source of intensity-dependent effect (e.g., array effects for column 1 and dye effects for column 2).
}
  \item{rm.adj}{
A logical scalar.  If set to 'TRUE' the influence of the adjustment variables will be removed from the normalized data. 
 }
  \item{verbose}{
A logical scalar.  If set to 'TRUE' various messages will be passed to the R console providing updates for each workflow.
}
}
\details{

This function implements supervised and unsupervised workflows for Affymetrix data.  The normalization is 
based on the supervised normalization of microarrays (snm) framework.  The unsupervised workflows 
use a model containing a probe intercept term and intensity-dependent array effects to normalize the data.  
Supervised workflows can be built by specifying the biological, adjustment, or intensity-dependent variables.

At this point the workflows will successfully process data from 35 of the most commonly used Affymetrix
gene expression platforms.  Other platforms can be supported upon request.

Note that for each platform we built probe sets by aligning probe sequences to the RefSeq database
available from NCBI. Each probe set measures one and only one gene, where a gene is defined as the set of transcripts
corresponding to a given ENTREZ Gene ID.  No two genes are measured by any two probe sets. Our motivation 
for this was to simplify downstream analyses, (including, but not limited to merging data across 
platforms/technologies), and to leverage the potential of these  data to identify splice-variant specific 
expression.  

To distinguish our probe sets from the traditional affymetrix probe sets we developed the following
naming convention: ENTREZ.GENE.ID_mg (for metaGenomics transcript).  Including the Entrez GENE ID in 
the probe set name is useful when trying to learn more about its corresponding gene and simplifies 
merging data across platforms/technologies.

More information is available on the wiki for the metaGenomics project:
http://sagebionetworks.jira.com/metaGenomics/WIKI/home

}
\value{

A list of lists, one for each platform processed by the workflow.  
The objects of the lists are named according to the platform type (e.g. 'hgu133a').  Each sublist contains
the normalized and summarized data as an eset and a set of statistics useful for downstream analysis.

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
function (rawDataDir, gene2row = NULL, bio.var = NULL, adj.var = NULL, 
    int.var = NULL, rm.adj = FALSE, verbose = FALSE) 
{
    env <- new.env()
    platforms <- loadAvailablePlatforms()
    cel.files <- list.celfiles(rawDataDir, full.names = T)
    if (length(cel.files) == 0) {
        msg <- "No CEL files detected\n"
        class(msg) <- "try-error"
        return(msg)
    }
    cdfs <- sapply(cel.files, whatcdf)
    eme <- sapply(unique(cdfs), function(x) {
        if (x \%in\% platforms & sum(cdfs == x) > 1) {
            if (verbose) {
                cat("Platform:", x, "\n")
            }
            abatch <- ReadAffy(filenames = names(which(cdfs == 
                x)))
            scan.date <- abatch@protocolData$ScanDate
            int <- exprs(abatch)
            if (verbose) {
                cat("Number of Samples:", ncol(int), "\n")
            }
            gene2row.cdf <- getGene2Row(gene2row, annotation(abatch))
            pms <- int[unlist(gene2row.cdf), ]
            pms[pms <= 1] <- 1
            data <- log2(pms)
            if (is.null(int.var)) {
                int.var <- data.frame(array = factor(1:ncol(data)))
            }
            if (verbose) {
                cat("Normalizing Data\n")
            }
            snm.fit <- snm(data, bio.var, adj.var, int.var, diagnose = FALSE, 
                rm.adj = rm.adj, verbose = FALSE)
            colnames(snm.fit$norm.dat) <- sampleNames(abatch)
            gene2row.tmp <- split(1:nrow(pms), rep(names(gene2row.cdf), 
                sapply(gene2row.cdf, length)))
            if (verbose) {
                cat("Summarizing Data\n")
            }
            fits <- fit.pset(snm.fit$norm.dat, gene2row.tmp)
            if (verbose) {
                cat("Building metaGEO Object\n")
            }
            retval <- buildMetaGeoObject(data = fits$estimated.rna.concentration, 
                singular.values = fits$singular.values, probe.weights = fits$probe.weights)
            retval@protocolData$ScanDate <- abatch@protocolData$ScanDate
            string <- annotation(abatch)
            assign(string, retval, envir = env)
            if (length(cdfs) > 1) {
                cat("---------------------------------------------------\n")
            }
        }
        else if (sum(cdfs == x) < 2) {
            retval <- paste("Too few cel files for", x, "platform.")
            cat("Too few cel files for platform", x, ".\n")
            if (length(cdfs) > 1) {
                cat("---------------------------------------------------\n")
            }
            assign(x, retval, envir = env)
        }
        else {
            retval <- paste("Platform ", x, " not yet supported.")
            cat("Platform ", x, " not yet supported.\n")
            if (length(cdfs) > 1) {
                cat("---------------------------------------------------\n")
            }
            assign(x, retval, envir = env)
        }
    })
    as.list(env)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
