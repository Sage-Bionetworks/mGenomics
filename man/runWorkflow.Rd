\name{runWorkflow}
\alias{runWorkflow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Implement unsupervised workflows on microarray data
}
\description{
This function implements workflows on microarray data from Affymetrix and Agilent technologies.  Workflows for other technologies
are under development and will be available soon.  
}
\usage{
runWorkflow(rawDataDir, workflow = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rawDataDir}{
Character string describing directory that contains raw data files.  These directories should contain raw data
from either Affymetrix (e.g. .CEL files) or Agilent (e.g. .txt files).  They can also contain files of other types.
}
  \item{workflow}{
Character string specifying a specific workflow for Affymetrix data sets.  Currently accepted values are 'snm' for the snm framework,
'rma' for RMA, 'dchip' for dCHIP, 'mas5' for mas5, and 'gcrma' for GCRMA. Note that the 'snm' framework will be used to normalize
any data from the agilent technology. 
implement both 
}
  \item{...}{
Parameters passed to technology specific workflows. For an example see help file for affyGeWorkflow.SNM.R
  }
}

\details{

This is a wrapper function to implement either supervised or unsupervised workflows on data for Affymeterix and Agilent 
technologies.  An unsupervised workflow is one that implements a minimal study-specific model (mSSM).  The mSSMs are
technology specific and are designed to remove commonly encountered confounders such as intensity-dependent effects for Affymetrix 
single channel data and both intensity and dye effects from two-color Agilent data.  For more information on the 
technology specific workflows please consult the help files for the affyGEWorkflow and agilentGeWorkflow functions.  Note
that this function also allows user to call traditional algorithms such as RMA, dChip, Mas5, and GCRMA.  

A supervised workflow is any that implements a study-specific model (SSM). SSMs are designed by an analyst to 
remove confounders specific to a given data set (e.g. batch effects).  The process of designing and implementing the SSM 
is based on the Supervised Normalization of Microarrays (snm) framework, allowing for a diverse set of normalization 
schemes that can handle data from single and double channel microarray technologies.  Note that if the specified directory 
contains data from multiple platforms, then the user cannot use this function to implement supervised workflows on 
each platform independently.  In this situation, an unsupervised workflow will be applied to the data from each platform.  
If the user wishes to run supervised workflows on each platform we suggest they split the data into multiple 
directories then call the function on each.
 
More information is available on the wiki for the metaGenomics project:
http://sagebionetworks.jira.com/metaGenomics/WIKI/home
  
}
\value{
A list containing one object for each platform processed by the workflow.  The objects of the lists are
named according to the platform type (e.g. 'hgu133a').
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

## If the rawdata is in a directory called foo, then:
## to run snm 
## fits <- runWorkflow("foo","snm")

## to run rma 
## fits <- runWorkflow("foo","rma")

## Other valid parameters are "dchip","mas5","gcrma"

# For users interested in processing Agilent data simply run:
# fits <- runWorkflow("foo")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line