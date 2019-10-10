\name{segmentPICS} 
\alias{segmentPICS}
\alias{segReadsGeneric}


\title{Segment the genome into candidate regions} 
\description{ 
    Pre-process bidirectional aligned reads data from a single ChIP-Seq experiment to detect candidate regions with a minimum number of forward and reverse reads. These candidate regions will then be processed by PICS.}
\usage{
	segmentPICS(data, dataC=NULL, map=NULL, minReads=2, minReadsInRegion=3,
    jitter=FALSE, dataType="TF", maxLregion=0, minLregion=100)
}


\arguments{
  \item{data}{A \code{linkS4class{GRanges}} object containing the IP reads. See details for more information on how to set up the data.}
  \item{dataC}{A \code{linkS4class{GRanges}} object containing the control reads. Set to NULL by default, i.e. no control.}
  \item{map}{A `RangedData' object containing the mappability profiles. Set to NULL by default, i.e. no profiles.}
  \item{minReads}{The minimum number of F/R reads to be present in the sliding window.}
  \item{minReadsInRegion}{The minimum number of F/R reads to be present in the region.}
  \item{jitter}{A logical value stating whether some noise should be added to the read locations. This is recommended if the read positions have lots of duplicates.}
  \item{dataType}{Type of experiment. "TF" or "H".}
  \item{maxLregion}{The maximum length.}
  \item{minLregion}{The minimum length.}
}

\value{
An object of class \code{segReadsList} containing the results for all regions pre-processed. 
}
\references{
X. Zhang, G. Robertson, M. Krzywinski, K. Ning, A. Droit, S. Jones, and R. Gottardo, ``PICS: Probabilistic Inference for ChIP-seq'' arXiv, 0903.3206, 2009.
}
\author{ Xuekui Zhang, Arnaud Droit <\email{arnaud.droit@crchuq.ualaval.ca}> and Raphael Gottardo <\email{rgottard@fhcrc.org}>}


\examples{
# Read data
path<-system.file("extdata",package="PICS")
## Note that the col name for the chromosome needs to be space and not chr
dataIP<-read.table(file.path(path,"Treatment_tags_chr21_sort.bed"),header=TRUE,colClasses=c("factor","integer","integer","factor"))
dataIP<-as(dataIP,"GRanges")

dataCont<-read.table(file.path(path,"Input_tags_chr21_sort.bed"),header=TRUE,colClasses=c("factor","integer","integer","factor"))
dataCont<-as(dataCont,"GRanges")

map<-read.table(file.path(path,"mapProfileShort"),header=TRUE,colClasses=c("factor","integer","integer","NULL"))
map<-as(map,"GRanges")
seg<-segmentPICS(dataIP, dataC=dataCont, map=map, minReads=1)
}

\seealso{
	\code{\linkS4class{segReadsList}}
}

\keyword{functions}
