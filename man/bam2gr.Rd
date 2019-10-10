\name{bam2gr}
\alias{bam2gr}

\title{pre-process bam files}

\description{
Reads a bam file using \code{Rsamtools} and extract the reads for each chromosome.
}

\usage{
bam2gr(bamFile, chr=NULL, PE=FALSE, verbose=FALSE)
}

\arguments{
\item{bamFile}{A character string, the name of the .bam file to read.}
\item{chr}{An optional character string. If specified, only the selected chromosome will be returned. Speed up the computation.}
\item{PE}{A \code{logical}. This should be set to TRUE for paired-end sequencing data.}
\item{verbose}{A \code{logical}. Print additional information about the data.}
}

\value{
Returns a \code{GRanges} of all the reads for each chromosome. 
}

\note{
The user might encounter a memory allocation error when using bam files of bigger sizes. Splitting the file by chromosome before calling \code{bam2gr} will solve this issue.

For Paired-End data, non matched reads are discarded.
}

\author{Renan Sauteraud}

\seealso{
\code{\link{segmentPICS}}
}
\keyword{functions}
