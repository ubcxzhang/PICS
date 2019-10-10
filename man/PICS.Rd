\name{pics}
\alias{pics}
\alias{PICS}
\alias{pics-class}
\alias{as.list,AlignedRead-method}
\alias{code}
\alias{code,pics-method}
\alias{code,picsError-method}
\alias{code,picsList-method}
\alias{coerce}
\alias{coerce,picsList,RangedData-method}
\alias{coerce,picsList,data.frame-method}
\alias{coerce,AlignedRead,GenomeData-method}
\alias{coerce,RangedData,GenomeData-method}
\alias{coerce,data.frame,GenomeData-method}
\alias{coerce,data.frame,GRanges-method}
\alias{coerce,data.frame,picsList-method}
\alias{coerce,data.frame,picsList-method}
\alias{wigDensity}
\alias{wigDensity,pics-method}
\alias{wigDensity,picsList-method}
\alias{wigDensity,picsError-method}
\alias{score}
\alias{scoreForward}
\alias{scoreForward,pics-method}
\alias{scoreForward,picsList-method}
\alias{scoreForward,picsError-method}
\alias{scoreForward,data.frame-method}
\alias{scoreReverse}
\alias{scoreReverse,pics-method}
\alias{scoreReverse,picsList-method}
\alias{scoreReverse,picsError-method}
\alias{scoreReverse,data.frame-method}
\alias{sigmaSqR}
\alias{sigmaSqR,pics-method}
\alias{sigmaSqR,picsList-method}
\alias{sigmaSqR,picsError-method}
\alias{sigmaSqR,data.frame-method}
\alias{sigmaSqF}
\alias{sigmaSqF,pics-method}
\alias{sigmaSqF,picsList-method}
\alias{sigmaSqF,picsError-method}
\alias{sigmaSqF,data.frame-method}
\alias{K}
\alias{K,pics-method}
\alias{K,picsList-method}
\alias{K,picsError-method}
\alias{plot,pics,segReads-method}
\alias{plot,picsError,segReads-method}
\alias{plot,picsList,segReadsList-method}
\alias{maxRange}
\alias{maxRange,pics-method}
\alias{maxRange,picsList-method}
\alias{maxRange,picsList-method}
\alias{maxRange,picsError-method}
\alias{minRange}
\alias{minRange,pics-method}
\alias{minRange,picsList-method}
\alias{minRange,picsList-method}
\alias{minRange,picsError-method}

% I think a few methods are missing, e.g. [[<- and [<-

\title{Estimation of binding site positions} 

\description{
This object contains Estimation of binding site positions and has the following  slots: 
segReadsList, dataType.
}

\usage{
PICS(segReadsList,dataType=NULL, paraEM=NULL, paraPrior=NULL, nCores=1)
}

\arguments{
\item{segReadsList}{This object contains segmentation of Genome}
\item{dataType}{The type of data you are  processing: specified `TF' for transcription factor.}
\item{paraEM}{A list of parameters for the EM algorithm as returned by the \code{setParaEm} function. The default parameters should be good enough for most usages.}
	\describe{
		\item{}{\code{minK}: an \code{integer}, default=1. The minimum number of binding events per region. If the value is 0, the minimum number is automatically calculated.}
		\item{}{\code{maxK}: an \code{integer}, default=15. The maximum number of binding events per region. If the value is 0, the maximum number is automatically calculated.}
		\item{}{\code{tol}: a \code{numeric}, default=1e-4. The tolerance for the EM algorithm.}
		\item{}{\code{B}: an \code{integer}, default=100. The maximum number of iterations to be used.}
		\item{}{\code{mSelect}: a character string specifying the information criteria to be used when selecting the number of binding events. Default="BIC"}
		\item{}{\code{mergePeaks}: a \code{logical} stating whether overlapping binding events should be picked. Default=TRUE}
		\item{}{\code{mapCorrect}: a \code{logical} stating whether mappability profiles should be incorporated in the estimation, i.e: missing reads estimated. Default=TRUE}
	}
\item{paraPrior}{A list of parameters for the prior distribution as returned by the \code{setParaPrior} function. The default parameters should be good enough for most usages.}
	\describe{
		\item{}{\code{xi}: an \code{integer}, default=200. The average DNA fragment size.}
		\item{}{\code{rho}: an \code{integer}, default=1. A variance parameter for the average DNA fragment size distribution.}
		\item{}{\code{alpha}: an \code{integer}, default=20. First hyperparameter of the inverse Gamma distribution for sigma^2 in the PICS model}
		\item{}{\code{beta}: an \code{integer}, default=40000. Second hyperparameter of the inverse Gamma distribution for sigma^2 in the PING model}
		\item{}{\code{lambda}: an \code{integer}, default=0. The precision of the prior for mu used for histone data.}
		\item{}{\code{dMu}: an \code{integer}, default=0. Our best guess for the distance between two neighboring nucleosomes.}
	}
\item{nCores}{An \code{integer}. The number of cores that should be used in parallel by the function.}
}

\section{Methods}{
  \describe{
 \item{code}{\code{signature(x = ``pics'')}: return the error code for each list element (i.e. candidate region) of a PICS object. If the string is empty, there were no errors.}
 \item{plot}{\code{signature(x = ``pics'')}: Plot all regions in the PICS object. This might be long, and should only be used to plot a few regions, so subset the object before plotting.}
 \item{sigmaSqR}{\code{signature(x = ``pics'')}: return the variance parameter of the reverse (R) distribution for each binding event.}
 \item{sigmaSqF}{\code{signature(x = ``pics'')}: return the variance parameter of the forward (F) distribution for each binding event.}
 \item{score}{\code{signature(x = ``pics'')}: return the score for each binding event.} 
 \item{scoreF}{\code{signature(x = ``pics'')}: return the score  of the forward (F)  for each binding event.}
 \item{scoreR}{\code{signature(x = ``pics'')}: return the score  of the forward (R)  for each binding event.}
\item{maxRange}{\code{signature(x = ``pics'')}: return the range maximum.}
\item{minRange}{\code{signature(x = ``pics'')}: return the range minimal.}
 \item{K}{\code{signature(x = ``pics'')}: subset PICS object.}
 \item{wigDensity}{\code{signature(x = ``pics'')}: return the density for each binding event.}
}
}

 \author{ Xuekui Zhang, Arnaud Droit <\email{arnaud.droit@crchuq.ualaval.ca}> and Raphael Gottardo <\email{rgottard@fhcrc.org}>}


\references{
X. Zhang, G. Robertson, M. Krzywinski, K. Ning, A. Droit, S. Jones, and R. Gottardo, ``PICS: Probabilistic Inference for ChIP-seq'' arXiv, 0903.3206, 2009. To appear in Biometrics.
}
\seealso{
\code{\link{pics}}
}
\keyword{models}
