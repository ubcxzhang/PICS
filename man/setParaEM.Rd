\name{setParaEM}
\alias{setParaEM}

\title{
	Function that returns a list of parameters for the EM algorithm that can be used as an argument of PICS.
}

\description{
	This function takes from 0 to 7 EM algorithm parameters as argument, check if they are valid and returns a list to be used in a call to PICS.
}

\usage{
	 setParaEM(minK=1,maxK=15,tol=1e-4,B=100,mSelect="BIC",mergePeaks=TRUE,mapCorrect=TRUE,dataType=NULL)
}

\arguments{
	\item{minK}{ An \code{integer}. The minimum number of binding events per region. If the value is 0, the minimum number is automatically calculated.}
	\item{maxK}{ An \code{integer}. The maximum number of binding events per region. If the value is 0, the maximum number is automatically calculated.}
	\item{tol}{ A \code{numeric}. The tolerance for the EM algorithm.}
	\item{B}{ An \code{integer}. The maximum number of iterations to be used.}
	\item{mSelect}{ A character string specifying the information criteria to be used when selecting the number of binding events.}
	\item{mergePeaks}{ A \code{logical} stating whether overlapping binding events should be picked.}
	\item{mapCorrect}{ A \code{logical} stating whether mappability profiles should be incorporated in the estimation, i.e: missing reads estimated.}
	\item{dataType}{ A \code{character}. If a dataType is set, the algorithm will use the default parameters for this type of data (all the previous arguments will be ignored).}
}

\value{
	Returns a list of parameters to be used in PICS.
}

\author{
Renan Sauteraud
}

\seealso{
	\code{PICS}
}

\keyword{functions}
