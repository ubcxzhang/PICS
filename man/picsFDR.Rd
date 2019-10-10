\name{picsFDR}
\alias{picsFDR}

\title{
	Estimate the FDR.
}

\description{
	Estimate the false detection rate for an object of class \code{pics} or \code{picsList}.
}

\usage{
	picsFDR(picsIP,picsCont,filter=list(delta=c(0,Inf),se=c(0,Inf),sigmaSqF=c(0,Inf),sigmaSqR=c(0,Inf)))
}

\arguments{
	\item{picsIP}{An object of class \code{pics} or \code{picsList} containing the informations for the IP reads}.
	\item{picsCont}{An object of class \code{pics} or \code{picsList} containing the informations for the control reads}.
	\item{filter}{filter}{A list of ranges for filtering regions based on \code{PICS} parameters. By default filter is set to 'NULL' and all regions are used.}
		\describe{
		\item{delta}{Length of the binding sites.}
		\item{se}{Standard error.}
		\item{sigmaSqF}{Forward peak variance}
		\item{sigmaSqR}{Reverse peak variance}
		}
}

\value{
	A 3 columns \code{data.frame} with the following columns:
	FDR, score, N.
}

\author{
	Xuekui Zhang
}

\seealso{
	\code{\linkS4class{picsList}}
	\code{\linkS4class{pics}}
}

\keyword{functions}