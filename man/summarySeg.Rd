\name{summarySeg}
\Rdversion{1.1}
\alias{summarySeg}

\title{
	Summarize a segReadsList object.
}

\description{
	Returns info about a segReadsList object in a \code{data.frame}
	containing the following informations:
	chr : chromosome id
	NF : number of forward reads
	NR : number of reverse reads
	L  : length of segment
	min: start location of segments
	max: end location of segments	
}

\usage{
	summarySeg(seg)
}

\arguments{
	\item{seg}{ An object of class \code{segReadsList}}
}

\value{
	A six columns \code{data.frame}.
}

\author{
	Xuekui Zhang
}

\seealso{
	\code{\linkS4class{segReadsList}}
}

\keyword{functions}