\name{segReadsList}
\docType{class}
\alias{segReadsList}
\alias{segReadsList-class}
\alias{[,segReadsList-method}
\alias{[,segReadsList,ANY,ANY-method}
\alias{[[,segReadsList-method}
\alias{[[,segReadsList,ANY,ANY-method}
\alias{length,segReadsList-method}

\title{Segment the genome into candidate regions} 

\description{
Pre-process bidirectional aligned reads data from a single ChIP-Seq experiment to detect candidate regions with a minimum number of forward and reverse reads. These candidate regions will then be processed by PICS.}




\section{Methods}{
  \describe{
 \item{[}{\code{signature(x = ``pics'')}: subset gadem object.}
 \item{[[}{\code{signature(x = ``pics'')}: subset gadem object.}
	}
}


\section{Methods}{
  \describe{
 \item{length}{\code{signature(x = ``pics'')}: subset PICS object.}
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
