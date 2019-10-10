\name{segReads}
\alias{segReads}
\alias{segReads-class}
\alias{map}
\alias{map,segReads-method}
\alias{map,segReadsList-method}
\alias{\%in\%,ANY,segReads-method}


\title{Segment the genome into candidate regions} 

\description{
Pre-process bidirectional aligned reads data from a single ChIP-Seq experiment to detect candidate regions with a minimum number of forward and reverse reads. These candidate regions will then be processed by PICS.
}


\section{Methods}{
  \describe{
 \item{map}{\code{signature(x = ``pics'')}: subset PICS object.}
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