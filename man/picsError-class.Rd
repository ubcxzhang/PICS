\name{picsError-class}
\docType{class}
\alias{picsError}
\alias{picsError-class}
\alias{chromosome,picsError-method}
\alias{mu,picsError-method}
\alias{delta,picsError-method}
\alias{w,picsError-method}
\alias{score,picsError-method}
\alias{se,picsError-method}
\alias{seF,picsError-method}
\alias{seR,picsError-method}
\alias{newPicsError}
\alias{newPicsError,picsError-method}

\title{The pics class}

\description{
This object is used to return an error code when the PICS function failed to return a valid set of estimates for a candidate regions. This could be due to non-convergence of the EM algorithm, a singular information matrix, or a number of reads below the limit specified by the user. All of these are typically due to too few reads in the region and do not affect the rest of the analysis, as such regions would most likely be labelled as false positives. 
}

\section{Accessors}{
All of the accessors defined for a `pics' object still work for a `picsError' object but will simply return a NULL pointer.
}

\section{Constructor}{
\describe{
newPicsError(string) where `string' is the error code.
}
}

\section{Constructor}{
newPicsError<-function(string)
\describe{
\item{string}{The mixture weights (a vector)}
}
}

\author{ Xuekui Zhang, Arnaud Droit <\email{arnaud.droit@crchuq.ualaval.ca}> and Raphael Gottardo <\email{rgottard@fhcrc.org}>}


\references{
X. Zhang, G. Robertson, M. Krzywinski, K. Ning, A. Droit, S. Jones, and R. Gottardo, ``PICS: Probabilistic Inference for ChIP-seq'' arXiv, 0903.3206, 2009. To appear in Biometrics.
}

\seealso{
\code{\link{pics}}
}


\examples{
# Here is an example on how to construct such a picsError object
# Typically, you would not do this manually, you would use the pics function to return a 'picsList' that contains a list of 'pics' or 'picsError' object.
# Contructor
myPicsError<-newPicsError("Singular information matrix")
# Accessors
# Get the standard error of Mu
se(myPicsError)
# Get the standard error of MuF
seF(myPicsError)
# Get the scores
score(myPicsError)
}


\keyword{models}