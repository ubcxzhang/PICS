\name{pics-class}
\docType{class}
\alias{newPics}
\alias{chromosome,pics-method}
\alias{chromosome,data.frame-method}
\alias{mu,pics-method}
\alias{mu,data.frame-method}
\alias{delta,pics-method}
\alias{delta,data.frame-method}
\alias{w,pics-method}
\alias{w,data.frame-method}
\alias{score,pics-method}
\alias{score,data.frame-method}
\alias{se,pics-method}
\alias{se,data.frame-method}
\alias{seF,pics-method}
\alias{seF,data.frame-method}
\alias{seR,pics-method}
\alias{seR,data.frame-method}
\alias{newPics,pics-method}
\alias{seF}
\alias{seR}
\alias{se}
\alias{w}
\alias{delta}
\alias{chromosome}
\alias{mu}

\title{The pics class}

\description{
This object is used to gather all parameters from fitting PICS to a single candidate region. 
The objet contains the following slots: 'estimates', 'infMat', 'Nmerged', 'converge', 'chr'. 'estimates' is a list containing all parameters estimates as well as standard errors. 'infMat' is the Cholesky decomposition of the information matrix, 'converge' is a logical value indicating whether the EM algorithm has converged, while 'chr' is a character string corresponding to a candidate region's chromosome. 'Nmerged' gives the number of binding events that were merged; binding events that overlap are merged (see the cited paper below for details).
}

\section{Accessors}{
The PICS package provide accessors to directly access to most of the parameters/standard errors and chromosome. In the code snippets below, 'x' is a 'pics' object.
\describe{
\item{`chromosome(x)'}{Gets the chromosome name of the candidate region.}
\item{`mu(x)'}{Gets the position estimates of all binding sites identified in the region.}
\item{`delta(x)'}{Gets the average fragment lengths of all binding sites identified in the region.}
\item{`sigmaSqF(x)'}{Gets the F peak variances of all binding sites identified in the region.}
\item{`sigmaSqR(x)'}{Gets the R peak variances of all binding sites identified in the region.}
\item{`seF(x)'}{Gets the standard errors of all binding site position estimates identified in the region.}
\item{`seF(x)'}{Gets the standard errors of all F peak modes identified in the region.}
\item{`seR(x)'}{Gets the standard errors of all R peak modes identified in the region.}
 \item{score}{\code{signature(x = "pics")}: return the score for each binding event.} 
 \item{scoreF}{\code{signature(x = "pics")}: return the score  of the forward (F)  for each binding event.}
 \item{scoreR}{\code{signature(x = "pics")}: return the score  of the forward (R)  for each binding event.}
}
}

\section{Constructor}{
newPics(w,mu,delta,sigmaSqF,sigmaSqR,seMu,seMuF,seMuR,score,Nmerged,converge,infMat,chr) construct a new 'pics' object with the following arguments:

\describe{
\item{w}{The mixture weights (a vector)}
\item{mu}{The binding site positions (a vector)}
\item{delta}{The DNA fragment lengths (a vector)}
\item{sigmaSqF}{The variance parameters for the forward distribution (vector)}
\item{sigmaSqR}{The variance parameters for the forward distribution (vector)}
\item{seMu}{The standard errors for mu (vector)}
\item{seMuF}{The standard errors for muF (vector)}
\item{seMuR}{The standard errors for muR (vector)}
\item{seMuR}{The standard errors for muR (vector)}
\item{score}{The scores for each binding event (vector)}
\item{Nmerged}{The number of peaks that got merged (integer)}
\item{converge}{A logical value, TRUE, if the EM as converged}
\item{infMat}{The information matrix}
\item{chr}{The chromosome for the region}
}
}

\author{ Xuekui Zhang, Arnaud Droit <\email{arnaud.droit@crchuq.ualaval.ca}> and Raphael Gottardo <\email{rgottard@fhcrc.org}>}


\references{
X. Zhang, G. Robertson, M. Krzywinski, K. Ning, A. Droit, S. Jones, and R. Gottardo, ``PICS: Probabilistic Inference for ChIP-seq'' arXiv, 0903.3206, 2009. To appear in Biometrics.
}
\seealso{
\code{\link{pics}}
\code{\link{picsError}}
}


\examples{
# Here is an example of how to construct such a region.
# Typically, you would not do this manually, you would use the pics function to return a 'picsList' that contains a list of 'pics' or a 'picsError' object.
w<-1
mu<-10000
delta<-150
sigmaSqF<-5000
sigmaSqR<-5000
seMu<-10
seMuF<-10
seMuR<-10
score<-5
Nmerged<-0
converge<-TRUE
chr<-"chr1"
range<-c(1000,2000)
# Contructor
#myPICS<-newPics(w,mu,delta,sigmaSqF,sigmaSqR,seMu,seMuF,seMuR,score,Nmerged,as.integer(range),chr)
}


\keyword{models}