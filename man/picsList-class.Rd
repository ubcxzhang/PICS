\name{picsList-class}
\docType{class}
\alias{picsList}
\alias{picsList-class}
\alias{chromosome,picsList-method}
\alias{mu,picsList-method}
\alias{delta,picsList-method}
\alias{w,picsList-method}
\alias{score,picsList-method}
\alias{density,picsList-method}
\alias{se,picsList-method}
\alias{seF,picsList-method}
\alias{seR,picsList-method}
\alias{[,picsList-method}
\alias{[,picsList,ANY,ANY-method}
\alias{[[,picsList-method}
\alias{[[,picsList,ANY,ANY-method}
\alias{newPicsList,picsList-method}
\alias{length}
\alias{length,picsList-method}
\alias{newPicsList}
\title{The pics class}

\description{
This object is used to gather all parameters from fitting PICS to multiple candidate regions (as returned by the `segmentReads' function). 
The objet contains the following slots: `List', `paraPrior', `paraEM', `minReads', `N', `Nc'. `List' is a list of `pics' or `picsError' objects. `paraPrior' is a list containing the hyperparameters used for the prior, `paraEM' is a list of convergence parameters for the EM, `minReads' is a list containing the minimum number of reads used to fit a region with `PICS', `N' is the total number of reads in the ChIP samples while `Nc' is the total number of reads in the control sample.
}

\section{Accessors}{
The PICS package provide accessors to directly access to most of the parameters/standard errors and chromosomes. 
In the code snippets below, `x' is a `picsList' object. For all accessors, the `picsError' objects are omitted, so that the accessors only return values for the `pics' objects (i.e. all valid binding events).
\describe{
\item{'chromosome(x)'}{Gets the chromosome names of all candidate regions.}
\item{'mu(x)'}{Gets the position estimates of all binding sites identified in all candidate regions.}
\item{'delta(x)'}{Gets the average fragment lengths of all binding sites identified in all candidate regions.}
\item{'sigmaSqF(x)'}{Gets the F peak variances of all binding sites identified in all candidate regions.}
\item{'sigmaSqR(x)'}{Gets the R peak variances of all binding sites identified in all candidate regions.}
\item{'seF(x)'}{Gets the standard errors of all binding site position estimates identified in all candidate regions.}
\item{'seF(x)'}{Gets the standard errors of all F peak modes identified in all candidate regions.}
\item{'seR(x)'}{Gets the standard errors of all R peak modes identified in all candidate regions.}
\item{'score(x)'}{Gets the scores of all binding events identified in all candidate regions.}
}
}

\section{Constructor}{
newPicsList(List, paraEM, paraPrior, minReads, N, Nc)
\describe{
\item{List}{The mixture weights (a vector)}
\item{paraEM}{The binding site positions (a vector)}
\item{paraPrior}{The DNA fragment lengths (a vector)}
\item{N}{The variance parameters for the forward distribution (vector)}
\item{Nc}{The variance parameters for the forward distribution (vector)}
}
}


\section{Methods}{
  \describe{
 \item{[}{\code{signature(x = ``pics'')}: subset PICS object.}
}
}


\section{Methods}{
  \describe{
 \item{length}{\code{signature(x = ``pics'')}: subset PICS object.}
}
}

\arguments{
		\item{object}{An object of class \code{pics}.}	
}


\section{Constructor}{
newPicsList<-function(List, paraEM, paraPrior, minReads, N, Nc) constructs a new `picsList' object with the following arguments. 
\describe{
\item{newPicsList}{}
\item{w}{The mixture weights (a vector)}
\item{mu}{The binding site positions (a vector)}
\item{delta}{The DNA fragment lengths (a vector)}
\item{sigmaSqF}{The variance parameters for the forward distribution (vector)}
\item{sigmaSqR}{The variance parameters for the reverse distribution (vector)}
\item{seMu}{The standard errors for mu (vector)}
\item{seMuF}{The standard errors for muF (vector)}
\item{seMuR}{The standard errors for muR (vector)}
\item{seMuR}{The standard errors for muR (vector)}
\item{score}{The scores for each binding event (vector)}
\item{Nmerged}{The number of peaks that were merged (integer)}
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
}


\examples{
# Here is an example of how to construct such a region
# Typically, you would not do this manually, you would use the pics function to return a 'picsList' that contains a list of 'pics' or 'picsError' object.
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
infMat<-matrix(0)
chr<-"chr1"
range<-c(1000,2000)
# Contructor
#myPICS1<-newPics(w,mu,delta,sigmaSqF,sigmaSqR,seMu,seMuF,seMuR,score,Nmerged,converge,infMat,as.integer(range),chr)
#myPICS2<-newPics(w,mu+1000,delta,sigmaSqF,sigmaSqR,seMu,seMuF,seMuR,score,Nmerged,converge,infMat,as.integer(range),chr)

#minReads<-list(perPeak=2,perRegion=5)
#paraPrior<-list(xi=200,rho=1,alpha=20,beta=40000)
#paraEM<-list(minK=1,maxK=15,tol=10e-6,B=100)
#N<-100
#Nc<-200

#mynewPicsList<-newPicsList(list(myPICS1,myPICS2), paraEM, paraPrior, minReads, as.integer(100), as.integer(200))
# Accessors
# Get the standard error of Mu
#se(mynewPicsList)
# Get the standard error of MuF
#seF(mynewPicsList)
# Get the scores
#score(mynewPicsList)
}


\keyword{models}
