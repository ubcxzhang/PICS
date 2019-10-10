\name{plot-FDR}
\docType{methods}
\alias{FDR}
\alias{plot,picsList,picsList-method}

\title{FDR plot for PICS} 

\description{
This method plots an FDR curve showing the FDR as a function of the PICS scores.
}

\usage{
\S4method{plot}{picsList,picsList}(x, y, filter=NULL, h=.1, \dots)
}

\arguments{
\item{x}{A \code{picsList} object as returned by the function \code{PICS} run on the treatment data.}
\item{y}{A \code{picsList} object as returned by the function \code{PICS} run on the control data.}
\item{filter}{A list of ranges for filtering regions based on \code{PICS} parameters. By default filter is set to 'NULL' and all regions are used.}
\item{h}{A value between 0 and 1, representing the desired FDR. This simply draws a horizontal line at the given value.}
\item{\dots}{Further graphical parameters passed to the generic function \code{plot}.}
}

% \note{
% }

\author{ Xuekui Zhang, Arnaud Droit <\email{arnaud.droit@crchuq.ualaval.ca}> and Raphael Gottardo <\email{rgottard@fhcrc.org}>}


\references{
X. Zhang, G. Robertson, M. Krzywinski, K. Ning, A. Droit, S. Jones, and R. Gottardo, "PICS: Probabilistic Inference for ChIP-seq," Biometrics, iss. In press, 2010.
}

\seealso{
\code{\link{PICS}}
}

%\examples{
%}

\keyword{graphs}