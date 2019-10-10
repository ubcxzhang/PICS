\name{makeRangedDataOutput}
\alias{makeRangedDataOutput}

\title{Create a RangedData object from a PICS output}
\description{
Create a list of `RangedData' objects from a `pics' object. The resulting RangedData object can then be analyzed with the `IRanges' packages and/or exported to bed/wig files with the `rtracklayer' package.
}
\usage{
makeRangedDataOutput(obj, type="fixed", filter=list(delta=c(0,Inf),se=c(0,Inf),sigmaSqF=c(0,Inf),sigmaSqR=c(0,Inf),score=c(0,Inf)),length=100)
}

\arguments{
  \item{obj}{An object of class `picsList' as returned by `PICS' when running it on the IP/Control data.}
  \item{type}{The type of intervals to be created. The different types are `bed', `wig', `ci' and `fixed'. See details for more info.}
  \item{filter}{A list of filters to be used before computing the FDR. By default all regions are included, see details for more info on how to specify the filters.}
  \item{length}{The length to be used for the fixed type `RangedData', see details.}
}
\details{
`bed' will generate intervals from the forward peak max to the reverse peak max. `wig' will generate a density profile for the forward and reverse reads. `bed' and `wig' types should be used to be exported to wig/bed files to be used with the UCSC genome browser. `ci' corresponds to the binding site estimates +/-3*se, while `fixed' corresponds to the binding site estimates +/-3*length. `bed' and `wig' files can be exported using the `export' function fo the `rtracklayer' package.
}
\value{
An object of type `RangedData'.
}
\references{
X. Zhang, G. Robertson, M. Krzywinski, K. Ning, A. Droit, S. Jones, and R. Gottardo, ``PICS: Probabilistic Inference for ChIP-seq'' arXiv, 0903.3206, 2009. To appear in Biometrics.
}
\author{ Xuekui Zhang, Arnaud Droit <\email{arnaud.droit@crchuq.ualaval.ca}> and Raphael Gottardo <\email{rgottard@fhcrc.org}>}
\seealso{export}
\examples{


\dontrun{
  rdBed<-makeRangedDataOutput(pics,type="bed",filter=list(delta=c(50,Inf),se=c(0,50),sigmaSqF=c(0,22500),sigmaSqR=c(0,22500),score=c(10,Inf)))
  export(rbBed,"myfile.bed")
  rdBed<-makeRangedDataOutput(pics,type="wig",filter=list(delta=c(50,Inf),se=c(0,50),sigmaSqF=c(0,22500),sigmaSqR=c(0,22500),score=c(10,Inf)))
  export(rbBed,"myfile.wig")}
  
}


\keyword{process}
\keyword{data}
