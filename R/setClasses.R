## =========================================================================##
## =========================================================================##
##                    Class definitions and contructors                     ##
## =========================================================================##
## =========================================================================##


setClass("segReads", representation(yF="numeric", yR="numeric", cF="numeric", cR="numeric", map="matrix",chr="character"),
		prototype(yF=numeric(0), yR=numeric(0), cF=numeric(0), cR=numeric(0), map=matrix(0,0,2),chr=character(0)))

## I will need to add some annotations latter on.
setClass("segReadsList", representation(List="list", paraSW="list", N="integer", Nc="integer"),
		prototype(List=list(List=list(0),paraSW=list(step=integer(0),width=integer(0),minReads=integer(0)),N=integer(0),Nc=integer(0))))

### Constructor

segReads<-function(yF, yR, cF, cR, map, chr)
{
  if(!is.vector(yF) || !is.vector(yR) || !is.numeric(yF) || !is.numeric(yR))
  {
    stop("Argument 'yF/yR' must be numeric vectors ", call.=FALSE)
  }
  if((!is.vector(cF) || !is.vector(cR) || !is.numeric(cF) || !is.numeric(cR)) & (!is.null(cF) || !is.null(cR)))
  {
    stop("Argument 'cF/cR' must be numeric vectors ", call.=FALSE)
  }
  
  if(!is.matrix(map))
  {
    stop("Argument 'map' must be a matrix ", call.=FALSE)
  }
  if(!is.character(chr))
  {
    stop("Argument 'chr' must be a character string", call.=FALSE)
  }

  new("segReads", yF=yF, yR=yR, cF=cF, cR=cR, map=map, chr=chr)
}

segReadsList<-function(List,paraSW,N,Nc)
{
  if(!is.list(paraSW) & !all(sapply(paraSW,"is.numeric")))
  {
    stop("Argument 'paraSW' must be a list of numeric arguments", call.=FALSE)
  }
  if(any(lapply(List,"class")!="segReads"))
  {
    stop("Argument 'List' must be a list of segReads arguments", call.=FALSE)
  }
  if(!is.integer(N) | !is.integer(Nc))
  {
    stop("Argument 'N' and 'Nc' must be integers", call.=FALSE)    
  }
  new("segReadsList", List=List, paraSW=paraSW, N=N, Nc=Nc)
}

setClass("pics", representation(estimates="list",score="numeric",scoreF="numeric",scoreR="numeric",Nmerged="numeric",converge="logical",range="numeric",chr="character"),
		prototype(estimates=list(w=numeric(0),mu=numeric(0),delta=numeric(0),sigmaSqF=numeric(0),sigmaSqR=numeric(0),seMu=numeric(0),seMuF=numeric(0),seMuR=numeric(0)),score=numeric(0),scoreF=numeric(0),scoreR=numeric(0),Nmerged=numeric(0),converge=logical(0),range=numeric(0),chr=character(0)))

setClass("picsError", representation(errorCode="character"),prototype(errorCode=character(0)))

## Should I add some annotations?
setClass("picsList", representation(List="list", paraPrior="list", paraEM="list", minReads="list", N="integer", Nc="integer"), 
		prototype(List=list(0),minReads=list(perPeak=integer(0),perRegion=integer(0)), paraPrior=list(xi=double(0),rho=double(0),alpha=double(0),beta=double(0)),paraEM=list(kMax=integer(0),B=integer(0),tol=double(0)),N=integer(0), Nc=integer(0)))

### Constructor
newPics<-function(w,mu,delta,sigmaSqF,sigmaSqR,seMu,seMuF,seMuR,score,scoreF,scoreR,Nmerged,converge,range,chr)
{
  if(!all(is.double(w)))
  {
    stop("Argument 'w' must be numeric ", call.=FALSE)
  }
  if(!all(is.double(mu)))
  {
    stop("Argument 'mu' must be numeric ", call.=FALSE)
  }
  if(!all(is.double(delta)))
  {
    stop("Argument 'delta' must be numeric ", call.=FALSE)
  }
  if(!all(is.double(sigmaSqF)) | !all(is.double(sigmaSqR)))
  {
    stop("Argument 'sigmaSqF/sigmaSqR' must be numeric ", call.=FALSE)
  }
#  if(!all(is.double(seMu)) | !all(is.double(seMuF)) | !all(is.double(seMuR)))
#  {
#    stop("Argument 'seMu/seMuF/seMuR' must be numeric ", call.=FALSE)
#  }
  if(!all(is.double(score)))
  {
    stop("Argument 'score' must be numeric ", call.=FALSE)
  }
  if(!is.numeric(Nmerged))
  {
    stop("Argument 'Nmerged' must be numeric ", call.=FALSE)
  }
  if(!is.logical(converge))
  {
    stop("Argument 'converge' must be logical ", call.=FALSE)
  }
  if(!is.character(chr))
  {
    stop("Argument 'chr' must be a character string", call.=FALSE)
  }
  # if(!all(is.numeric(range)))
  # {
  #   stop("Argument 'range' must be numeric ", call.=FALSE)
  # }  
  new("pics", estimates=list(w=w,mu=mu,delta=delta,sigmaSqF=sigmaSqF,sigmaSqR=sigmaSqR,seMu=seMu,seMuF=seMuF,seMuR=seMuR),converge=converge,score=score,scoreF=scoreF,scoreR=scoreR,Nmerged=Nmerged,range=range,chr=chr)
}

# In case the algorithm does not converge
newPicsError<-function(string)
{
  if(!is.character(string))
  {
    stop("Argument 'errorCode' must be of class character", call.=FALSE)
  }
  new("picsError", errorCode=string)
}

newPicsList<-function(List, paraEM, paraPrior, minReads, N, Nc)
{
  if(!is.list(paraEM) & !all(sapply(paraEM,"is.numeric")))
  {
    stop("Argument 'paraEM' must be a list of numeric arguments", call.=FALSE)
  }
  if(!is.list(paraPrior) & !all(sapply(paraPrior,"is.numeric")))
  {
    stop("Argument 'paraPrior' must be a list of numeric arguments", call.=FALSE)
  }
  if(!is.list(minReads) & !all(sapply(minReads,"is.numeric")))
  {
    stop("Argument 'minReads' must be a list of numeric arguments", call.=FALSE)
  }
  if(!all((lapply(List,"class")=="pics" | lapply(List,"class")=="picsError")))
  {
    stop("Argument 'List' must be a list of 'pics' or 'picsError' arguments", call.=FALSE)
  }
  if(!is.integer(N) | !is.integer(Nc))
  {
    stop("Argument 'N' and 'Nc' must be integers", call.=FALSE)    
  }
  new("picsList", List=List, paraEM=paraEM, paraPrior=paraPrior, minReads=minReads, N=N, Nc=Nc)
}
### Define the classes ##
 setClass("segReadsPE", 
        contains="segReads",
        representation(yFm="numeric", yRm="numeric", cFm="numeric", cRm="numeric"),
	prototype(yFm=numeric(0), yRm=numeric(0), cFm=numeric(0), cRm=numeric(0)))
 setClass("segReadsListPE",
        contains="segReadsList",
        representation(NFm="integer", NRm="integer", NcFm="integer", NcRm="integer"),
 	prototype(list(List=list(0), paraSW=list(islandDepth=integer(0),min_cut=integer(0),max_cut=integer(0), xi=0), NFm=integer(0), NRm=integer(0), NcFm=integer(0), NcRm=integer(0))))
#        representation(yF="numeric", yR="numeric", yFm="numeric", yRm="numeric", cF="numeric",cR="numeric", cFm="numeric", cRm="numeric", map="matrix",chr="character"),
#	prototype(yF=numeric(0), yR=numeric(0), yFm=numeric(0), yRm=numeric(0), cF=numeric(0), cR=numeric(0), cFm=numeric(0), cRm=numeric(0), map=matrix(0,0,2),chr=character(0)))
# setClass("segReadsListPE",
#        representation(List="list", paraSW="list", N="integer", NFm="integer", NRm="integer", Nc="integer", NcFm="integer", NcRm="integer"),
# 	prototype(list(List=list(0), paraSW=list(islandDepth=integer(0),min_cut=integer(0),max_cut=integer(0)), N=integer(0), NFm=integer(0), NRm=integer(0), Nc=integer(0), NcFm=integer(0), NcRm=integer(0))))
 
 ## Constructor ##
 segReadsPE <- function(yF, yR, yFm, yRm, cF, cR, cFm, cRm, map, chr)
 {
 	if(!is.vector(yF) || !is.vector(yR) || !is.numeric(yF) || !is.numeric(yR))
 	{
 		stop("Argument 'yF/yR' must be numeric vectors ", call.=FALSE)
 	}
 	if(!is.vector(yFm) || !is.vector(yRm) || !is.numeric(yFm) || !is.numeric(yRm))
 	{
 		stop("Argument 'yFm/yRm' must be numeric vectors ", call.=FALSE)
 	}
 	
 	if((!is.vector(cF) || !is.vector(cR) || !is.numeric(cF) || !is.numeric(cR)) & (!is.null(cF) || !is.null(cR)))
 	{
 		stop("Argument 'cF/cR' must be numeric vectors ", call.=FALSE)
 	}
 	if((!is.vector(cFm) || !is.vector(cRm) || !is.numeric(cFm) || !is.numeric(cRm)) & (!is.null(cFm) || !is.null(cRm)))
 	{
 		stop("Argument 'cFm/cRm' must be numeric vectors ", call.=FALSE)
 	}
 	
 	if(!is.matrix(map))
 	{
 		stop("Argument 'map' must be a matrix ", call.=FALSE)
 	}	
 	new("segReadsPE", yF=yF, yR=yR, yFm=yFm, yRm=yRm, cF=cF, cR=cR, cFm=cFm, cRm=cRm, map=map, chr=chr)
 }
 
 segReadsListPE<-function(List, paraSW, N, NFm, NRm, Nc, NcFm, NcRm)
 {
 	if(!is.list(paraSW) & !all(sapply(paraSW,"is.numeric")))
 	{
 		stop("Argument 'paraSW' must be a list of numeric arguments", call.=FALSE)
 	}
 	if(any(lapply(List,"class")!="segReadsPE"))
 	{
 		stop("Argument 'List' must be a list of segReadsPE arguments", call.=FALSE)
 	}
 	if(!is.integer(N) | !is.integer(Nc))
 	{
 		stop("Argument 'N' and 'Nc' must be integers", call.=FALSE)    
 	}
 	if(!is.integer(NFm) | !is.integer(NRm))
 	{
 		stop("Argument 'NFm' and 'NRm' must be integers", call.=FALSE)    
 	}
 	if(!is.integer(NcFm) | !is.integer(NcRm))
 	{
 		stop("Argument 'NcFm' and 'NcRm' must be integers", call.=FALSE)    
 	}
 	
 	new("segReadsListPE", List=List, paraSW=paraSW, N=N, NFm=NFm, NRm=NRm, Nc=Nc, NcFm=NcFm, NcRm=NcRm)
 }
