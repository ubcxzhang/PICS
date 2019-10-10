.useRanges <- function() {
  .Deprecated(msg = "Storing reads only by their start positions is deprecated. Please use a range representation like GRanges.")
}

  
## All the setAs method force an object to belong to a chosen class
## setAs(from, to, def)
##value<-as(object, Class)
##setMethod("unique", "GenomeData",
##function(x,incomparables = FALSE, ...)
##{
##  GenomeData(lapply(x,function(x){lapply(x,unique)}))
##})

setAs("picsList", "RangedData",
function(from) 
{
  makeRangedDataOutput(from, type="bed", filter=list(delta=c(50,300),se=c(0,50),sigmaSqF=c(0,22500),sigmaSqR=c(0,22500),score=c(1,Inf)),length=100)
}
)

##Transform a list of reads into a GRanges object
#setAs("list", "GRanges", function(from, chrs)
#{
#  idx<-which(names(reads)=="P")
#  if(is.null(idx))
#    stop("The list to coerce should have an attribute 'P'")
#  idx2 <- ( from$P$"pos.-" >  from$P$"pos.+")
#  PE_data<-from$P[idx2,]
#  #The reads are saved in a file containing 2 R objects: reads & chrs
#  GRanges(IRanges(start=PE_data$"pos.+" , end=PE_data$"pos.-"), strand="*", seqnames=chrs)
#}


#Transform a .bed like data.frame into a GRanges object
setAs("data.frame", "GRanges",function(from)
{
  if(length(from)<4)
    stra<-"*"
  else
    stra<-from[,4]
  GRanges(ranges=IRanges(start=from[,2],end=from[,3]), seqnames=from[,1], strand=stra)
})


#I need a concatenate method to add pics objects


#setAs("RangedData", "GenomeData",
#function(from)
#{
#  .useRanges()
#  readStart <- ifelse(strand(from) == "-",end(from),start(from))
#  alignLocs <-
#  split(data.frame(position = readStart, strand = strand(from)),space(from)[drop=TRUE])
#  GenomeData(lapply(alignLocs,function(df) with(df, split(position, strand))[c("-", "+")]))
#}
#)
#
#setAs("data.frame", "GenomeData",            
#function(from) 
#{
#  from<-as(from,"RangedData")
#  readStart <- ifelse(strand(from) == "-",end(from),start(from))
#  alignLocs <-
#  split(data.frame(position = readStart, strand = strand(from)),
#  space(from)[drop=TRUE])
#  GenomeData(lapply(alignLocs,function(df) with(df, split(position, strand))[c("-", "+")]))
#}
#)

setAs("picsList", "data.frame",
function(from)
{
  ans <- data.frame(ID=rep(1:length(from),K(from)),chr=chromosome(from),w=w(from), mu=mu(from),
  delta=delta(from), sigmaSqF=sigmaSqF(from), sigmaSqR=sigmaSqR(from),se=se(from),
  score=score(from), scoreF=scoreForward(from),scoreR=scoreReverse(from),
  minRange=minRange(from), maxRange=maxRange(from))
  ans$chr	<- as.character(ans$chr)
  ans		<- ans[is.finite(ans$mu),]
  return(ans)
}
)


## show and summary methods
setMethod("show", "segReads",
          function(object)
      {
          cat("Object of class ",as.character(class(object)),"\n")
          cat("This object has the following slots: \n")
          cat(paste(names(getSlots(class(object))),collapse=", "),"\n")
          #cat("yR, yF, cF, cR, map\n")
      })

setMethod("show", "segReadsList",
          function(object)
      {
          cat("Object of class",as.character(class(object)),"\n")
          cat("This object has the following slots: \n")
          cat(paste(names(getSlots(class(object))),collapse=", "),"\n")
          #cat("List, paraSW, N, Nc\n")
          cat("List is a list of 'segReads' ojects, each of which has the following slots:\n")
          cat("yR, yF, cR, cF, map, chr\n")
      })

setMethod("show", "pics",
      function(object)
      {
        cat("Object of class ",class(object),"\n")
        cat("This object has the following slots: \n")
        cat(paste(names(getSlots(class(object))),collapse=", "),"\n")
        #cat("estimates, score, scoreF, scoreR, Nmerged, converge, chr, range\n")     
        })

setMethod("show", "picsError",
          function(object)
          {
            cat("Object of class ",class(object),"\n")
            cat("This object has the following slot: \n")
            cat(paste(names(getSlots(class(object))),collapse=", "),"\n")
            #cat("errorCode\n")     
          })

setMethod("show", "picsList",
          function(object)
          {
            cat("Object of class ",class(object),"\n")
            cat("This object has the following slots: \n")
            cat(paste(names(getSlots(class(object))),collapse=", "),"\n")
            #cat("List, paraEM, paraPrior, minReads, N, Nc\n")
            cat("List is a list of 'pics' or picsError ojects\n")
          })


#setGeneric("score", function(x, ...) standardGeneric("score"))
setMethod("score", "pics",
          function(x)
          {
            return(x@score)
})
          
setMethod("score", "picsError",
          function(x)
          {
            return(NULL)
})

setMethod("score", "picsList",
          function(x)
          {
            ans<-.Call("getScore", x@List, PACKAGE="PICS");
            return(ans)
            
          })

setMethod("score", "data.frame",
		function(x)
		{
		  return(x$score)
		})


setGeneric("minRange", function(x, ...) standardGeneric("minRange"))
setMethod("minRange", "pics",
          function(x)
          {
            return(x@range[1])
})
          
setMethod("minRange", "picsError",
          function(x)
          {
            return(NULL)
})

setMethod("minRange", "picsList",
          function(x)
          {
            ans<-.Call("getMin", x@List, PACKAGE="PICS");
            return(ans)
          })

setGeneric("maxRange", function(x, ...) standardGeneric("maxRange"))
setMethod("maxRange", "pics",
          function(x)
          {
            return(x@range[2])
})
          
setMethod("maxRange", "picsError",
          function(x)
          {
            return(NULL)
})

setMethod("maxRange", "picsList",
          function(x)
          {
            ans<-.Call("getMax", x@List, PACKAGE="PICS");
            return(ans)
          })


setGeneric("scoreReverse", function(x, ...) standardGeneric("scoreReverse"))
setMethod("scoreReverse", "pics",
          function(x)
          {
            return(x@scoreR)
})
          
setMethod("scoreReverse", "picsError",
          function(x)
          {
            return(NULL)
})

setMethod("scoreReverse", "picsList",
          function(x)
          {
            ans<-.Call("getScoreR", x@List, PACKAGE="PICS");
            return(ans)            
})

setMethod("scoreReverse", "data.frame",
		function(x)
	  	{
		  return(x$scoreR)
})
          
setGeneric("scoreForward", function(x, ...) standardGeneric("scoreForward"))
setMethod("scoreForward", "pics",
          function(x)
          {
            return(x@scoreF)
})
          
setMethod("scoreForward", "picsError",
          function(x)
          {
            return(NULL)
})

setMethod("scoreForward", "picsList",
          function(x)
          {
            ans<-.Call("getScoreF", x@List, PACKAGE="PICS");
            return(ans)
            
})

setMethod("scoreForward", "data.frame",
		function(x)
		{
			return(x$scoreF)
})

setGeneric("chromosome", function(x, ...) standardGeneric("chromosome"))
setMethod("chromosome", "pics",
          function(x)
          {
            return(rep(x@chr,length(x@estimates$w)))
})
setMethod("chromosome", "picsError",
          function(x)
          {
            return(NULL)
          }
)
setMethod("chromosome", "picsList",
          function(x)
          {
            ans<-.Call("getChr", x@List, PACKAGE="PICS");
            return(ans)
          }
)

setMethod("chromosome", "data.frame",
		function(x)
		{
			return(x$chr)
		})



setGeneric("map", function(x, ...) standardGeneric("map"))
setMethod("map", "segReads",
          function(x)
          {
            if(is.null(x@map) | (nrow(x@map)==0))
            {
              return(0);
            }
            else
            {
              n<-nrow(x@map)
              m<-min(x@yF[1],x@yR[1],x@map[1,1]);M<-max(tail(x@yF,1),tail(x@yR,1),x@map[n,2]);
              return(sum(diff(t(x@map)))/max(M-m,1));
            }
})

setMethod("map", "segReadsList",
          function(x)
          {
            ans<-.Call("getMap", x@List, PACKAGE="PICS");
            return(ans)
          }
)


setGeneric("se", function(x, ...) standardGeneric("se"))
setMethod("se", "pics",
          function(x)
          {
            return(x@estimates$seMu)
})

setMethod("se", "picsError",
          function(x)
          {
            return(NULL)
})

setMethod("se", "picsList",
          function(x)
          {
              ans<-.Call("getVector", x@List, as.integer(5), PACKAGE="PICS");
              return(ans)
          }
)

setMethod("se", "data.frame",
		function(x)
		{
			return(x$se)
		})

setGeneric("seF", function(x, ...) standardGeneric("seF"))
setMethod("seF", "pics",
          function(x)
          {
            return(x@estimates$seMuF)
})

setMethod("seF", "picsError",
          function(x)
          {
            return(NULL)
})

setMethod("seF", "picsList",
          function(x)
          {
              ans<-.Call("getVector", x@List, as.integer(6), PACKAGE="PICS");
              return(ans)
          }
)

setMethod("seF", "data.frame",
		function(x)
		{
			return(x$seF)
		})

setGeneric("seR", function(x, ...) standardGeneric("seR"))
setMethod("seR", "pics",
          function(x)
          {
            return(x@estimates$seMuR)
})

setMethod("seR", "picsError",
          function(x)
          {
            return(NULL)
})


setMethod("seR", "picsList",
          function(x)
          {
              ans<-.Call("getVector", x@List, as.integer(7), PACKAGE="PICS");
              return(ans)
          }
)

setMethod("seR", "data.frame",
		function(x)
		{
			return(x$seR)
		})


setGeneric("sigmaSqF", function(x, ...) standardGeneric("sigmaSqF"))
setMethod("sigmaSqF", "pics",
          function(x)
          {
            return(x@estimates$sigmaSqF)
})

setMethod("sigmaSqF", "picsError",
          function(x)
          {
            return(NULL)
})

setMethod("sigmaSqF", "picsList",
          function(x)
          {
              # temp<-lapply(x@List,"mu")
              # return(unlist(temp))
              ans<-.Call("getVector", x@List, as.integer(3), PACKAGE="PICS");
              return(ans)
          }
)

setMethod("sigmaSqF", "data.frame",
		function(x)
		{
			return(x$sigmaSqF)
		})

setGeneric("sigmaSqR", function(x, ...) standardGeneric("sigmaSqR"))
setMethod("sigmaSqR", "pics",
          function(x)
          {
            return(x@estimates$sigmaSqR)
})

setMethod("sigmaSqR", "picsError",
          function(x)
          {
            return(NULL)
})

setMethod("sigmaSqR", "picsList",
          function(x)
          {
              # temp<-lapply(x@List,"mu")
              # return(unlist(temp))
              ans<-.Call("getVector", x@List, as.integer(4), PACKAGE="PICS");
              return(ans)
          }
)

setMethod("sigmaSqR", "data.frame",
		function(x)
		{
			return(x$sigmaSqR)
		})

setGeneric("delta", function(x, ...) standardGeneric("delta"))
setMethod("delta", "pics",
          function(x)
          {
            return(x@estimates$delta)
})

setMethod("delta", "picsError",
          function(x)
          {
            return(NULL)
})

setMethod("delta", "picsList",
          function(x)
          {
              # temp<-lapply(x@List,"mu")
              # return(unlist(temp))
              ans<-.Call("getVector", x@List, as.integer(2), PACKAGE="PICS");
              return(ans)
          }
)

setMethod("delta", "data.frame",
		function(x)
		{
			return(x$delta)
		})

setGeneric("mu", function(x, ...) standardGeneric("mu"))
setMethod("mu", "pics",
          function(x)
          {
            return(x@estimates$mu)
})

setMethod("mu", "picsError",
          function(x)
          {
            return(NULL)
})

setMethod("mu", "picsList",
          function(x)
          {
              # temp<-lapply(x@List,"mu")
              # return(unlist(temp))
              ans<-.Call("getVector", x@List, as.integer(1), PACKAGE="PICS");
              return(ans)
          }
)

setMethod("mu", "data.frame",
		function(x)
		{
			return(x$mu)
		})

setGeneric("w", function(x, ...) standardGeneric("w"))
setMethod("w", "pics",
          function(x)
          {
            return(x@estimates$w)
})

setMethod("w", "picsError",
          function(x)
          {
            return(NULL)
})

setMethod("w", "picsList",
          function(x)
          {
              ans<-.Call("getVector", x@List, as.integer(0), PACKAGE="PICS");
              return(ans)
          }
)

setGeneric("K", function(x, ...) standardGeneric("K"))
setMethod("K", "pics",
          function(x)
          {
            return(length(x@estimates$w))
})

setMethod("K", "picsError",
          function(x)
          {
            return(0)
})

setMethod("K", "picsList",
          function(x)
          {
              ans<-.Call("getK", x@List, PACKAGE="PICS");
              return(ans)
          }
)

setGeneric("code", function(x, ...) standardGeneric("code"))
setMethod("code", "pics",
          function(x)
          {
            return("")
})

setMethod("code", "picsError",
          function(x)
          {
            return(x@errorCode)
})

setMethod("code", "picsList",
          function(x)
          {
              temp<-lapply(x@List,"code")
              return(unlist(temp))
          }
)

setMethod("length", "picsList",
          function(x)
          {
            return(length(x@List))
})

setMethod("length", "segReadsList",
          function(x)
          {
            return(length(x@List))
})

setGeneric("wigDensity", function(x, ...) standardGeneric("wigDensity"))
setMethod("wigDensity", "pics",
          function(x,strand="+",step=10,sum=FALSE,filter=NULL,scale=TRUE)
          {
            
            # Check that all filters are passed
            missingNames<-!c("delta","sigmaSqF","sigmaSqR","se","seF","seR","score")%in%names(filter)
            filter[c("delta","sigmaSqF","sigmaSqR","se","seF","seR","score")[missingNames]]<-list(c(0,Inf))
            if(strand=="+")
            {
              strand<-1
            }
            else if(strand=="-")
            {
              strand<--1
            }
            else if(strand=="*")
            {
              strand<-0
            }
            else
            {
              stop("Strand must be either '+', '-' or '*'")
            }            
            strand<-as.double(paste(strand,"1",sep=""))
            ans<-.Call("getDensity", x, strand, step, filter, sum, scale, PACKAGE="PICS")
            return(ans);
          }
)

setMethod("wigDensity", "picsList",
          function(x,strand="+",step=10,sum=FALSE,filter=NULL,scale=TRUE)
          {
            # Check that all filters are passed
            missingNames<-!c("delta","sigmaSqF","sigmaSqR","se","seF","seR","score")%in%names(filter)
            filter[c("delta","sigmaSqF","sigmaSqR","se","seF","seR","score")[missingNames]]<-list(c(0,Inf))

            if(strand=="+")
            {
              strand<-1
            }
            else if(strand=="-")
            {
              strand<--1
            }
            else if(strand=="*")
            {
              strand<-0
            }
            else
            {
              stop("Strand must be either '+', '-' or '*'")
            }
            ans<-.Call("getDensityList", x, strand, step, filter, sum, scale, PACKAGE="PICS")
            return(ans);
          }
)

setMethod("wigDensity", "picsError",
          function(x,strand=NULL,step=NULL,sum=NULL,filter=NULL)
          {
            return(NULL)
          }
)
setMethod("summary", "segReadsList",
          function(object)
      {
          cat("** Experiment information ** \n")
          cat("Chromosomes interogated: ")
          cat(unique(unlist(lapply(object@List,function(obj){obj@chr}))),"\n")
          cat("Number of reads")
          cat(" in IP: ",object@N," and in control: ",object@Nc,"\n")
          cat("** Segmentation parameters ** \n")
          cat("The following settings were used:\n")          
          cat("  Sliding window half width: ", object@paraSW$width,"\n")
          cat("  Step size: ", object@paraSW$step,"\n")          
          cat("  Minimum number of reads: ", object@paraSW$minReads,"\n")
          cat("** Segmentation summary ** \n")                    
          cat("Number of segmented regions:",length(object@List),"\n")
          cat("Summary on the number of Forward/Reverse reads per region:\n")
          cat("  Forward:\n") 
          cat("  ")
          tempF<-lapply(object@List,function(obj){length(obj@yF)})
          print(summary(as.integer(unlist(tempF))))
          cat("  Reverse:\n") 
          cat("  ")
          tempR<-lapply(object@List,function(obj){length(obj@yR)})
          print(summary(as.integer(unlist(tempR))))
          cat("Summary on the number of control Forward/Reverse reads per region:\n")
          cat("  Forward:\n") 
          cat("  ")
          tempF<-lapply(object@List,function(obj){length(obj@cF)})
          print(summary(as.integer(unlist(tempF))))
          cat("  Reverse:\n") 
          cat("  ")
          tempR<-lapply(object@List,function(obj){length(obj@cR)})
          print(summary(as.integer(unlist(tempR))))                    
          tempMap<-map(object)
          cat("** Mappability summary **\n")
          cat("Non mappable intervals cover an average ", mean(unlist(tempMap)),"% of all regions \n")                  
      })


setMethod("summary", "segReads",
      function(object)
      {
        m<-min(object@yF[1],object@yR[1])
        M<-max(tail(object@yF,1),tail(object@yR,1))
        cat("** Region summary ** \n")
        cat("Summary on Forward reads:\n")
        print(summary(object@yF,digits=100))
        cat("Summary on Reverse reads:\n")
        print(summary(object@yR,digits=100))
        cat("Summary on control Forward reads:\n")
        print(summary(object@cF,digits=100))
        cat("Summary on control Reverse reads:\n")
        print(summary(object@cR,digits=100))
        cat("Non mappable intervals cover ", sum(diff(t(object@map)))/(M-m),"% of the region \n")
      })


setMethod("[","segReadsList",
		function(x,i, j,..., drop=FALSE)
		{
			if(missing(i))
			{
				return(x)
			}
			if(!missing(j))
			{
			  stop("incorrect number of dimensions")
			}
      else
      {
        segReadsList(x@List[i],x@paraSW,x@N,x@Nc)
      }
		})

setMethod("[[","segReadsList",
    function(x, i, j, ..., exact = TRUE)
    {
      if(length(i) != 1)
      {
        stop("subscript out of bounds (index must have length 1)")
      }
      if(missing(i))
      {
        return(x)
      }
      if(!missing(j))
      {
        stop("incorrect number of dimensions")
      }
      x@List[[i]]
})


setMethod("[","segReadsListPE",
                function(x,i, j,..., drop=FALSE)
                {
                        if(missing(i))
                        {
                                return(x)
                        }
                        if(!missing(j))
                        {
                          stop("incorrect number of dimensions")
                        }
      else
      {
        segReadsListPE(x@List[i],x@paraSW,x@N,x@NFm,x@NRm,x@Nc,x@NcFm,x@NcRm)
      }
                })

#Same as for segReadsList. Should not be needed
setMethod("[[","segReadsListPE",
    function(x, i, j, ..., exact = TRUE)
    {
      if(length(i) != 1)
      {
        stop("subscript out of bounds (index must have length 1)")
      }
      if(missing(i))
      {
        return(x)
      }
      if(!missing(j))
      {
        stop("incorrect number of dimensions")
      }
      x@List[[i]]
})









setMethod("[","picsList",
		function(x,i, j,..., drop=FALSE)
		{
			if(missing(i))
			{
				return(x)
			}
			if(!missing(j))
			{
			  stop("incorrect number of dimensions")
			}
      else
      {
        newPicsList(x@List[i], x@paraEM, x@paraPrior, x@minReads, x@N, x@Nc)        
      }
		})

setMethod("[[","picsList",
    function(x, i, j, ..., exact = TRUE)
    {
      if(length(i) != 1)
      {
        stop("subscript out of bounds (index must have length 1)")
      }
      if(missing(i))
      {
        return(x)
      }
      if(!missing(j))
      {
        stop("incorrect number of dimensions")
      }
      x@List[[i]]
})


setMethod("summary", "picsList",
          function(object)
          {
            cat("** Experiment information ** \n")
            cat("Chromosomes interogated: ")
            cat(unique(chromosome(object)),"\n")
            cat("Number of reads:")
            cat("In IP: ",object@N," in control: ",object@Nc,"\n")
            cat("** Prior parameters ** \n")
            cat("The following settings were used:\n")          
            cat("  Hyper parameters for the fragment length distribution:\n")
            cat("  xi, rho, alpha, beta: ", object@paraPrior$xi,",", object@paraPrior$rho, ",", object@paraPrior$alpha, ",", object@paraPrior$beta,"\n")          
            cat("** Score summary ** \n")                    
            print(summary(score(object)))
            cat("** Fragment length distribution summary ** \n")                    
            print(summary(delta(object)))
            cat("** Summary on the number of binding events per candidate region** \n")
            summary(K(object))
      })


setMethod("summary", "pics",
          function(object)
          {
            cat("** Score ** \n")                    
            cat(score(object),"\n")
            cat("** Fragment length estimate ** \n")                    
            cat(delta(object),"\n")
            cat("** Number of binding events in the candidate region** \n")
            cat(K(object),"\n")
})


setMethod("plot", signature("pics", "segReads"),
function(x, y, addKernel=FALSE, addNucleosome=FALSE, addSe=TRUE, main=NULL, ...)
{
  #Set outer and figure margins to reduce gap between plots
  if(addNucleosome)
  {
    nG<-4
  }
  else
  {
    nG<-2
  }
  par(oma=c(2.5,5,5,5),mar=c(0,5,0,0),cex.lab=2)
  layout(matrix(1:nG,ncol=1), heights = c(.5,.2,.1,.1,.1))

  step<-5
  .densityMix<-function(x,para)
  {
    v<-4
    w<-para$w
    mu<-para$mu
    sigmaSq<-para$sigmaSq
    sigma<-sqrt(sigmaSq)
    xNorm<-outer(-mu,x,"+")/sigma 
    return(colSums(w*dt(xNorm,df=v)/sigma))
  }

  yF<-y@yF
  yR<-y@yR
  cF<-y@cF
  cR<-y@cR        
  map<-y@map
  m<-min(yF[1],yR[1])-100
  M<-max(tail(yF,1),tail(yR,1))+100

  paraR<-list(w=x@estimates$w, mu=x@estimates$mu+x@estimates$delta/2, sigmaSq=x@estimates$sigmaSqR)
  paraF<-list(w=x@estimates$w, mu=x@estimates$mu-x@estimates$delta/2, sigmaSq=x@estimates$sigmaSqF)

  dR<-.densityMix(seq(m,M,step),paraR)
  dF<-.densityMix(seq(m,M,step),paraF)
  maxRange<-max(c(dF,dR))
  plot(seq(m,M,step),dF,xlim=c(m,M),ylim=c(0,maxRange),lty=2,type="l",xlab="",ylab="density",xaxt='n',axes=FALSE)
  title(main=main,outer=TRUE,cex.main=2)
  axis(2)
  axis(1)


  lines(seq(m,M,step),dR,lty=2,col=2)

  # if(length(map)>0)
  # {
    #   nMap<-nrow(map)
    #   for(i in 1:nMap)
    #   {
      #     segments(map[i,1], 0, map[i,2], 0,lwd=3,col=3)
      #   }
      # }

      # Add kernel density estimate
      if((addKernel==TRUE) & (length(yF)>1 & length(yR)>1))
      {
        dkF<-density(yF)
        dkR<-density(yR)
        lines(dkF,lty=3)
        lines(dkR,col=2,lty=3)
      }

      #Add single components and se's
      K<-length(x@estimates$w)
      for(k in 1:K)
      {
        paraR<-list(w=x@estimates$w[k], mu=x@estimates$mu[k]+x@estimates$delta[k]/2, sigmaSq=x@estimates$sigmaSqR[k])
        paraF<-list(w=x@estimates$w[k], mu=x@estimates$mu[k]-x@estimates$delta[k]/2, sigmaSq=x@estimates$sigmaSqF[k])

        dsR<-.densityMix(seq(m,M,step),paraR)
        dsF<-.densityMix(seq(m,M,step),paraF)

        lines(seq(m,M,step),dsF,lty=1)
        lines(seq(m,M,step),dsR,col=2,lty=1)
      }

      stripchart(yF[1],pch=">",method="overplot",cex=2,at=.55,add=FALSE,axes=FALSE,xlim=c(m,M),ylim=c(0,1))        
      if(length(map)>0)
      {
        nMap<-nrow(map)
        symbols((map[,1]+map[,2])/2,rep(.35,nMap),rectangle=cbind(map[,2]-map[,1],rep(.6,nMap)), inches=FALSE, bg=grey(.6), fg=0, add=TRUE,xlim=c(m,M),ylim=c(0,1))
      }

      stripchart(yF,pch=">",method="overplot",cex=2,at=.55,axes=FALSE,xlim=c(m,M),ylim=c(0,1),add=TRUE)
      mtext("IP",cex=1.2,side=2,las=2,at=.5)                  
      stripchart(yR,pch="<",method="overplot",cex=2,at=.45,col=2,add=TRUE)

      abline(h=.35,lty=3)
      if(addSe)
      {
        points(x@estimates$mu,rep(.35,K),pch="+",cex=2)
        if (any(x@estimates$seMu!=0))
        {
          points(x@estimates$mu-2*x@estimates$seMu,rep(.35,K),pch="[",cex=1)
          points(x@estimates$mu+2*x@estimates$seMu,rep(.35,K),pch="]",cex=1)
          segments(x@estimates$mu-2*x@estimates$seMu,rep(.35,K),x@estimates$mu+2*x@estimates$seMu,rep(.35,K),lwd=1,lty=rep(1,K))
        }
      }

      if(length(cF)>0)
      {
        stripchart(cF,pch=">",method="overplot",at=0.25,cex=2,add=TRUE,xlim=c(m,M),ylab="Cont.",axes=FALSE)
      }
      if(length(cR)>0)
      {
        stripchart(cR,pch="<",method="overplot",at=0.15,cex=2,col=2,add=TRUE)
      }
      mtext("Cont.",cex=1.2,side=2,las=2,at=.2)

      if(addNucleosome)
      {
        plot(c(m,M),c(0,1),axes=FALSE,col=0,ylim=c(0,1),xlim=c(m,M),ylab="")
        if(addSe)
        {
          symbols(x@estimates$mu,rep(.5,K),rec=matrix(rep(c(147,.8),K),ncol=2,byrow=TRUE), inches=FALSE, bg=grey(.5*pmin(se(x)/50,1)), fg=0, add=TRUE,xlim=c(m,M),ylim=c(0,1))
        }
        else
        {
          symbols(x@estimates$mu,rep(.5,K),rec=matrix(rep(c(147,.8),K),ncol=2,byrow=TRUE), inches=FALSE, fg=0, bg=1, add=TRUE,xlim=c(m,M),ylim=c(0,1))            
        }
        mtext("Nucl.",cex=1.2,side=2,las=2,at=.5)
      }
})


setMethod("plot", signature("picsError", "segReads"),
          function(x, y, addKernel=FALSE, main=NULL, ...)
      {
        par(oma=c(2.5,5,5,5),mar=c(0,5,0,0))
        layout(matrix(1:2,ncol=1), heights = c(.2,.1))
        
        yF<-y@yF
        yR<-y@yR
        cF<-y@cF
        cR<-y@cR
        map<-y@map
        m<-min(yF[1],yR[1])-100
        M<-max(tail(yF,1),tail(yR,1))+100

        stripchart(yF,pch=">",method="overplot",cex=2,at=.5,add=FALSE,axes=FALSE,xlim=c(m,M),ylab="Cont | Inp.",ylim=c(0,1))
        stripchart(yR,pch="<",method="overplot",cex=2,at=.5,col=2,add=TRUE)

        abline(h=.35,lty=3)

        # Add kernel density estimate
        if((addKernel==TRUE) & (length(yF)>1 & length(yR)>1))
        {
          dkF<-density(yF,bw=75)
          dkR<-density(yR,bw=75)
          plot(dkF,lty=3)
          lines(dkR,col=2,lty=3)
        }

        if(length(cF)>0)
        {
          stripchart(cF,pch=">",method="overplot",at=0.2,cex=2,add=TRUE,xlim=c(m,M),ylab="Cont.",axes=FALSE)
        }
        if(length(cR)>0)
        {
          stripchart(cR,pch="<",method="overplot",at=0.2,cex=2,col=2,add=TRUE)
        }
})

setMethod("plot", signature("picsList", "segReadsList"),
          function(x, y, regionIndex=NULL, addKernel=FALSE, addNucleosome=FALSE, addSe=TRUE,main=NULL, ...)
{
  if(is.null(main))
  {
    setMain<-TRUE
  }
  if(is.null(regionIndex))
  {
    regionIndex<-1:length(x@List)
  }
  for(i in regionIndex)
  {    
    if(setMain)
    {
      main<-paste(as.character(i)," (",y@List[[i]]@chr,")",sep="")
    }
    if(class(x@List[[i]])!="picsError")
    {
      plot(x@List[[i]],y@List[[i]],addKernel=addKernel, addNucleosome=addNucleosome, addSe=addSe,main=main,...)
    }
    else
    {
      plot(x@List[[i]],y@List[[i]],addKernel=addKernel, main=paste(as.character(i)," (",y@List[[i]]@chr,")",sep=""),...)
      warning("Object of class picsError, no PICS density displayed")
    }
  }
})

setMethod("plot", signature("picsList", "picsList"),
          function(x, y, filter=NULL, h=.1, ...)
{
  FDR<-picsFDR(x,y,filter=filter)
  arg<-list(...)
  par(mar=c(4, 4, 4.5, 4) + 0.1)
  # points(FDR[,2],FDR[,3]/max(FDR[,3]),xaxt="n",yaxt="n",lty=3,col=3,pch=2)
  if(length(arg$xlim)!=2)
  {
    xlim<-range(FDR[,2])
  }
  else
  {
    xlim<-c(max(arg$xlim[1],min(FDR[,2])),min(arg$xlim[2],max(FDR[,2])))
  }
  plot(FDR[,2],FDR[,1],xlab="score",ylab="FDR",panel.first=grid(nx=50),...)
  xx<-FDR[FDR[,2]>xlim[1] & FDR[,2]<xlim[2],2]
  yy<-FDR[FDR[,2]>xlim[1] & FDR[,2]<xlim[2],3]
  xx<-xx[seq(1,length(xx),length.out=10)]
  yy<-yy[seq(1,length(yy),length.out=10)]
  axis(3,at=xx,labels=yy)
  mtext("# regions", side = 3, line = 3, ...)
  FDRex<-FDR[FDR[,1]>0,]
  notDup<-rev(!duplicated(rev(FDRex[,1])))
  lines(FDRex[notDup,2],FDRex[notDup,1],col=2,lty=2,lwd=1.5)
  abline(h=h,lw=1.5,col="grey")  
})
