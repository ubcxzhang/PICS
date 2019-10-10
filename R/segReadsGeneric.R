#perform the segmentation depending on the package
segReadsGeneric<-function(data, dataC=NULL, map=NULL, minReads=2, minReadsInRegion=3, jitter=FALSE, maxLregion=0,minLregion=100,
			step=20, width=250, package="PICS")
{
	#maxLregion save max allowable region length, if it is not positive, that means no upper bound the regions length
	
	## Check that we have the right data type
	if(!is(data,"GRanges"))
	{
		stop("The input data should be 'GRanges' object. Provided: ", class(data))
	}
	## Check that we have the same number of chromosomes
	if(!is.null(dataC))
	{
		if(length(levels(seqnames(dataC)))!=length(levels(seqnames(data))))
		{
			stop("Your IP and control data do not have the same number of chromosomes. IP: ", length(levels(seqnames(dataC)))," Control: ",length(levels(seqnames(data))))
		}
	}

	## Total number of reads per sample
	lIP<-length(data)
	
	if(is.null(minReads))
	{
		#Done once per chr
		chrs<-levels(seqnames(data))
		length<-vector('list',length(chrs))
		names(length)<-chrs
		for(cc in chrs)
		{
			length[[cc]]<-diff(c(
							min(start(data[strand(data)=="+"])[1], end(data[strand(data)=="-"])[1]),
							max(tail(start(data[strand(data)=="+"]),1),tail(end(data[strand(data)=="-"]),1))
					))
		}
		minReads<-(ceiling(as.numeric(lIP)/(2*as.numeric(length))*width))
		minReads<-as.integer(names(which.max(table(minReads))))
		minReads<-min(minReads,5)
		minReads<-max(minReads,2)
		print(paste("We automatically calculated minReads, which is ", minReads,".", sep=""))
	}
	
	minReadsInRegion=max(minReadsInRegion,minReads)
	
	paraSW<-list(step=as.integer(step), width=as.integer(width), minReads=as.integer(minReads))
	if(!is.null(map) & !is(map,"GRanges"))
	{
		stop("Map should be a 'GRanges' object. Provided:", class(map))
	}
	else if(is.null(map))
	{
		start<-NULL
		end<-NULL
	}
	else
	{
		map<-map[as.character(seqnames(map)) %in% levels(seqnames(data))]
		chrs<-levels(seqnames(data))
		start<-end<-vector('list',length(chrs))
		names(start) <- names(end) <- chrs
		for(cc in chrs)
		{
			start[[cc]]<-start(map[seqnames(map)==cc])
			end[[cc]]<-end(map[seqnames(map)==cc])
		} 
	}
	
	if (maxLregion>0) maxStep=(maxLregion-2*paraSW$width)/paraSW$step else maxStep=0
	
	## Prepare C input:
	data<-.formatCInput(data)
	if(!is.null(dataC))
	{
		lCont<-length(dataC) #before data transformation
		dataC<-.formatCInput(dataC)
	}
	#If no control, build an empty object of the same size
	else
	{
		dataC<-vector('list',length(data))
		names(dataC)<-names(data)
		for(cc in names(data))
		{
			dataC[[cc]]<-vector('list',2)
			names(data[[cc]])<-c("+","-")
		}
		lCont<-0
	}
	
	## Perform the segmentation
	newSegReadsList<-.Call("segReadsAll", data, dataC, start, end, as.integer(jitter), paraSW , as.integer(maxStep), as.integer(minLregion), pPackage=package, PACKAGE="PICS")
	
	
	newSegReadsList<-unlist(newSegReadsList,recursive=FALSE,use.names=FALSE)
	if(is.null(newSegReadsList))
	{
		stop("No Candidate regions found, you should decrease 'minReads'")
	}
	#newSet<-segReadsList(newSegReadsList,paraSW,as.integer(sum(unlist(lIP))),as.integer(sum(unlist(lCont))))
	newSet<-segReadsList(newSegReadsList,paraSW,as.integer(lIP),as.integer(lCont))
	
	ttt=summarySeg(newSet)
	indrm=((ttt$L<minLregion)|(ttt$NF<minReadsInRegion)|(ttt$NR<minReadsInRegion))
	newSet@List=newSet@List[!indrm]
	
	return(newSet)
}

## summary a segmentList object return some information of each segment as a dataframe
## returned info including 
# chr: chromosome id
# NF : number of forward reads
# NR : number of reverse reads
# L  : length of segment
# min: start location of segments
# max: end location of segments
summarySeg <- function(seg)
{
	temp<-.Call("getSegL", seg@List, PACKAGE="PICS");
	ans <- data.frame(chr=temp[[1]],NF=temp[[2]],NR=temp[[3]],L=temp[[4]],min=temp[[5]],max=temp[[6]])
	ans$chr <- as.character(ans$chr)
	return(ans)
}


## Input: a GRanges object
## Output: a list: list$chr$strand
.formatCInput<-function(GRObject)
{
	chrs<-levels(seqnames(GRObject))
	lData<-vector('list',length(chrs))
	names(lData)<-chrs
	for(cc in chrs)
	{
		GRccObject<-GRObject[seqnames(GRObject)==cc]
		lData[[cc]]<-vector('list',2)
		names(lData[[cc]])<-c("+","-")
		lData[[cc]][["+"]]<-start(GRccObject[strand(GRccObject)=="+"])
		lData[[cc]][["-"]]<-end(GRccObject[strand(GRccObject)=="-"])
	}
	return(lData)
}

###
# bam2gr:
#   INPUT: A bam file with paired-end sequencingdata
#          An optional chr if only selected chr are needed
#   OUTPUT: A GRanges object that can be used in segmentPING
###
bam2gr<-function(bamFile, chr=NULL, PE=FALSE, verbose=FALSE){
  paras <- ScanBamParam(what=c("qname", "rname", "strand", "pos", "mapq", "qwidth"), flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=FALSE))
  bga<-readGAlignments(bamFile, use.names=TRUE, param=paras)
  if(verbose){
    cat(length(bga)," Reads in '",bamFile,"'","\n", sep="")
  }
  hiQScoreIdx<-which(elementMetadata(bga)$mapq>10)
  if(verbose){
    cat(length(bga)-length(hiQScoreIdx)," Reads with low quality scores filtered out","\n")
  }
  bga<-bga[hiQScoreIdx]#filter out elt with low quality scores
  gr<-GRanges()
  
  if(!is.null(chr)){
    chrs<-chr
  } else{
    chrs<-unique(as.character(runValue(seqnames(bga))))
  }
  
  for(i in 1:length(chrs)){   
    cat("Chromosome ", chrs[[i]], "\n")
    bga2<-bga[seqnames(bga)==chrs[i]]
    if(isTRUE(PE)){   
      #change names
      qname<-elementMetadata(bga2)$qname
      qname<-substring(qname,15)
      qname<-gsub(pattern="/3", replacement="", qname)
      elementMetadata(bga2)$qname<-qname
      #merge pairs
      asdf<-as(bga2, "data.frame")
      if(verbose){
        df<-reshape(asdf, timevar="strand", idvar="qname", direction="wide")
      } else{
        suppressWarnings(df<-reshape(asdf, timevar="strand", idvar="qname", direction="wide"))
      }
      df2<-df[,c("start.+", "end.-")]
      rownames(df2)<-df[,"qname"]
      badReads<-which(df2$`start.+`>df2$`end.-`)
      if(length(badReads)>0){
        df2<-df2[-badReads,]
      }
      #Split PE and SE
      idx <- is.na(df2[,c("start.+","end.-")])
      reads <- list(P=df2[!(idx[,1]|idx[,2]),], yFm=df2[idx[,2],], yRm=df2[idx[,1],])
      #Build object
      ir<-IRanges(start=reads$P$`start.+`, end=reads$P$`end.-`)
      gr<-c(gr,GRanges(ranges=ir, seqnames=chrs[i])) #GR with PE (no SE)
    } else{
      ir<-IRanges(start=elementMetadata(bga2)$pos, width=elementMetadata(bga2)$qwidth)
      gr<-c(gr,GRanges(ranges=ir, strand=elementMetadata(bga2)$strand, seqnames=chrs[i]))
    }
  }
  return(gr)
}



