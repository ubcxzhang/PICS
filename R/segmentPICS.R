segmentPICS<-function(data, dataC=NULL, map=NULL, minReads=2, minReadsInRegion=3, jitter=FALSE, dataType="TF",maxLregion=0,minLregion=100)
{
  #Paras depends on the datatype
  step=20
  if(dataType=="TF")  width=250
  if(dataType=="H")   width=150
  
  newSet<-segReadsGeneric(data, dataC=dataC, map=map, minReads=minReads, minReadsInRegion=minReadsInRegion, jitter=jitter, maxLregion=maxLregion,minLregion=minLregion, step=step, width=width, package="PICS")
  return(newSet)
}

###
# bam2gr:
#   INPUT: A bam file with paired-end sequencingdata
#          An optional chr if only selected chr are needed
#   OUTPUT: A GRanges object that can be used in segmentPING
###
bam2gr<-function(bamFile, chr=NULL, PE=FALSE, verbose=FALSE)
{
  #paras <- ScanBamParam(what=c("qname", "rname", "strand", "pos", "mapq", "qwidth"), flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=FALSE))
  paras <- ScanBamParam(what=c("qname", "rname", "strand","mapq"), flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=FALSE))
  #bga<-readGAlignments(bamFile, use.names=TRUE, param=paras)
  bga<-readGAlignments(bamFile, use.names=FALSE, param=paras)
  if(verbose){
    cat(length(bga)," Reads in '",bamFile,"'","\n", sep="")
  }
  hiQScoreIdx<-which(elementMetadata(bga)$mapq>10)
  if(verbose){
    cat(length(bga)-length(hiQScoreIdx)," Reads with low quality scores filtered out","\n")
  }
  bga<-bga[hiQScoreIdx]#filter out elt with low quality scores

  if(isTRUE(PE)){
    qname<-elementMetadata(bga)$qname
    qname<-substring(qname,15)
    qname<-gsub(pattern="/3", replacement="", qname)
    elementMetadata(bga)$qname<-qname
    #merge pairs
    asdf<-as(bga, "data.frame")
    if(verbose){
      df<-reshape(asdf, timevar="strand", idvar="qname", direction="wide")
    } else{
      suppressWarnings(df<-reshape(asdf, timevar="strand", idvar="qname", direction="wide"))
    }   
    df2<-df[,c("start.+", "end.-","rname.+")]
    colnames(df2)<-c("start", "end", "chr")
    rownames(df2)<-df[,"qname"]
    badReads<-which(df2$start>df2$end)
    if(length(badReads)>0){
      df2<-df2[-badReads,]
    }   
    #Split PE and SE
    idx <- is.na(df2[,c("start","end")])
    reads <- df2[!(idx[,1]|idx[,2]),]
    yFm <- df2[idx[,2],] #Forward reads, missing matching reverse
    yRm <- df2[idx[,1],] #Reverse reads
    gr<-GRanges(ranges=IRanges(start=reads$start, end=reads$end), strand="*", seqnames=reads$chr)
  } else{
    gr<-GRanges(ranges=IRanges(start=start(bga), end=end(bga)), strand=strand(bga), seqnames=seqnames(bga))
  }
  return(gr)
}

