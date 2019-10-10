###################################################
###  Identify Candidate region 
###################################################
## PE.RD is a file including perfect sequenced pairs and one-side sequenced pairs ##
## This data consists of start, end position and strand ##
candidate.region <- function(PE.RD, islandDepth, min_cut, max_cut){
	min_region <- min(start(PE.RD))
    	max_region <- max(end(PE.RD))
	
## get the coverage for each base
	cvg <- coverage(PE.RD)#, shift=-min_region+1, width=(max_region-min_region+1)) ##R# why width?
## keep only the covered regions
	read.region <- slice(cvg, lower=1)
## get the highest covered regions
	max.depth <- viewMaxs(read.region)
## Identify candidate regions having at least islandDepth height ##
	cand.region <- read.region[max.depth>=islandDepth]
##  cand.region.temp <- peakSummary(cand.region)
	cand_IR <- IRanges(start=start(cand.region), end=end(cand.region))

## get index of regions with width > max_cut    
	max_which <- which(width(cand_IR)>max_cut)
	if(length(max_which)>0){
		over_cand_IR <- cand_IR[max_which,]
		# cand_recu <- NULL
		cand_recu <- IRanges(NULL)
#		overlaps <-  findOverlaps(PE.RD,over_cand_IR)@matchMatrix
		overlaps<- as.matrix(findOverlaps(PE.RD,over_cand_IR))
		idx <- split(overlaps[,1], overlaps[,2])
		for(j in 1:length(over_cand_IR)){
			# PE.RD.Cand <- PE.RD[(PE.RD %in% over_cand_IR[j,]) ==TRUE,] # %in% stop working after PING is loaded
			# PE.RD.Cand <- PE.RD[(end(PE.RD) >= start(over_cand_IR)[j]) & (start(PE.RD) <= end(over_cand_IR)[j]), ]
			PE.RD.Cand <- PE.RD[idx[[j]],]
			temp_cand_recu <- cand_recursive(PE.RD.Cand, over_cand_IR[j,], islandDepth, min_region, max_region, min_cut, max_cut)
			# cand_recu <- rbind(cand_recu, temp_cand_recu)
			cand_recu <- c(cand_recu, temp_cand_recu)
		}
		# cand_recu <- IRanges(start=cand_recu$start, end=cand_recu$end)
		cand_IR <-c(cand_IR[-max_which, ], cand_recu)
	}
	return(cand_IR)
}


####################################################################################
###  Cut the wide ranged candidate region into max_cut ranged candidate regions 
## This function is for each candidate region having wide range 
####################################################################################
cand_recursive <- function(PE.RD.Cand, over_cand_IR, islandDepth, min_region, max_region, min_cut, max_cut) 
{
	
	temp_over_IR <- over_cand_IR
	width_temp <- max(width(temp_over_IR))
	islandDepth <- islandDepth + 1
	st <- start(temp_over_IR)
	ed <- end(temp_over_IR)
	while(width_temp>max_cut){
		temp_cand <- s_cut(PE.RD.Cand, temp_over_IR, islandDepth, min_region=min_region, max_region=max_region)
		temp_cand_IR<-temp_cand$cand_IR_sep

		st <- c(st, start(temp_cand_IR))
		ed <- c(ed, end(temp_cand_IR))
				
		temp_over_IR <- temp_cand_IR[width(temp_cand_IR)>max_cut,]
		if(length(temp_over_IR)>0){
			width_temp <- max(width(temp_over_IR))	
		}else{
			width_temp <- max(width(temp_cand_IR))
		}

		islandDepth <- islandDepth + 1
		#cat("islandDepth=", islandDepth)
		
	}	
	cand_IR_sep <- IRanges(start=sort(st), end=sort(ed))
	each_cand_IR <- disjoin(cand_IR_sep)
	
	cand_IR <- each_cand_IR[width(each_cand_IR)>min_cut]
# 	browser()
# 	return(as.data.frame(cand_IR))
}

s_cut <- function(PE.RD.Cand, over_cand_IR, islandDepth, min_region, max_region)
{
	# PE_reads <- PE.RD.Cand[(PE.RD.Cand %in% over_cand_IR) ==TRUE,]  # %in% stop working after PING is loaded
#	overlaps <- findOverlaps(PE.RD.Cand, over_cand_IR)@matchMatrix
	overlaps <- as.matrix(findOverlaps(PE.RD.Cand, over_cand_IR))
	PE_reads <- PE.RD.Cand[overlaps[,1],]
	
	cvg_temp <- coverage(PE_reads)##, width=(max_region-min_region+1))

	s_temp <- slice(cvg_temp, lower=islandDepth)	
	st <- start(s_temp)
	ed <- end(s_temp)
	
	cand_IR_sep <- IRanges(start=sort(st), end=sort(ed))

	return(list(cand_IR_sep=cand_IR_sep, minS=unique(viewMins(s_temp))))
}


#############################################################
###  Identify reads included in each candidate regions 
#############################################################
segChrRead <- function(candidate_RD, PE.RD, PEMF.RD, PEMR.RD , PEC.RD=NULL, PECMF.RD=NULL, PECMR.RD=NULL, map.Start, map.End, chr)
{
	seg <- vector("list", length(candidate_RD))
	map.IR <- IRanges(start=map.Start, end=map.End)
 	index_map_all <-  as.matrix(findOverlaps(candidate_RD, map.IR))
#	index_map_all <-  findOverlaps(candidate_RD, map.IR)@matchMatrix
	index_map_cand <- index_map_all[,1]
	index_map_map <- index_map_all[,2]

	flag <- rep("FALSE", length(candidate_RD))
	flag[index_map_cand] <- "TRUE"
	for (i in 1: length(candidate_RD)){
		target <- candidate_RD[i,]
		
		## Pull out which rows are included in the candidate region ##
 		index_PE <- as.matrix(findOverlaps(target, PE.RD))[,2]
 		index_PEMF <- as.matrix(findOverlaps(target, PEMF.RD))[,2]
 		index_PEMR <- as.matrix(findOverlaps(target, PEMR.RD))[,2]
#		index_PE <- findOverlaps(target, PE.RD)@matchMatrix[,2]
#		index_PEMF <- findOverlaps(target, PEMF.RD)@matchMatrix[,2]
#		index_PEMR <- findOverlaps(target, PEMR.RD)@matchMatrix[,2]

		yF <- start(PE.RD)[index_PE]
		yR <- end(PE.RD)[index_PE]
		
		# I modified this part
		ord <- order(yF,yR)
		yF <- yF[ord]
		yR <- yR[ord]
		
		if(length(index_PEMF)>0){
			yFm <- start(PEMF.RD)[index_PEMF]
		}else{
			yFm <- numeric(0)
		}
		if(length(index_PEMR)>0){
			yRm <- start(PEMR.RD)[index_PEMR]
		}else{
			yRm <- numeric(0)
		}
				
		if(!is.null(PEC.RD)){
 			index_PEC <- as.matrix(findOverlaps(target, PEC.RD))[,2]
 			index_PECMF <- as.matrix(findOverlaps(target, PECMF.RD))[,2]
 			index_PECMR <- as.matrix(findOverlaps(target, PECMR.RD))[,2]
#			index_PEC <-   findOverlaps(target, PEC.RD )@matchMatrix[,2]
#			index_PECMF <- findOverlaps(target, PECM.RD)@matchMatrix[,2]
#			index_PECMR <- findOverlaps(target, PECM.RD)@matchMatrix[,2]
			
			cF <- start(PEC.RD)[index_PEC]
			cR <- end(PEC.RD)[index_PEC]

			if(length(index_PECMF)>0){
				cFm <- start(PECMF.RD)[index_PECMF]
			}else{
				cFm <- numeric(0)
			}
			if(length(index_PECMR)>0){
				cRm <- start(PECMR.RD)[index_PECMR]	
			}else{
				cRm <- numeric(0)
			}
		}

		## resize the mappability region ##
		if (flag[i]==TRUE){
			 index_map_within <-  as.matrix(findOverlaps(target, map.IR[index_map_map[index_map_cand==i],], type="within"))[,2]
#			index_map_within <-  findOverlaps(target, map.IR[index_map_map[index_map_cand==i],], type="within")@matchMatrix[,2]
			if(length(index_map_within)>0){
				map <- cbind(start(map.IR[index_map_map[index_map_cand==i],]), end(map.IR[index_map_map[index_map_cand==i],]))
			}else{
				index_map <-  index_map_map[index_map_cand==i]
				if(length(index_map_within)>0){
					index_map_out <- index_map[-index_map_within]
				}else{
					index_map_out <- index_map
				}
				if(length(index_map_out)>0){
					temp_map_IR <- map.IR[index_map_out,]
		
					for (k in 1:length(index_map_out)){
						if(start(temp_map_IR[k,]) < start(target)){
							start(temp_map_IR[k,]) <- start(target)
						}
						if(end(temp_map_IR[k,]) > end(target)){
							end(temp_map_IR[k,]) <- end(target)
						}
					}
					map <- cbind(start(temp_map_IR), end(temp_map_IR))
				}
			}
		}else{
			map <- matrix(0, 0,2)
		}
		seg[[i]] <- segReadsPE(as.numeric(yF), as.numeric(yR), yFm=numeric(0), yRm=numeric(0), cF=numeric(0), cR=numeric(0), cFm=numeric(0), cRm=numeric(0), map=map, chr)
	}	
	return(seg)
}
