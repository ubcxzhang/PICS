### R code from vignette source 'PICS.Rnw'

###################################################
### code chunk number 1: Loading PICS
###################################################
library(PICS)


###################################################
### code chunk number 2: Reading the data
###################################################
path<-system.file("extdata",package="PICS")
dataIP<-read.table(file.path(path,"Treatment_tags_chr21_sort.bed"),header=TRUE,
colClasses=c("factor","integer","integer","factor"))
dataIP<-as(dataIP,"GRanges")

dataCont<-read.table(file.path(path,"Input_tags_chr21_sort.bed"),header=TRUE
,colClasses=c("factor","integer","integer","factor"))
dataCont<-as(dataCont,"GRanges")


###################################################
### code chunk number 3: Reading the mappability profile
###################################################
map<-read.table(file.path(path,"mapProfileShort"),header=TRUE,
colClasses=c("factor","integer","integer","NULL"))
map<-as(map,"GRanges")


###################################################
### code chunk number 4: Genome segmentation
###################################################
seg<-segmentPICS(dataIP, dataC=dataCont, map=map, minReads=1)


###################################################
### code chunk number 5: Cluster initialization
###################################################
library(parallel)


###################################################
### code chunk number 6: PICS analysis
###################################################
pics<-PICS(seg, nCores=2)


###################################################
### code chunk number 7: FDR estimation
###################################################
segC<-segmentPICS(dataCont,dataC=dataIP,map=map,minReads=1)
picsC<-PICS(segC)
fdr<-picsFDR(pics,picsC,filter=list(delta=c(50,Inf),se=c(0,50),
sigmaSqF=c(0,22500),sigmaSqR=c(0,22500)))


###################################################
### code chunk number 8: plot-FDR1
###################################################
plot(pics,picsC,xlim=c(2,8),ylim=c(0,.2),filter=list(delta=c(50,300),
se=c(0,50),sigmaSqF=c(0,22500),sigmaSqR=c(0,22500)),type="l")


###################################################
### code chunk number 9: plot-FDR2
###################################################
plot(fdr[,c(3,1)])


###################################################
### code chunk number 10: Create RangedData data object of enriched regions
###################################################
myFilter=list(delta=c(50,300),se=c(0,50),sigmaSqF=c(0,22500),sigmaSqR=c(0,22500))
rdBed<-makeRangedDataOutput(pics,type="bed",filter=c(myFilter,list(score=c(1,Inf))))


###################################################
### code chunk number 11: Create the bed file (eval = FALSE)
###################################################
## library(rtracklayer)
## export(rdBed,"myfile.bed")


###################################################
### code chunk number 12: Create the wig file (eval = FALSE)
###################################################
## rdBed<-makeRangedDataOutput(pics,type="wig",filter=c(myFilter,list(score=c(1,Inf))))
## export(rdBed,"myfile.wig")


