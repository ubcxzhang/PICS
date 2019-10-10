picsFDR<-function(picsIP,picsCont,filter=list(delta=c(0,Inf),se=c(0,Inf),sigmaSqF=c(0,Inf),sigmaSqR=c(0,Inf)))
{
  score1<-score(picsIP)
  delta1<-delta(picsIP)
  se1<-se(picsIP)
  sf1<-sigmaSqF(picsIP)
  sr1<-sigmaSqR(picsIP)

  score2<-score(picsCont)
  delta2<-delta(picsCont)
  se2<-se(picsCont)
  sf2<-sigmaSqF(picsCont)
  sr2<-sigmaSqR(picsCont)

  # Normalizing ratio to account for the fact that the control and treatment samples might not have
  # the same number of reads
  # ratio<-picsIP@Nc/picsIP@N
  # I now standardized the fdr when I compute scores.
  ratio<-1

  if(!is.null(filter))
  {
    ### Filter based on delta
    ind1.1<-delta1>filter$delta[1] & delta1<filter$delta[2]
    ind1.2<-sf1>filter$sigmaSqF[1] & sf1<filter$sigmaSqF[2]
    ind1.3<-sr1>filter$sigmaSqR[1] & sr1<filter$sigmaSqR[2]
    ind1.4<-se1>filter$se[1] & se1<filter$se[2]

    ind1<-ind1.1&ind1.2&ind1.3&ind1.4

    ind2.1<-delta2>filter$delta[1] & delta2<filter$delta[2]
    ind2.2<-sf2>filter$sigmaSqF[1] & sf2<filter$sigmaSqF[2]
    ind2.3<-sr2>filter$sigmaSqR[1] & sr2<filter$sigmaSqR[2]
    ind2.4<-se2>filter$se[1] & se2<filter$se[2]

    ind2<-ind2.1&ind2.2&ind2.3&ind2.4
  }
  else
  {
    ind1<-rep(TRUE,length(score1))
    ind2<-rep(TRUE,length(score2))
  }
  score1<-score1[ind1]
  score2<-score2[ind2]
  scoreVector<-seq(min(score2,ratio*score1,na.rm=TRUE),max(ratio*score1,na.rm=TRUE),.1)
  FDR<-rep(0,length(scoreVector))
  N<-rep(0,length(scoreVector))
  FDR[1]<-1
  N[1]<-length(score1)

  n<-length(FDR)
  i<-2
  
  while((i<n) & (FDR[i-1]>0))
  {
    # I make the FDR monotonic
    # FDR[i]<-min(1,max(FDR[i+1],sum(score2>scoreVector[i],na.rm=TRUE)/sum(ratio*score1>scoreVector[i],na.rm=TRUE)))
    FDR[i]<-min(1,min(FDR[i-1],sum(score2>scoreVector[i],na.rm=TRUE)/sum(score1>scoreVector[i],na.rm=TRUE)))    
    N[i]<-sum(ratio*score1>scoreVector[i],na.rm=TRUE)
    i<-i+1
  }  
  data.frame(FDR=FDR[1:i],score=scoreVector[1:i],N=N[1:i])
}#end of function FDRs
