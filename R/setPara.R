## Functions to return a list of parameters to be used by PICS functions
#The default paraneters are set for PICS
setParaEM<-function(minK=1,maxK=15,tol=1e-4,B=100,mSelect="BIC",mergePeaks=TRUE,mapCorrect=TRUE,dataType=NULL)
{
  if(!is.null(dataType))
  {
    if(tolower(dataType)=="mnase" | tolower(dataType)=="sonicated")
	{
		#message("Using the default paraEM for MNase/sonicated data")
		minK=0;maxK=0;tol=1e-4;B=100;mSelect="AIC3";mergePeaks=TRUE;mapCorrect=TRUE;
	}
	#else if(tolower(dataType)=="sonicated")
	#{
	#	message("Using the default paraEM for ChIP-Seq data")
	#	minK=1;maxK=15;tol=1e-4;B=100;mSelect="BIC";mergePeaks=TRUE;mapCorrect=TRUE;
	#}
	else
	{
		stop("Invalid dataType")
	}
  }
#  if(dataType!="TF" & dataType!="H")
#  {
#    stop("Object 'dataType' must be either 'TF' or 'H'")
#  }
  if(!is.finite(tol) & tol<=0 & tol>1)
  {
    stop("'tol' must be a positive number between 0 and 1")
  }
  if(!is.finite(maxK) & maxK<=0 & round(maxK)!=maxK)
  {
    stop("'maxK' must be a positive integer")
  }
  if(!is.finite(B) & B<1 & round(B)!=B)
  {
    stop("'B' must be a positive integer")
  }
  if(!is.character(mSelect) & ((mSelect!="BIC") | (mSelect!="AIC3")))
  {
    stop("'lambda' must be a positive number")
  }
  if(!is.logical(mergePeaks))
  {
    stop("'mergePeaks' must be a logical value")
  }
  if(!is.logical(mapCorrect))
  {
    stop("'mapCorrect' must be a logical value")
  }
  return(list(minK=minK,maxK=maxK,tol=tol,B=B,mSelect=mSelect,mergePeaks=mergePeaks,mapCorrect=mapCorrect))
}

#default for PICS
setParaPrior<-function(xi=200,rho=1,alpha=20,beta=40000,lambda=0,dMu=0, dataType=NULL, PExi=0)
{
  if(!is.null(dataType))
  {
	  if(tolower(dataType)=="mnase")
	  {
		  message("Using the default paraPrior for MNase data, for sonicated data use the argument dataType='sonicated'")
		  xi=150;rho=3;alpha=20;beta=20000;lambda=-0.000064;dMu=200; #Xuekui's seal of approval
	  }
	  else if(tolower(dataType)=="sonicated")
	  {
		  message("Using the default paraPrior for sonicated data")
		  xi=150;rho=1.2;alpha=10;beta=20000;lambda=-0.000064;dMu=200; #From Xuekui's email
	  }
	  else
	  {
		  stop("Invalid dataType")
	  }
  }
  if(!is.finite(xi))
  {
    stop("'xi' must be a numeric value")
  }
  if(!is.finite(rho) & rho<=0)
  {
    stop("'rho' must be a positive number")
  }
  if(!is.finite(alpha) & alpha<=0)
  {
    stop("'alpha' must be a positive number")
  }
  if(!is.finite(beta) & beta<=0)
  {
    stop("'beta' must be a positive number")
  }
  if(!is.finite(lambda) & lambda<=0)
  {
    stop("'lambda' must be a positive number")
  }
  if(!is.finite(dMu) & dMu<=0)
  {
    stop("'dMu' must be a positive number")
  }
  if(PExi>0)
  {
    xi<-PExi
  }
  return(list(xi=xi,rho=rho,alpha=alpha,beta=beta,lambda=lambda,dMu=dMu))
}
