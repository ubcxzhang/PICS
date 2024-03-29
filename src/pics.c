#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h>
#include <Rdefines.h>
#include <stdbool.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <float.h>
//GSL stuff
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>

double fun0(double x);
double fun1(double x, int nu);
double fun2(double x);
void ECM1(int nu, SEXP R, SEXP F, SEXP para, double xi, double alpha, double betap, double rho, SEXP a, SEXP b, double tol, double cst, double lambda, double dMu);
void ECM2(int nu, SEXP R, SEXP F, SEXP para, double xi, double alpha, double betap, double rho, SEXP a, SEXP b, double tol, double cst, double lambda, double dMu);
SEXP iterEM(SEXP iMax, SEXP nu, SEXP yR, SEXP yF, SEXP para, SEXP xi, SEXP alpha, SEXP betap, SEXP rho, SEXP a, SEXP b, SEXP tol, SEXP cst, SEXP lambda, SEXP dMu);
double logDensityMix(double *Y, double *w, double *mu, double *sigmaSq, int K, int N);
double plogp(double y, void * params);
double intplogp (int start, int end, double *w, double *mu, double *sigmaSq, int K);
SEXP BIC(SEXP nuC, SEXP R, SEXP F, SEXP para, SEXP dMuC, SEXP lambdaC, SEXP a, SEXP b, SEXP mselect);
SEXP initPara(SEXP F, SEXP R, SEXP kk);
SEXP fitModel(SEXP kk, SEXP iMax, SEXP tol, SEXP mselect, SEXP yR, SEXP yF, SEXP a, SEXP b, SEXP xi, SEXP alpha, SEXP betap, SEXP rho, SEXP lambda, SEXP dMu, SEXP cst, SEXP nu, SEXP minReadPerPeak);
SEXP fitModelK(SEXP kk, SEXP iMax, SEXP tol, SEXP mselect, SEXP yR, SEXP yF, SEXP a, SEXP b, SEXP xi, SEXP alpha, SEXP betap, SEXP rho, SEXP lambda, SEXP dMu, SEXP cst, SEXP nu, SEXP minReadPerPeak);
void printPara(SEXP para);
int mergePeak(SEXP para, gsl_matrix* infMat, gsl_vector* se, gsl_vector* seF, gsl_vector* seR, int *K, double nu, double nSe, double minSpacingPeaks, int dataType);
int getInfMat(SEXP R,SEXP F, SEXP para, SEXP a, SEXP b, double rho, double xi, double alpha, double cst, double lambda, double nu, gsl_matrix* infMat, gsl_vector* se, gsl_vector* seF, gsl_vector* seR);
SEXP fitModelAllk(SEXP segReads, SEXP paraEM, SEXP paraPrior, SEXP minReads, SEXP N, SEXP Nc, SEXP chr);
SEXP fitPICS(SEXP segReadsList, SEXP paraEM, SEXP paraPrior, SEXP minReads);


double fun0(double x)
{
	double y = pow(x,5)/5-pow(x,3)*2/3+x;
	return(y);
}

double fun1(double x, int nu)
{
	double y = pow(1+x*x/nu,-(nu+1.0)/2.0);
	return(y);
}

double fun2(double x)
{
	double 	y=pow(x,3)/3-pow(x,5)/5;
	return(y);
}

  //Main PICS function
SEXP fitPICS(SEXP segReadsList, SEXP paraEM, SEXP paraPrior, SEXP minReads)
{
  int l=0;
  SEXP picsList, segReads, List;
  SEXP N, Nc, chr;
  int nProtected=0;
  
  PROTECT(List=GET_SLOT(segReadsList,install("List"))); nProtected++;
  PROTECT(picsList=NEW_LIST(length(List))); nProtected++;
  
  PROTECT(N=GET_SLOT(segReadsList,install("N"))); nProtected++;
  PROTECT(Nc=GET_SLOT(segReadsList,install("Nc"))); nProtected++;
  
  for(l=0;l<length(List);l++)
  {
    // I have added the option to interrupt R
    R_CheckUserInterrupt();
    segReads=VECTOR_ELT(List,l);
    chr=GET_SLOT(segReads,install("chr"));    
    SET_VECTOR_ELT(picsList,l,fitModelAllk(segReads, paraEM, paraPrior, minReads, N, Nc, chr));
  }
  UNPROTECT(nProtected);
  return(picsList);
}

SEXP fitModelAllk(SEXP segReads, SEXP paraEM, SEXP paraPrior, SEXP minReads, SEXP N, SEXP Nc, SEXP chr)
{
  int Nnucle=0, minK=0, maxKK=0, i=0, k=0, kBest=0, K=0;
  SEXP yF, yR, cF, cR, map, kk, classDef;
  SEXP score, scoreF, scoreR;
  double range=0;
  char myError[]="";
  int nProtected=0;
  SEXP a, b, ans, myPics;
  SEXP Para, estimates, converge, Nmerged, names, name, wOut,muOut,sfOut,srOut,deltaOut,seOut, seOutF,seOutR, rangeOut;
  SEXP mapCorrect, mergePeaks;
  SEXP cst, nu;
  double calpha=0, sumF=0, sumR=0, sumcF=0, sumcR=0;
  double *w, *mu, *delta, *sF, *sR;
  int flag=0;
  gsl_matrix *infMat;
  gsl_vector *se, *seF, *seR;
  SEXP minReadPerPeak;
  SEXP mselect, code;
  
  
  /** Create shortcuts to the respective slots **/
  yF=GET_SLOT(segReads,install("yF"));
  yR=GET_SLOT(segReads,install("yR"));
  cF=GET_SLOT(segReads,install("cF"));
  cR=GET_SLOT(segReads,install("cR"));
  map=GET_SLOT(segReads,install("map"));
  
    //extract the corresponding variables from paraEM
  mapCorrect=VECTOR_ELT(paraEM,6);
  mergePeaks=VECTOR_ELT(paraEM,5);
  
  /** constant **/
  PROTECT(cst=NEW_NUMERIC(1));nProtected++;
  REAL(cst)[0]=gsl_sf_gamma(3.5)/gsl_sf_gamma(3.0)/M_SQRTPI;
  PROTECT(nu=NEW_INTEGER(1));nProtected++;
  INTEGER(nu)[0]=4;
  PROTECT(minReadPerPeak=VECTOR_ELT(minReads,0));nProtected++;
  
  /** range of the data used to infer the min/max number of components **/
  range=(fmax2(REAL(yF)[length(yF)-1],REAL(yR)[length(yR)-1])-fmin2(REAL(yF)[0],REAL(yR)[0]));
  if(REAL(VECTOR_ELT(paraPrior, 5))[0]>0) // automatically decide minK and maxKK for histone data
  {
    Nnucle=(int)(range/200.); //number of nuclesomes in this region
    minK=imax2(1,(int)(Nnucle/3.)); //assume max 66% of region are NFR, Slide-window size =500 < d_mu * 3
    maxKK=imax2((int)(Nnucle*1.5+1),minK+1);
    // 
    // minK=1;
    // maxKK=3;
    
  }
  else
  {
    minK=(int)(REAL(VECTOR_ELT(paraEM, 0))[0]);
    maxKK=(int)fmin2(REAL(VECTOR_ELT(paraEM, 1))[0],range/(REAL(VECTOR_ELT(paraPrior, 0))[0]+4.*sqrt(REAL(VECTOR_ELT(paraPrior, 3))[0]/REAL(VECTOR_ELT(paraPrior, 2))[0]))+1.0);
      //Here we make sure that maxKK>minK
    maxKK=imax2(minK,maxKK);
  }
  
  /** List of the number of components **/
  PROTECT(kk=NEW_INTEGER(maxKK-minK+1));nProtected++;
  
  for(i=0;i<(maxKK-minK+1);i++)
  {
    INTEGER(kk)[i]=minK+i;
  }
  
  /** Initialize the mappability intervals **/
    //Only if we have set mapCorrect to TRUE and we are given a mappability profile
  PROTECT(a=R_NilValue);nProtected++;
  PROTECT(b=R_NilValue);nProtected++;
  
  if(LOGICAL(mapCorrect)[0] & (length(map)!=0))
  {
    PROTECT(a=NEW_INTEGER(length(map)/2));nProtected++;
    PROTECT(b=NEW_INTEGER(length(map)/2));nProtected++;
    for(i=0;i<(length(map)/2);i++)
    {
        // Check that I am getting the elements in the correct order
      INTEGER(a)[i]=INTEGER(map)[i];
      INTEGER(b)[i]=INTEGER(map)[(length(map)/2)+i];
    }
  }
  
  PROTECT(mselect=NEW_INTEGER(1)); nProtected++;
  if(strcmp(CHAR(STRING_PTR(VECTOR_ELT(paraEM, 4))[0]),"BIC")==0)
  {
    INTEGER(mselect)[0]=1;
  }
  else if(strcmp(CHAR(STRING_PTR(VECTOR_ELT(paraEM, 4))[0]),"AIC")==0)
  {
    INTEGER(mselect)[0]=2;
  }
  else if(strcmp(CHAR(STRING_PTR(VECTOR_ELT(paraEM, 4))[0]),"AIC3")==0)
  {
    INTEGER(mselect)[0]=3;
  }
  
  PROTECT(ans=fitModelK(kk, VECTOR_ELT(paraEM, 3), VECTOR_ELT(paraEM, 2), mselect, yR, yF, a, b, VECTOR_ELT(paraPrior, 0), VECTOR_ELT(paraPrior, 2), VECTOR_ELT(paraPrior, 3), VECTOR_ELT(paraPrior, 1), VECTOR_ELT(paraPrior, 4), VECTOR_ELT(paraPrior, 5), cst, nu, minReadPerPeak));
  nProtected++;

  /** Check if we have had an error, if yes we return an object of class picsError **/
  if(strcmp(CHAR(STRING_PTR(VECTOR_ELT(ans,3))[0]),myError) != 0)
  {
    classDef=MAKE_CLASS("picsError");
    PROTECT(myPics=NEW_OBJECT(classDef));
    nProtected++;
    SET_SLOT(myPics,mkChar("errorCode"),VECTOR_ELT(ans,3));
    UNPROTECT(nProtected);
    return(myPics);
  }
  Para=VECTOR_ELT(ans,0);
    
    // best BIC (length of w)
  kBest=length(VECTOR_ELT(Para,0));
  K=kBest;
    //allocate the memory for the information matrix
    //Initialized to zero too
  infMat=gsl_matrix_calloc(5*kBest-1,5*kBest-1);
    //initialize all se's to zero
  se=gsl_vector_calloc(kBest);  
  seF=gsl_vector_calloc(kBest);
  seR=gsl_vector_calloc(kBest);
  
    // Compute teh information matrix
  flag=getInfMat(yR, yF, Para, a, b, REAL(VECTOR_ELT(paraPrior, 1))[0], REAL(VECTOR_ELT(paraPrior, 0))[0], REAL(VECTOR_ELT(paraPrior, 2))[0], REAL(cst)[0], REAL(VECTOR_ELT(paraPrior, 4))[0], INTEGER(nu)[0], infMat, se, seF, seR);
  
  //  //We have an error and dMu==0
  //if((flag!=0) & (REAL(VECTOR_ELT(paraPrior, 5))[0]==0))
  //for histone, there are not many (3 out of 0.2million regions in whole genome) singular problem, and I looked at the regions with that problem, which give us nothing more than noise

  //We have an error for either TF or histone
  if((flag!=0))	
  {
    classDef=MAKE_CLASS("picsError");
    PROTECT(myPics=NEW_OBJECT(classDef));nProtected++;
    PROTECT(code=NEW_CHARACTER(1));nProtected++;
    STRING_PTR(code)[0]=mkChar("Singular information matrix");
    SET_SLOT(myPics,mkChar("errorCode"),code);
    UNPROTECT(nProtected);
    return(myPics);
  }
  else
  {
      // If we could not invert the information matrix, I reset it to zero so that all se's will be zero
      // Another possibility would be to stabilize the information matrix to make it invertible
    gsl_matrix_set_zero(infMat);
  }
  
  
  /** Merge Peaks **/
    // We want to merge peaks
  if(LOGICAL(mergePeaks)[0])
  {
    /** If dMu>0 **/
    if(REAL(VECTOR_ELT(paraPrior, 5))[0]>0) //Histone case
    {
        // Xuekui will take care of this
      mergePeak(Para, infMat, se, seF, seR, &K, INTEGER(nu)[0], 3, 120, 2);
    }
    else //TF case
    {
        mergePeak(Para, infMat, se, seF, seR, &K, INTEGER(nu)[0], 3, 0, 1);
    }
    //We encountered an error when merging
    if(flag!=0)
    {
      classDef=MAKE_CLASS("picsError");
      PROTECT(myPics=NEW_OBJECT(classDef));
      nProtected++;
      PROTECT(code=NEW_CHARACTER(1));
      nProtected++;
      STRING_PTR(code)[0]=mkChar("Singular information matrix");
      SET_SLOT(myPics,mkChar("errorCode"),code);
      UNPROTECT(nProtected);
      return(myPics);
    }
  }
  
  /** Enrichment score **/
  PROTECT(score=NEW_NUMERIC(K));nProtected++;
  PROTECT(scoreF=NEW_NUMERIC(K));nProtected++;
  PROTECT(scoreR=NEW_NUMERIC(K));nProtected++;
  
  calpha=1.5;//percentile cutoff
  
  /** Shortcut to the parameters **/
  w=REAL(VECTOR_ELT(Para,0));
  mu=REAL(VECTOR_ELT(Para,1));
  delta=REAL(VECTOR_ELT(Para,2));
  sF=REAL(VECTOR_ELT(Para,3));
  sR=REAL(VECTOR_ELT(Para,4));
  
  
  for(k=0;k<K;k++)
  {
    sumF=0;
    for(i=0;i<length(yF);i++)
    {
      if(abs(REAL(yF)[i]-mu[k]+delta[k]/2.)/sqrt(sF[k])<calpha)
      {
        sumF++;
      }
    }
    sumR=0;
    for(i=0;i<length(yR);i++)
    {
      if(abs(REAL(yR)[i]-mu[k]-delta[k]/2.)/sqrt(sR[k])<calpha)
      {
        sumR++;
      }
    }

    /** If we have control data, we can normalize the scores **/
    // INTEGER(Nc)[0]==0;
    if(INTEGER(Nc)[0]>0)
    { 
      sumcF=0;
      for(i=0;i<length(cF);i++)
      {
        if(abs(REAL(cF)[i]-mu[k]+delta[k]/2.)/sqrt(sF[k])<calpha)
        {
          sumcF++;
        }
      }
      sumcR=0;
      for(i=0;i<length(cR);i++)
      {
        if(abs(REAL(cR)[i]-mu[k]-delta[k]/2.)/sqrt(sR[k])<calpha)
        {
          sumcR++;
        }
      }
      REAL(scoreF)[k]=sumF/(INTEGER(N)[0]*fmax2(sumcF,2.)/INTEGER(Nc)[0]);
      REAL(scoreR)[k]=sumR/(INTEGER(N)[0]*fmax2(sumcR,2.)/INTEGER(Nc)[0]);
    }
    else
    {
      REAL(scoreF)[k]=sumF;
      REAL(scoreR)[k]=sumR;
    }
      // Sum of the scores
    REAL(score)[k]=REAL(scoreF)[k]+REAL(scoreR)[k];
  }
  
  classDef=MAKE_CLASS("pics");
  PROTECT(myPics=NEW_OBJECT(classDef));nProtected++;
  PROTECT(estimates=NEW_LIST(8));nProtected++;
  
  PROTECT(wOut=NEW_NUMERIC(K));nProtected++;
  PROTECT(muOut=NEW_NUMERIC(K));nProtected++;
  PROTECT(deltaOut=NEW_NUMERIC(K));nProtected++;
  PROTECT(sfOut=NEW_NUMERIC(K));nProtected++;
  PROTECT(srOut=NEW_NUMERIC(K));nProtected++;
  PROTECT(seOut=NEW_NUMERIC(K));nProtected++;
  PROTECT(seOutF=NEW_NUMERIC(K));nProtected++;
  PROTECT(seOutR=NEW_NUMERIC(K));nProtected++;
  
  
  for(k=0;k<K;k++)
  {
    REAL(wOut)[k]=w[k];
    REAL(muOut)[k]=mu[k];
    REAL(deltaOut)[k]=delta[k];
    REAL(sfOut)[k]=sF[k];
    REAL(srOut)[k]=sR[k];
    REAL(seOut)[k]=gsl_vector_get(se,k);
    REAL(seOutF)[k]=gsl_vector_get(seF,k);
    REAL(seOutR)[k]=gsl_vector_get(seR,k);
  }
  
  SET_VECTOR_ELT(estimates,0,wOut); //w
  SET_VECTOR_ELT(estimates,1,muOut); //mu
  SET_VECTOR_ELT(estimates,2,deltaOut); //delta
  SET_VECTOR_ELT(estimates,3,sfOut); //sF^2
  SET_VECTOR_ELT(estimates,4,srOut); //sR^2
  
  SET_VECTOR_ELT(estimates,5,seOut); //seMu
  SET_VECTOR_ELT(estimates,6,seOutF); //seMuF
  SET_VECTOR_ELT(estimates,7,seOutR); //seMuR
  
  PROTECT(names = allocVector(STRSXP, 8)); nProtected++;
  SET_STRING_ELT(names, 0, mkChar("w"));
  SET_STRING_ELT(names, 1, mkChar("mu"));
  SET_STRING_ELT(names, 2, mkChar("delta"));
  SET_STRING_ELT(names, 3, mkChar("sigmaSqF"));
  SET_STRING_ELT(names, 4, mkChar("sigmaSqR"));
  SET_STRING_ELT(names, 5, mkChar("seMu"));
  SET_STRING_ELT(names, 6, mkChar("seMuF"));
  SET_STRING_ELT(names, 7, mkChar("seMuR"));
  setAttrib(estimates, R_NamesSymbol, names);  
  
  SET_SLOT(myPics,mkChar("estimates"),estimates);
  SET_SLOT(myPics,mkChar("score"),score);
  SET_SLOT(myPics,mkChar("scoreF"),scoreF);
  SET_SLOT(myPics,mkChar("scoreR"),scoreR);
  
  
    //set the name of the estimates list
  
  PROTECT(Nmerged=NEW_NUMERIC(1));nProtected++;
  REAL(Nmerged)[0]=kBest-K;
  SET_SLOT(myPics,mkChar("Nmerged"),Nmerged);
  
  
  PROTECT(converge=allocVector(LGLSXP,1));nProtected++;
    //Extract the logical from the best answer
  LOGICAL(converge)[0]=LOGICAL(VECTOR_ELT(ans,2))[0];
  SET_SLOT(myPics,mkChar("converge"),converge);
  
    //Range of the data
  PROTECT(rangeOut=NEW_NUMERIC(2));nProtected++;
  REAL(rangeOut)[0]=fmin2(REAL(yF)[0],REAL(yR)[0]);
  REAL(rangeOut)[1]=fmax2(REAL(yF)[length(yF)-1],REAL(yR)[length(yR)-1]);
  SET_SLOT(myPics,mkChar("range"),rangeOut);
  
  PROTECT(name=NEW_CHARACTER(1));nProtected++;
  
    // Set the chromosome name
  SET_STRING_ELT(name,0,STRING_PTR(chr)[0]);    
  SET_SLOT(myPics,mkChar("chr"),name);
  
  UNPROTECT(nProtected);
  return(myPics);
}



SEXP fitModelK(SEXP kk, SEXP iMax, SEXP tol, SEXP mselect, SEXP yR, SEXP yF, SEXP a, SEXP b, SEXP xi, SEXP alpha, SEXP betap, SEXP rho, SEXP lambda, SEXP dMu, SEXP cst, SEXP nu, SEXP minReadPerPeak)
{
	int  nProtected=0, kmax=length(kk), k;
	double bestBIC=GSL_NEGINF;
	SEXP ans, temp, firstFit;
	bool decreaseBIC=FALSE; // indicator of if previous bic decreased
	bool tempFlag =FALSE; // indicator of if current bic decreased
	bool isTF =(REAL(dMu)[0]==0); // indicator of if analyzing TF data
	bool finiteBIC =FALSE; // indicator of if any k with finite BIC
	bool okFit =FALSE; // indicator of if any successfully fitted model (infnite bic is allowed)
  char myError[]="";
  SEXP nComp,iiMax,iminReadPerPeak;
  
  PROTECT(nComp=NEW_INTEGER(1));nProtected++;
  PROTECT(iiMax=NEW_INTEGER(1));nProtected++;
  INTEGER(iiMax)[0]=(int)REAL(iMax)[0];  
  PROTECT(iminReadPerPeak=NEW_INTEGER(1));nProtected++;
  INTEGER(iminReadPerPeak)[0]=(int)REAL(minReadPerPeak)[0];
  
	// prepare the output 'ans'
	// I DONT THINK WE NEED TO ALLOCATE THE MEMORY HERE!
  //  PROTECT(ans = allocVector(VECSXP,  4)); nProtected++; 

	for (k=0; k < kmax; k++)
	{

    INTEGER(nComp)[0]=INTEGER(kk)[k];
		//Rprintf("start fit %i mixtures \n", INTEGER(VECTOR_ELT(kk, k))[0]);
		PROTECT(temp = fitModel(nComp,  iiMax,  tol,  mselect,  yR,  yF,  a,  b,  xi,  alpha,  betap,  rho,  lambda,  dMu,  cst,  nu,  iminReadPerPeak));nProtected++;

		//Rprintf("end fit %i mixtures \n", INTEGER(VECTOR_ELT(kk, k))[0]);
		
		if (k==0)
		{ // save it as output in case no model is fitted without error message
			firstFit=temp;
		}
		
		if(strcmp(CHAR(STRING_PTR(VECTOR_ELT(temp, 3))[0]),myError))
		{ // successfully fitted a model, with possible infinite BIC
			okFit=TRUE;
		}
		//Rprintf("k= %d, bic=%lf, bestBIC=%lf \n", k+1, REAL(VECTOR_ELT(temp, 1))[0], bestBIC);
		
		if ( REAL(VECTOR_ELT(temp, 1))[0] > bestBIC ) 
		{ //update 'ans' with new result with better bic
			ans = temp;
			tempFlag  = FALSE;
			finiteBIC = TRUE;
			bestBIC  = REAL(VECTOR_ELT(temp, 1))[0];
			//Rprintf("update new result to ANS \n");
		}
		else if ( ! isTF )
		{    //do nothing for Histone data
			//Rprintf("continue the loop do nothing \n");
			continue;
		}else if (decreaseBIC) { //break loop as see bic decrease twice and is analyzing TF data
			//Rprintf("break the loop\n");
			break;
		}else {//change the flag
			//Rprintf("change decreasing BIC flag \n");
			tempFlag = TRUE;
		}
		decreaseBIC = tempFlag;
		//Rprintf("decreaseBIC= %i", tempFlag);
	}
 
	if (finiteBIC) {
		UNPROTECT(nProtected);
		//Rprintf("return with finite BIC \n\n\n");
		return(ans);
	}else {
		ans=firstFit;
		if (okFit) { // "no finite BIC"
			SET_STRING_ELT(VECTOR_ELT(ans, 3), 0, mkChar("No finite BIC values"));
		}
		UNPROTECT(nProtected);
		//Rprintf("return with no finite BIC \n\n\n");
		return(ans);
	}
}

SEXP fitModel(SEXP kk, SEXP iMax, SEXP tol, SEXP mselect, SEXP yR, SEXP yF, SEXP a, SEXP b, SEXP xi, SEXP alpha, SEXP betap, SEXP rho, SEXP lambda, SEXP dMu, SEXP cst, SEXP nu, SEXP minReadPerPeak)
{
	//paraEM : SEXP iMax, SEXP tol, SEXP mselect
	//segReads: SEXP yR, SEXP yF, SEXP a, SEXP b,
	//paraPrior: SEXP a, SEXP b, SEXP para, SEXP xi, SEXP alpha, SEXP betap, SEXP rho, SEXP lambda, SEXP dMu,
	//constant: SEXP cst, SEXP nu
	//# peaks: SEXP kk
	int NF=length(yF), NR=length(yR), Nmin=imin2(NF,NR), Nmax=imax2(NF,NR), K=INTEGER(kk)[0], nProtected=0, i=0, Kmax=Nmin/(INTEGER(minReadPerPeak)[0]);
	double minweight; 
	SEXP para, iter, ans, bic, converge, err, names;
  
//  Rprintf("I am in fitModel\n");

    // prepare the output 'ans'
	PROTECT(ans		 = allocVector(VECSXP,  4)); nProtected++; 
	PROTECT(bic		 = allocVector(REALSXP, 1)); nProtected++; REAL(bic)[0] = GSL_NEGINF;	
	PROTECT(converge = allocVector(LGLSXP,  1)); nProtected++; LOGICAL(converge)[0]	= FALSE; //not converged
	PROTECT(err		 = allocVector(STRSXP,  1)); nProtected++; SET_STRING_ELT(err, 0, mkChar("Not enough reads")); //error message, 0: "no error", 1: "not enough read", 2: "Numeric error in EM"
	
	PROTECT(names = allocVector(STRSXP, 4)); nProtected++; 
	SET_STRING_ELT(names, 0, mkChar("para"));
	SET_STRING_ELT(names, 1, mkChar("bic"));
	SET_STRING_ELT(names, 2, mkChar("converge"));
	SET_STRING_ELT(names, 3, mkChar("error"));
	
  
	//Handle the case of 'not enough reads'
	if (K>Kmax) { 
		//Rprintf("not enough reads, K=%i/%i=%i is more than upper bound %i \n", Nmin, (INTEGER(minReadPerPeak)[0],K, Kmax);
		SET_VECTOR_ELT(ans, 0, R_NilValue);	//para
		SET_VECTOR_ELT(ans, 1, bic);
		SET_VECTOR_ELT(ans, 2, converge);
		SET_VECTOR_ELT(ans, 3, err);
		setAttrib(ans, R_NamesSymbol, names);
		UNPROTECT(nProtected);
		//Rprintf("not enough reads, return from fitModel \n");
		return(ans);
	}
	
	//Rprintf("Initial parameter\n");
	para=initPara(yF, yR, kk);
	//printPara(para);
	//Rprintf("Start iterEM\n");
	iter=iterEM(iMax, nu, yR, yF, para, xi, alpha, betap, rho, a, b, tol, cst, lambda, dMu);
	//printPara(para);
	//Rprintf("iterEM ok \n");

	//Handle the case of 'smallest weight too small'	
	double *w=REAL(VECTOR_ELT(para, 0));
	minweight=w[0];
	for (i=1; i<K; i++) {
		if (w[i]<minweight) { minweight=w[i]; }
	}
	if (minweight<(1.0/Nmax)) {
		SET_VECTOR_ELT(ans, 0, R_NilValue);	//para
		SET_VECTOR_ELT(ans, 1, bic);
		SET_VECTOR_ELT(ans, 2, converge);
		SET_VECTOR_ELT(ans, 3, err);
		setAttrib(ans, R_NamesSymbol, names);
		UNPROTECT(nProtected);
		//Rprintf("too small weights\n");
		return(ans);
	}
	
	//normal return the final result
	LOGICAL(converge)[0] = (INTEGER(iter)[0]<=INTEGER(iMax)[0]); 
	bic = BIC(nu, yR, yF, para, dMu, lambda, a, b, mselect);

	SET_STRING_ELT(err, 0, mkChar(""));
	SET_VECTOR_ELT(ans, 0, para);//para
	SET_VECTOR_ELT(ans, 1, bic);
	SET_VECTOR_ELT(ans, 2, converge);
	SET_VECTOR_ELT(ans, 3, err);
	setAttrib(ans, R_NamesSymbol, names);
	UNPROTECT(nProtected);
	return(ans);
}	

SEXP iterEM(SEXP iMax, SEXP nu, SEXP yR, SEXP yF, SEXP para, SEXP xi, SEXP alpha, SEXP betap, SEXP rho, SEXP a, SEXP b, SEXP tol, SEXP cst, SEXP lambda, SEXP dMu)
{
	int p = length(VECTOR_ELT(para, 1));
	int i=0, j=0, index[p], Max=INTEGER(iMax)[0];
	//Rprintf("p= %i mixtures, imax=%i \n", p, Max);
	double oldMu[p], w[p], mu[p], delta[p], sF[p], sR[p], sumDiff=0;
	SEXP ans;
	
	/** Turn off the error handler **/
	gsl_set_error_handler_off();
	
	/** Initialize oldMu */
	for(j=0;j<p;j++)
	{
		oldMu[j]=0;
		index[j]=j;
	}
	
	
	for(i=1;i<=Max;i++)
	{
		/** Initialize the difference between new and current values */
		sumDiff=0;
		
		for(j=0;j<p;j++)
		{
			oldMu[j]=REAL(VECTOR_ELT(para, 1))[j];
		}
		
		/** Update the parameters with the two ECM **/
		ECM1(INTEGER(nu)[0], yR, yF, para, REAL(xi)[0], REAL(alpha)[0], REAL(betap)[0], REAL(rho)[0], a, b, REAL(tol)[0], REAL(cst)[0], REAL(lambda)[0], REAL(dMu)[0]);
		ECM2(INTEGER(nu)[0], yR, yF, para, REAL(xi)[0], REAL(alpha)[0], REAL(betap)[0], REAL(rho)[0], a, b, REAL(tol)[0], REAL(cst)[0], REAL(lambda)[0], REAL(dMu)[0]);	  
		
		/** Sort mu and the index for the new estimates*/
		rsort_with_index(REAL(VECTOR_ELT(para, 1)), index, p);
		/** Order all the parameters accordiing to mu */
		for(j=0;j<p;j++)
		{
			w[j]=REAL(VECTOR_ELT(para, 0))[index[j]];
			delta[j]=REAL(VECTOR_ELT(para, 2))[index[j]];
			sF[j]=REAL(VECTOR_ELT(para, 3))[index[j]];
			sR[j]=REAL(VECTOR_ELT(para, 4))[index[j]];
		}
		
		for(j=0;j<p;j++)
		{
			REAL(VECTOR_ELT(para, 0))[j]=w[j];
			REAL(VECTOR_ELT(para, 2))[j]=delta[j];
			REAL(VECTOR_ELT(para, 3))[j]=sF[j];
			REAL(VECTOR_ELT(para, 4))[j]=sR[j];
		}
		
		/** Compute the difference between old and new mu **/
		for(j=0;j<p;j++)
		{
			sumDiff+=fabs(oldMu[j]-REAL(VECTOR_ELT(para, 1))[j]);
			/** Reset the index **/
			index[j]=j;
		}

		/** If the difference is smaller than the tolerance, we stop **/
		if(sumDiff<REAL(tol)[0])
		{
			break;
		}
	}
	/** Return the current index to monitor the convergence **/
	PROTECT(ans=allocVector(INTSXP, 1));
	INTEGER(ans)[0]=i;
	UNPROTECT(1);
	return(ans);
}

void ECM1(int nu, SEXP R, SEXP F, SEXP para, double xi, double alpha, double betap, double rho, SEXP a, SEXP b, double tol, double cst, double lambda, double dMu)
{
	int i=0,j=0, k=0, K=length(VECTOR_ELT(para, 0)), NF=length(F), NR=length(R);
	double *w=REAL(VECTOR_ELT(para, 0)), *mu=REAL(VECTOR_ELT(para, 1)), *delta=REAL(VECTOR_ELT(para, 2)), *sigmaSqF=REAL(VECTOR_ELT(para, 3)), *sigmaSqR=REAL(VECTOR_ELT(para, 4));
	double *yR=REAL(R), *yF=REAL(F), yNormF, yNormR; 
	gsl_vector *sumF=gsl_vector_calloc(NF), *sumR=gsl_vector_calloc(NR);
	gsl_vector *muF=gsl_vector_calloc(K), *muR=gsl_vector_calloc(K);
	gsl_matrix *rF=gsl_matrix_calloc(K,NF), *rR=gsl_matrix_calloc(K,NR);
	gsl_matrix *ruyF=gsl_matrix_calloc(K,NF), *ruyR=gsl_matrix_calloc(K,NR);
	gsl_vector *OneF=gsl_vector_calloc(NF), *OneR=gsl_vector_calloc(NR);
	gsl_vector *chiF=gsl_vector_calloc(K), *chiR=gsl_vector_calloc(K);
	gsl_vector *sF=gsl_vector_calloc(K), *sR=gsl_vector_calloc(K);
	gsl_vector *nF=gsl_vector_calloc(K), *nR=gsl_vector_calloc(K), *ppp=gsl_vector_calloc(K);
	gsl_vector_view yFv, yRv;
	gsl_vector *bb=gsl_vector_calloc(2*K), *sol=gsl_vector_calloc(2*K);
	gsl_matrix *AA=gsl_matrix_calloc(2*K,2*K), *tmp=NULL;
	gsl_permutation *perm=gsl_permutation_alloc(2*K);
	gsl_matrix_view AA_view;
	int s,status1,status2;
	
	/** Initialize the vector of Ones **/
	gsl_vector_set_all(OneF, 1.0);
	gsl_vector_set_all(OneR, 1.0);
	
	/** Turn off the error handler **/
	gsl_set_error_handler_off();
	
	for(j=0;j<K;j++)
	{
		gsl_vector_set(muF,j,mu[j]-delta[j]/2.0);
		gsl_vector_set(muR,j,mu[j]+delta[j]/2.0);
		
		for(i=0;i<NF;i++)
		{
			yNormF=(yF[i]-gsl_vector_get(muF,j))/sqrt(sigmaSqF[j]);
			gsl_matrix_set(rF,j,i,w[j]*gsl_ran_tdist_pdf(yNormF,nu)/sqrt(sigmaSqF[j]));
			/** I track the sum of the posterior prob so that I can normalize them **/
			gsl_vector_set(sumF,i, gsl_vector_get(sumF,i)+gsl_matrix_get(rF,j,i));
		}
		for(i=0;i<NR;i++)
		{
			yNormR=(yR[i]-gsl_vector_get(muR,j))/sqrt(sigmaSqR[j]);
			gsl_matrix_set(rR,j,i,w[j]*gsl_ran_tdist_pdf(yNormR,nu)/sqrt(sigmaSqR[j]));
			/** I track the sum of the posterior prob so that I can normalize them **/
			gsl_vector_set(sumR,i, gsl_vector_get(sumR,i)+gsl_matrix_get(rR,j,i));
		}
	}
	
	for(i=0;i<NF;i++)
	{
		for(j=0;j<K;j++)
		{
			yNormF=(yF[i]-gsl_vector_get(muF,j))/sqrt(sigmaSqF[j]);
			
			/** Normalize **/
			gsl_matrix_set(rF,j,i,gsl_matrix_get(rF,j,i)/gsl_vector_get(sumF,i));
			gsl_matrix_set(ruyF,j,i, gsl_matrix_get(rF,j,i)*(nu+1.)/(nu+yNormF*yNormF));
		}
	}
	
	for(i=0;i<NR;i++)
	{		
		for(j=0;j<K;j++)
		{
			yNormR=(yR[i]-gsl_vector_get(muR,j))/sqrt(sigmaSqR[j]);
			/** Normalize **/
			gsl_matrix_set(rR,j,i,gsl_matrix_get(rR,j,i)/gsl_vector_get(sumR,i));
			gsl_matrix_set(ruyR,j,i, gsl_matrix_get(rR,j,i)*(nu+1.)/(nu+yNormR*yNormR));
		}
	}

	if(length(a)==0)
	{
		gsl_blas_dgemv(CblasNoTrans, 1.0, rF, OneF, 0, chiF);
		gsl_blas_dgemv(CblasNoTrans, 1.0, rR, OneR, 0, chiR);
		
		yFv=gsl_vector_view_array(yF, NF);
		
		yRv=gsl_vector_view_array(yR, NR);
		

		gsl_blas_dgemv(CblasNoTrans, 1.0, ruyF, &yFv.vector, 0, sF);
		gsl_blas_dgemv(CblasNoTrans, 1.0, ruyR, &yRv.vector, 0, sR);
		gsl_blas_dgemv(CblasNoTrans, 1.0, ruyF, OneF, 0, nF);
		gsl_blas_dgemv(CblasNoTrans, 1.0, ruyR, OneR, 0, nR);
		
		
		for(j=0;j<K;j++)
		{
			w[j]=(gsl_vector_get(chiF,j)+gsl_vector_get(chiR,j))/(NF+NR);
			gsl_vector_set(sF,j,gsl_vector_get(sF,j)/sigmaSqF[j]);
			gsl_vector_set(sR,j,gsl_vector_get(sR,j)/sigmaSqR[j]);
			gsl_vector_set(nF,j,gsl_vector_get(nF,j)/sigmaSqF[j]);
			gsl_vector_set(nR,j,gsl_vector_get(nR,j)/sigmaSqR[j]);
			gsl_vector_set(ppp,j,2*rho*(1.0/sigmaSqF[j]+1.0/sigmaSqR[j]));
			gsl_vector_set(bb,j,gsl_vector_get(sF,j)+gsl_vector_get(sR,j));
			gsl_vector_set(bb,K+j,gsl_vector_get(sF,j)-gsl_vector_get(sR,j)-gsl_vector_get(ppp,j)*xi);
		}

		if(K==1)
		{ 
			//Rprintf("k= %i \n",k);
			for(j=0;j<K;j++)
			{
				gsl_matrix_set(AA,j,0,gsl_vector_get(nF,j)+gsl_vector_get(nR,j));
				gsl_matrix_set(AA,j,1,0.5*(gsl_vector_get(nR,j)-gsl_vector_get(nF,j)));
				gsl_matrix_set(AA,K+j,0,gsl_vector_get(nF,j)-gsl_vector_get(nR,j));
				gsl_matrix_set(AA,K+j,1,-0.5*(gsl_vector_get(nR,j)+gsl_vector_get(nF,j))-gsl_vector_get(ppp,j));
			}
		}
		else
		{
			//Rprintf("k= %i \n",k);
			tmp=gsl_matrix_calloc(K,K);
			for(j=0;j<K;j++)
			{
				gsl_matrix_set(tmp,j,j,4*lambda);
			}
			gsl_matrix_set(tmp,0,0,2*lambda);
			gsl_matrix_set(tmp,K-1,K-1,2*lambda);
			gsl_matrix_set(tmp,0,1,-2*lambda);
			gsl_matrix_set(tmp,1,0,-2*lambda);
			if(K>2)
			{
				for(j=1;j<K;j++)
				{
					/** Off diagonals set to -2*lambda **/
					gsl_matrix_set(tmp,j,j-1,-2*lambda);
					gsl_matrix_set(tmp,j-1,j,-2*lambda);
				}
			}
			
			for(j=0;j<K;j++)
			{
				gsl_matrix_set(AA, j, j, gsl_vector_get(nF,j)+gsl_vector_get(nR,j));
				gsl_matrix_set(AA, j, K+j, 0.5*(gsl_vector_get(nR,j)-gsl_vector_get(nF,j)));
				gsl_matrix_set(AA, K+j, j, gsl_vector_get(nF,j)-gsl_vector_get(nR,j));
				gsl_matrix_set(AA, K+j, K+j, -0.5*(gsl_vector_get(nR,j)+gsl_vector_get(nF,j))-gsl_vector_get(ppp,j));
			}
			gsl_vector_set(bb, 0, gsl_vector_get(bb, 0) - 2*dMu*lambda);
			gsl_vector_set(bb, K-1, gsl_vector_get(bb, K-1) + 2*dMu*lambda);
			/** Create a submatrix view from (0,0) of size (K,K) **/
			AA_view=gsl_matrix_submatrix(AA, 0, 0, K, K);
			gsl_matrix_add(&AA_view.matrix, tmp);
			/** Free the memory for tmp **/
			gsl_matrix_free(tmp);
		}
		
		/** Solve the system via LU decomposition **/
		status1=gsl_linalg_LU_decomp(AA, perm, &s);
		status2=gsl_linalg_LU_solve (AA, perm, bb, sol);
		
		/** Copy the new values **/
		/** Check that we could invert the matrix, otherwise we do not update the values **/
		if(status1==0 & status2==0)
		{
			for(j=0;j<K;j++)
			{
				mu[j]=gsl_vector_get(sol,j);
				delta[j]=gsl_vector_get(sol,K+j);
			}
		}
	}
	else /** Missing data **/
	{
		int J=length(a);
		double aNormF, aNormR, bNormF, bNormR, P0F=1, P0R=1, PhiF, PhiR, PhiFw, PhiRw; 
		gsl_matrix *cdfAF=gsl_matrix_calloc(K,J), *cdfBF=gsl_matrix_calloc(K,J), *cdfAR=gsl_matrix_calloc(K,J),*cdfBR=gsl_matrix_calloc(K,J);
		gsl_matrix *aSinF=gsl_matrix_calloc(K,J), *bSinF=gsl_matrix_calloc(K,J), *aSinR=gsl_matrix_calloc(K,J),*bSinR=gsl_matrix_calloc(K,J);
		gsl_matrix *H3F=gsl_matrix_calloc(K,J), *H1F=gsl_matrix_calloc(K,J),  *H0F=gsl_matrix_calloc(K,J);
		gsl_matrix *H3R=gsl_matrix_calloc(K,J), *H1R=gsl_matrix_calloc(K,J),  *H0R=gsl_matrix_calloc(K,J);
		gsl_vector *PjF=gsl_vector_calloc(J), *PjR=gsl_vector_calloc(J), *OneJ=gsl_vector_calloc(J);
		gsl_vector *rowSumH0F=gsl_vector_calloc(K), *rowSumH0R=gsl_vector_calloc(K);
		gsl_vector *rowSumH1F=gsl_vector_calloc(K), *rowSumH1R=gsl_vector_calloc(K);
		gsl_vector *rowSumH3F=gsl_vector_calloc(K), *rowSumH3R=gsl_vector_calloc(K);
		int *aP=INTEGER(a), *bP=INTEGER(b);
		
		/** Initialize the vector of Ones **/
		gsl_vector_set_all(OneJ, 1.0);
		
		
		for(i=0;i<J;i++)
		{
			for(j=0;j<K;j++)
			{
				aNormF=(aP[i]-gsl_vector_get(muF,j))/sqrt(sigmaSqF[j]);
				aNormR=(aP[i]-gsl_vector_get(muR,j))/sqrt(sigmaSqR[j]);
				bNormF=(bP[i]-gsl_vector_get(muF,j))/sqrt(sigmaSqF[j]);
				bNormR=(bP[i]-gsl_vector_get(muR,j))/sqrt(sigmaSqR[j]);
				
				gsl_matrix_set(cdfAF,j,i,gsl_cdf_tdist_P(aNormF,nu));
				gsl_matrix_set(cdfAR,j,i,gsl_cdf_tdist_P(aNormR,nu));
				gsl_matrix_set(cdfBF,j,i,gsl_cdf_tdist_P(bNormF,nu));
				gsl_matrix_set(cdfBR,j,i,gsl_cdf_tdist_P(bNormR,nu));
				
				
				gsl_matrix_set(aSinF,j,i,sin(atan(aNormF/2)));
				gsl_matrix_set(aSinR,j,i,sin(atan(aNormR/2)));
				gsl_matrix_set(bSinF,j,i,sin(atan(bNormF/2)));
				gsl_matrix_set(bSinR,j,i,sin(atan(bNormR/2)));	
				
				
				gsl_matrix_set(H3F,j,i,gsl_matrix_get(cdfBF,j,i)-gsl_matrix_get(cdfAF,j,i));
				gsl_matrix_set(H3R,j,i,gsl_matrix_get(cdfBR,j,i)-gsl_matrix_get(cdfAR,j,i));
				
				gsl_matrix_set(H0F,j,i,cst*(fun0(gsl_matrix_get(bSinF,j,i)) - fun0(gsl_matrix_get(aSinF,j,i))));
				gsl_matrix_set(H0R,j,i,cst*(fun0(gsl_matrix_get(bSinR,j,i)) - fun0(gsl_matrix_get(aSinR,j,i))));
				
				
				gsl_matrix_set(H1F,j,i, (fun1(aNormF,4) - fun1(bNormF,4)) * cst/(nu+1));
				gsl_matrix_set(H1R,j,i, (fun1(aNormR,4) - fun1(bNormR,4)) * cst/(nu+1));
				
				gsl_vector_set(PjF,i, gsl_vector_get(PjF,i)+w[j]*gsl_matrix_get(H3F,j,i));
				gsl_vector_set(PjR,i, gsl_vector_get(PjR,i)+w[j]*gsl_matrix_get(H3R,j,i));
			}
		}
		
		for(i=0;i<J;i++)
		{
			P0F -= gsl_vector_get(PjF,i);
			P0R -= gsl_vector_get(PjR,i);
		}
		PhiF = NF/P0F;
		PhiR = NR/P0R;
		
		
		gsl_blas_dgemv(CblasNoTrans, 1.0, rF, OneF, 0, chiF);
		gsl_blas_dgemv(CblasNoTrans, 1.0, rR, OneR, 0, chiR);
		
		yFv=gsl_vector_view_array(yF, NF);
		yRv=gsl_vector_view_array(yR, NR);
		
		gsl_blas_dgemv(CblasNoTrans, 1.0, ruyF, &yFv.vector, 0, sF);
		gsl_blas_dgemv(CblasNoTrans, 1.0, ruyR, &yRv.vector, 0, sR);
		gsl_blas_dgemv(CblasNoTrans, 1.0, ruyF, OneF, 0, nF);
		gsl_blas_dgemv(CblasNoTrans, 1.0, ruyR, OneR, 0, nR);
		
		gsl_blas_dgemv(CblasNoTrans, 1.0, H0F, OneJ, 0, rowSumH0F);
		gsl_blas_dgemv(CblasNoTrans, 1.0, H0R, OneJ, 0, rowSumH0R);
		gsl_blas_dgemv(CblasNoTrans, 1.0, H1F, OneJ, 0, rowSumH1F);
		gsl_blas_dgemv(CblasNoTrans, 1.0, H1R, OneJ, 0, rowSumH1R);
		gsl_blas_dgemv(CblasNoTrans, 1.0, H3F, OneJ, 0, rowSumH3F);
		gsl_blas_dgemv(CblasNoTrans, 1.0, H3R, OneJ, 0, rowSumH3R);	  
		
		for(j=0;j<K;j++)
		{
			PhiFw = PhiF*w[j];
			PhiRw = PhiR*w[j];
			
			gsl_vector_set(chiF,j,gsl_vector_get(chiF,j)+ PhiFw*gsl_vector_get(rowSumH3F,j));
			gsl_vector_set(chiR,j,gsl_vector_get(chiR,j)+ PhiRw*gsl_vector_get(rowSumH3R,j));
			gsl_vector_set(sF,j,(gsl_vector_get(sF,j)+PhiFw*(gsl_vector_get(muF,j)*gsl_vector_get(rowSumH0F,j)+2*sqrt(sigmaSqF[j])*gsl_vector_get(rowSumH1F,j)))/sigmaSqF[j]);
			gsl_vector_set(sR,j,(gsl_vector_get(sR,j)+PhiRw*(gsl_vector_get(muR,j)*gsl_vector_get(rowSumH0R,j)+2*sqrt(sigmaSqR[j])*gsl_vector_get(rowSumH1R,j)))/sigmaSqR[j]);
			gsl_vector_set(nF,j,(gsl_vector_get(nF,j)+PhiFw*gsl_vector_get(rowSumH0F,j))/sigmaSqF[j]);
			gsl_vector_set(nR,j,(gsl_vector_get(nR,j)+PhiRw*gsl_vector_get(rowSumH0R,j))/sigmaSqR[j]);
			gsl_vector_set(ppp,j,2*rho*(1.0/sigmaSqF[j]+1.0/sigmaSqR[j]));
			gsl_vector_set(bb,j,gsl_vector_get(sF,j)+gsl_vector_get(sR,j));
			gsl_vector_set(bb,K+j,gsl_vector_get(sF,j)-gsl_vector_get(sR,j)-gsl_vector_get(ppp,j)*xi);
			
			w[j]=(gsl_vector_get(chiF,j)+gsl_vector_get(chiR,j))/(PhiF+PhiR);
		}
		
		if(K==1)
		{
			for(j=0;j<K;j++)
			{
				gsl_matrix_set(AA,j,0,gsl_vector_get(nF,j)+gsl_vector_get(nR,j));
				gsl_matrix_set(AA,j,1,0.5*(gsl_vector_get(nR,j)-gsl_vector_get(nF,j)));
				gsl_matrix_set(AA,K+j,0,gsl_vector_get(nF,j)-gsl_vector_get(nR,j));
				gsl_matrix_set(AA,K+j,1,-0.5*(gsl_vector_get(nR,j)+gsl_vector_get(nF,j))-gsl_vector_get(ppp,j));
			}
		}
		else
		{
			tmp=gsl_matrix_calloc(K,K);
			for(j=0;j<K;j++)
			{
				gsl_matrix_set(tmp,j,j,4*lambda);
			}
			gsl_matrix_set(tmp,0,0,2*lambda);
			gsl_matrix_set(tmp,K-1,K-1,2*lambda);
			gsl_matrix_set(tmp,0,1,-2*lambda);
			gsl_matrix_set(tmp,1,0,-2*lambda);
			if(K>2)
			{
				for(j=1;j<K;j++)
				{
					/** Off diagonals set to -2*lambda **/
					gsl_matrix_set(tmp,j,j-1,-2*lambda);
					gsl_matrix_set(tmp,j-1,j,-2*lambda);
				}
			}
			
			for(j=0;j<K;j++)
			{
				gsl_matrix_set(AA, j, j, gsl_vector_get(nF,j)+gsl_vector_get(nR,j));
				gsl_matrix_set(AA, j, K+j, 0.5*(gsl_vector_get(nR,j)-gsl_vector_get(nF,j)));
				gsl_matrix_set(AA, K+j, j, gsl_vector_get(nF,j)-gsl_vector_get(nR,j));
				gsl_matrix_set(AA, K+j, K+j, -0.5*(gsl_vector_get(nR,j)+gsl_vector_get(nF,j))-gsl_vector_get(ppp,j));
			}
			gsl_vector_set(bb, 0, gsl_vector_get(bb, 0) - 2*dMu*lambda);
			gsl_vector_set(bb, K-1, gsl_vector_get(bb, K-1) + 2*dMu*lambda);
			/** Create a submatrix view from (0,0) of size (K,K) **/
			AA_view=gsl_matrix_submatrix(AA, 0, 0, K, K);
			gsl_matrix_add(&AA_view.matrix, tmp);
			/** Free the memory for tmp **/
			gsl_matrix_free(tmp);
		}
		/** Solve the system via LU decomposition **/
		status1=gsl_linalg_LU_decomp(AA, perm, &s);
		status2=gsl_linalg_LU_solve (AA, perm, bb, sol);
		
		/** Copy the new values **/
		/** Check that we could invert the matrix, otherwise we do not update the values **/
		if(status1==0 & status2==0)
		{
			for(j=0;j<K;j++)
			{
				mu[j]=gsl_vector_get(sol,j);
				delta[j]=gsl_vector_get(sol,K+j);
			}
		}	  
		
		/** Free the memory **/ 
		gsl_vector_free(PjF);gsl_vector_free(PjR); gsl_vector_free(OneJ);
		gsl_vector_free(rowSumH0F);  gsl_vector_free(rowSumH1F);      gsl_vector_free(rowSumH3F);
		gsl_vector_free(rowSumH0R);  gsl_vector_free(rowSumH1R);      gsl_vector_free(rowSumH3R);
		gsl_matrix_free(cdfAF); gsl_matrix_free(cdfBF); gsl_matrix_free(cdfAR); gsl_matrix_free(cdfBR);
		gsl_matrix_free(aSinF); gsl_matrix_free(bSinF); gsl_matrix_free(aSinR); gsl_matrix_free(bSinR);
		gsl_matrix_free(H0F); gsl_matrix_free(H1F); gsl_matrix_free(H3F);
		gsl_matrix_free(H0R); gsl_matrix_free(H1R); gsl_matrix_free(H3R);
	}
	
	
	/** Free the memory **/
	gsl_permutation_free(perm);
	gsl_vector_free(sol);gsl_vector_free(bb);
	gsl_vector_free(ppp);
	gsl_matrix_free(AA);
	gsl_vector_free(nF);gsl_vector_free(nR);
	gsl_vector_free(sF);gsl_vector_free(sR);
	gsl_vector_free(sumF);gsl_vector_free(sumR);
	gsl_vector_free(chiF);gsl_vector_free(chiR);
	gsl_vector_free(OneF);gsl_vector_free(OneR);
	gsl_vector_free(muF);gsl_vector_free(muR);
	gsl_matrix_free(rF);gsl_matrix_free(rR);
	gsl_matrix_free(ruyF);gsl_matrix_free(ruyR);
}

void ECM2(int nu, SEXP R, SEXP F, SEXP para, double xi, double alpha, double betap, double rho, SEXP a, SEXP b, double tol, double cst, double lambda, double dMu)
{
	int i=0,j=0, k=0, K=length(VECTOR_ELT(para, 0)), NF=length(F), NR=length(R);
	double *w=REAL(VECTOR_ELT(para, 0)), *mu=REAL(VECTOR_ELT(para, 1)), *delta=REAL(VECTOR_ELT(para, 2)), *sigmaSqF=REAL(VECTOR_ELT(para, 3)), *sigmaSqR=REAL(VECTOR_ELT(para, 4));
	double *yR=REAL(R), *yF=REAL(F), yNormF, yNormR; 
	gsl_vector *sumF=gsl_vector_calloc(NF), *sumR=gsl_vector_calloc(NR);
	gsl_vector *muF=gsl_vector_calloc(K), *muR=gsl_vector_calloc(K);
	gsl_matrix *rF=gsl_matrix_calloc(K,NF), *rR=gsl_matrix_calloc(K,NR);
	gsl_matrix *ruyF=gsl_matrix_calloc(K,NF), *ruyR=gsl_matrix_calloc(K,NR);
	gsl_vector *OneF=gsl_vector_calloc(NF), *OneR=gsl_vector_calloc(NR);
	gsl_vector *chiF=gsl_vector_calloc(K), *chiR=gsl_vector_calloc(K); 
	double chiSum, etaF, etaR, etaDiff, dd, aaF, aaR, ggF, ggR, cc, ee;
	
	/** Initialize the vector of Ones **/
	gsl_vector_set_all(OneF, 1.0);
	gsl_vector_set_all(OneR, 1.0);
	
	for(j=0;j<K;j++)
	{
		gsl_vector_set(muF,j,mu[j]-delta[j]/2);
		gsl_vector_set(muR,j,mu[j]+delta[j]/2);
		
		for(i=0;i<NF;i++)
		{
			yNormF=(yF[i]-gsl_vector_get(muF,j))/sqrt(sigmaSqF[j]);
			gsl_matrix_set(rF,j,i,w[j]*gsl_ran_tdist_pdf(yNormF,nu)/sqrt(sigmaSqF[j]));
			gsl_vector_set(sumF,i, gsl_vector_get(sumF,i)+gsl_matrix_get(rF,j,i));
		}
		for(i=0;i<NR;i++)
		{
			yNormR=(yR[i]-gsl_vector_get(muR,j))/sqrt(sigmaSqR[j]);
			gsl_matrix_set(rR,j,i,w[j]*gsl_ran_tdist_pdf(yNormR,nu)/sqrt(sigmaSqR[j]));
			gsl_vector_set(sumR,i, gsl_vector_get(sumR,i)+gsl_matrix_get(rR,j,i));
		}
	}
	
	for(i=0;i<NF;i++)
	{
		for(j=0;j<K;j++)
		{
			yNormF=(yF[i]-gsl_vector_get(muF,j))/sqrt(sigmaSqF[j]);
			gsl_matrix_set(rF,j,i,gsl_matrix_get(rF, j, i)/gsl_vector_get(sumF,i));
			gsl_matrix_set(ruyF,j,i, gsl_matrix_get(rF, j, i)*(nu+1.)/(nu+yNormF*yNormF));
		}
	}
	
	for(i=0;i<NR;i++)
	{
		for(j=0;j<K;j++)
		{
			yNormR=(yR[i]-gsl_vector_get(muR,j))/sqrt(sigmaSqR[j]);
			gsl_matrix_set(rR,j,i,gsl_matrix_get(rR, j, i)/gsl_vector_get(sumR,i));
			gsl_matrix_set(ruyR,j,i, gsl_matrix_get(rR, j, i)*(nu+1.)/(nu+yNormR*yNormR));
		}
	}
	
	if(length(a)==0)
	{
		gsl_blas_dgemv(CblasNoTrans, 1.0, rF, OneF, 0, chiF);
		gsl_blas_dgemv(CblasNoTrans, 1.0, rR, OneR, 0, chiR);
		cc=2*alpha-1;
		for(j=0;j<K;j++)
		{
			chiSum=gsl_vector_get(chiF,j)+gsl_vector_get(chiR,j);
			/** Initialize the sums **/
			etaF=0;etaR=0;
			for(i=0;i<NF;i++)
			{
				etaF+=gsl_matrix_get(ruyF,j,i)*(yF[i]-gsl_vector_get(muF,j))*(yF[i]-gsl_vector_get(muF,j));
			}
			for(i=0;i<NR;i++)
			{
				etaR+=gsl_matrix_get(ruyR,j,i)*(yR[i]-gsl_vector_get(muR,j))*(yR[i]-gsl_vector_get(muR,j));
			}
			dd=rho*(delta[j]-xi)*(delta[j]-xi)+2*betap;
			
      sigmaSqF[j]=(etaF+dd)/(cc+gsl_vector_get(chiF,j));
      sigmaSqR[j]=(etaR+dd)/(cc+gsl_vector_get(chiR,j));
		}
	}
	else 
	{	  
		int J=length(a);
		double aNormF, aNormR, bNormF, bNormR, P0F=1, P0R=1, PhiF, PhiR, PhiFw, PhiRw; 
		gsl_matrix *cdfAF=gsl_matrix_calloc(K,J), *cdfBF=gsl_matrix_calloc(K,J), *cdfAR=gsl_matrix_calloc(K,J),*cdfBR=gsl_matrix_calloc(K,J);
		gsl_matrix *aSinF=gsl_matrix_calloc(K,J), *bSinF=gsl_matrix_calloc(K,J), *aSinR=gsl_matrix_calloc(K,J),*bSinR=gsl_matrix_calloc(K,J);
		gsl_matrix *H3F=gsl_matrix_calloc(K,J), *H2F=gsl_matrix_calloc(K,J);
		gsl_matrix *H3R=gsl_matrix_calloc(K,J), *H2R=gsl_matrix_calloc(K,J);
		gsl_vector *PjF=gsl_vector_calloc(J), *PjR=gsl_vector_calloc(J), *OneJ=gsl_vector_calloc(J);
		gsl_vector *rowSumH2F=gsl_vector_calloc(K), *rowSumH2R=gsl_vector_calloc(K);
		gsl_vector *rowSumH3F=gsl_vector_calloc(K), *rowSumH3R=gsl_vector_calloc(K);
		int  *aP=INTEGER(a), *bP=INTEGER(b); 	  
		
		/** Initialize the vector of Ones **/
		gsl_vector_set_all(OneJ, 1.0);
		
		
		for(i=0;i<J;i++)
		{
			for(j=0;j<K;j++)
			{
				aNormF=(aP[i]-gsl_vector_get(muF,j))/sqrt(sigmaSqF[j]);
				aNormR=(aP[i]-gsl_vector_get(muR,j))/sqrt(sigmaSqR[j]);
				bNormF=(bP[i]-gsl_vector_get(muF,j))/sqrt(sigmaSqF[j]);
				bNormR=(bP[i]-gsl_vector_get(muR,j))/sqrt(sigmaSqR[j]);
				
				gsl_matrix_set(cdfAF,j,i,gsl_cdf_tdist_P(aNormF,nu));
				gsl_matrix_set(cdfAR,j,i,gsl_cdf_tdist_P(aNormR,nu));
				gsl_matrix_set(cdfBF,j,i,gsl_cdf_tdist_P(bNormF,nu));
				gsl_matrix_set(cdfBR,j,i,gsl_cdf_tdist_P(bNormR,nu));			  
				
				
				gsl_matrix_set(aSinF,j,i,sin(atan(aNormF/2)));
				gsl_matrix_set(aSinR,j,i,sin(atan(aNormR/2)));
				gsl_matrix_set(bSinF,j,i,sin(atan(bNormF/2)));
				gsl_matrix_set(bSinR,j,i,sin(atan(bNormR/2)));	
				
				
				gsl_matrix_set(H3F,j,i,gsl_matrix_get(cdfBF,j,i)-gsl_matrix_get(cdfAF,j,i));
				gsl_matrix_set(H3R,j,i,gsl_matrix_get(cdfBR,j,i)-gsl_matrix_get(cdfAR,j,i));
				
				
				gsl_matrix_set(H2F,j,i,cst*(fun2(gsl_matrix_get(bSinF,j,i)) - fun2(gsl_matrix_get(aSinF,j,i))));
				gsl_matrix_set(H2R,j,i,cst*(fun2(gsl_matrix_get(bSinR,j,i)) - fun2(gsl_matrix_get(aSinR,j,i))));
				
				gsl_vector_set(PjF,i, gsl_vector_get(PjF,i)+w[j]*gsl_matrix_get(H3F,j,i));
				gsl_vector_set(PjR,i, gsl_vector_get(PjR,i)+w[j]*gsl_matrix_get(H3R,j,i));
			}
		}
		
		for(i=0;i<J;i++)
		{
			P0F -= gsl_vector_get(PjF,i);
			P0R -= gsl_vector_get(PjR,i);
		}
		PhiF = NF/P0F;
		PhiR = NR/P0R;
		
		
		gsl_blas_dgemv(CblasNoTrans, 1.0, rF, OneF, 0, chiF);
		gsl_blas_dgemv(CblasNoTrans, 1.0, rR, OneR, 0, chiR);
		gsl_blas_dgemv(CblasNoTrans, 1.0, H2F, OneJ, 0, rowSumH2F);
		gsl_blas_dgemv(CblasNoTrans, 1.0, H2R, OneJ, 0, rowSumH2R);
		gsl_blas_dgemv(CblasNoTrans, 1.0, H3F, OneJ, 0, rowSumH3F);
		gsl_blas_dgemv(CblasNoTrans, 1.0, H3R, OneJ, 0, rowSumH3R);	  
		
		
		cc=2*alpha-1;
		for(j=0;j<K;j++)
		{
			PhiFw = PhiF*w[j];
			PhiRw = PhiR*w[j];
			
			gsl_vector_set(chiF,j,gsl_vector_get(chiF,j)+ PhiFw*gsl_vector_get(rowSumH3F,j));
			gsl_vector_set(chiR,j,gsl_vector_get(chiR,j)+ PhiRw*gsl_vector_get(rowSumH3R,j));
			chiSum=gsl_vector_get(chiF,j)+gsl_vector_get(chiR,j);
			
			/** Initialize the sums **/
			etaF=0;etaR=0;
			for(i=0;i<NF;i++)
			{
				etaF+=gsl_matrix_get(ruyF,j,i)*(yF[i]-gsl_vector_get(muF,j))*(yF[i]-gsl_vector_get(muF,j));
			}
			for(i=0;i<NR;i++)
			{
				etaR+=gsl_matrix_get(ruyR,j,i)*(yR[i]-gsl_vector_get(muR,j))*(yR[i]-gsl_vector_get(muR,j));
			}
			etaF += PhiFw*sigmaSqF[j]*nu*gsl_vector_get(rowSumH2F,j);
			etaR += PhiRw*sigmaSqR[j]*nu*gsl_vector_get(rowSumH2R,j);
			dd=rho*(delta[j]-xi)*(delta[j]-xi)+2*betap;
			
      // if(priortype==1) 
      // {/* update sigma using prior on sigmaSqF^(-2),sigmaSqR^(-2) separately*/
				sigmaSqF[j]=(etaF+dd)/(cc+gsl_vector_get(chiF,j));
				sigmaSqR[j]=(etaR+dd)/(cc+gsl_vector_get(chiR,j));
      // }else {/* update sigma using prior on sigmaSqF^(-2)+sigmaSqR^(-2)*/
      //  etaDiff=etaF-etaR;
      //  // ddetaF and ddetaR 
      //  etaF+=dd;
      //  etaR+=dd;
      //  aaF=etaDiff*(cc+gsl_vector_get(chiF,j))+etaF*chiSum;
      //  aaR=-etaDiff*(cc+gsl_vector_get(chiR,j))+etaR*chiSum;
      //  ggF=2*etaF*etaDiff;
      //  ggR=2*etaR*etaDiff;
      //  ee=sqrt((-cc*etaDiff+etaR*gsl_vector_get(chiF,j))*(-cc*etaDiff+etaR*gsl_vector_get(chiF,j))
      //      +2*(cc*etaDiff+etaR*gsl_vector_get(chiF,j))*etaF*gsl_vector_get(chiR,j)+etaF*etaF*gsl_vector_get(chiR,j)*gsl_vector_get(chiR,j));
      //  sigmaSqF[j]=ggF/(aaF-ee);
      //  sigmaSqR[j]=ggR/(-aaR+ee);
      // }
		}
		/** Free the memory **/ 
		gsl_vector_free(PjF);gsl_vector_free(PjR); gsl_vector_free(OneJ); 
		gsl_vector_free(rowSumH2F); gsl_vector_free(rowSumH3F);
		gsl_vector_free(rowSumH2R); gsl_vector_free(rowSumH3R);
		gsl_matrix_free(cdfAF); gsl_matrix_free(cdfBF); gsl_matrix_free(cdfAR); gsl_matrix_free(cdfBR);
		gsl_matrix_free(aSinF); gsl_matrix_free(bSinF); gsl_matrix_free(aSinR); gsl_matrix_free(bSinR);
		gsl_matrix_free(H2F); gsl_matrix_free(H3F);
		gsl_matrix_free(H2R); gsl_matrix_free(H3R);
	}
	
	gsl_vector_free(sumF);gsl_vector_free(sumR);
	gsl_vector_free(chiF);gsl_vector_free(chiR);
	gsl_vector_free(OneF);gsl_vector_free(OneR);
	gsl_vector_free(muF);gsl_vector_free(muR);
	gsl_matrix_free(rF);gsl_matrix_free(rR);
	gsl_matrix_free(ruyF);gsl_matrix_free(ruyR);
}

double logDensityMix(double *Y, double *w, double *mu, double *sigmaSq, int K, int N){
	int nu=4, i, j;
	double ans=0.0, sigma, tmp, yNorm;
	
	for(i=0; i<N; i++){
		tmp=0;
		for (j=0; j<K; j++){
			sigma=sqrt(sigmaSq[j]);
			yNorm=(Y[i]-mu[j])/sigma;
			tmp += w[j]*gsl_ran_tdist_pdf(yNorm,nu)/sigma;
		}
		ans += log(tmp);
	}
	return(ans);
}

double plogp(double y, void * params){
	double *paras = (double *) params;
	int i, K=(int) paras[0];
	double w[K], mu[K], sigmaSq[K];

	for (i=0; i<K; i++) 
	{
		w[i]= paras[i+1];
		mu[i]=paras[i+K+1];
		sigmaSq[i]= paras[i+2*K+1];
	}
	
	int nu=4, j;
	double ans=0.0, sigma, tmp=0.0, yNorm;
	
	for (j=0; j<K; j++){
		sigma=sqrt(sigmaSq[j]);
		yNorm=(y-mu[j])/sigma;
		tmp += w[j]*gsl_ran_tdist_pdf(yNorm,nu)/sigma;			
	}
	ans = tmp*log(tmp);
	return(ans);
}

double intplogp (int start, int end, double *w, double *mu, double *sigmaSq, int K)
{
  gsl_integration_workspace * wk = gsl_integration_workspace_alloc(1000);
  double result, error;
  int i=0;

  double paras[3*K+1];
  paras[0]= (double)K;
  for (i=0; i<K; i++) {
    paras[i+1]=w[i];
    paras[K+i+1]=mu[i];
    paras[2*K+i+1]=sigmaSq[i];
  }

  gsl_function FFF;
  FFF.function = &plogp;
  FFF.params = &paras;

  gsl_integration_qags (&FFF, start, end, 0, 1e-7, 1000, wk, &result, &error); 

  gsl_integration_workspace_free (wk);

  return result;
}

SEXP BIC(SEXP nuC, SEXP R, SEXP F, SEXP para, SEXP dMuC, SEXP lambdaC, SEXP a, SEXP b, SEXP mselect){
	int  nu=INTEGER(nuC)[0], K=length(VECTOR_ELT(para, 0)), NF=length(F), NR=length(R), i=0, j=0, k=0, N, type=INTEGER(mselect)[0];
	double *w=REAL(VECTOR_ELT(para, 0)), *mu=REAL(VECTOR_ELT(para, 1)), *delta=REAL(VECTOR_ELT(para, 2)), *sigmaSqF=REAL(VECTOR_ELT(para, 3)), *sigmaSqR=REAL(VECTOR_ELT(para, 4));
	double *yR=REAL(R), *yF=REAL(F), dMu=REAL(dMuC)[0], lambda=REAL(lambdaC)[0]; 
	double bic, tmp, muF[K], muR[K], penalty, distance;
	
	
	
	for (i=0; i<K; i++) {
		muF[i] = mu[i] - delta[i]/2.0;
		muR[i] = mu[i] + delta[i]/2.0;
	}
	
	
	N=NF+NR;
	switch (type) {
		case 1: //BIC
			penalty=(5.0*K-1.0)*log(N)/2.0;
			break;
		case 2: //AIC
			penalty=(5*K-1);
			break;
		case 3://Patial prior of mu "sum_i(diff(mu)_i-dMu)^2"
			tmp=0;
			for (i=1; i<K; i++) {
				distance = mu[i] - mu[(i-1)] - dMu;
				tmp += pow(distance,2);
			}		
			penalty = -lambda * tmp;
			break;
		default:
			//puts("Only 1,2,3 are allowed.");
			break;
	}
  // Rprintf("type= %d , penalty= %lf \n", type, penalty);
	
	bic = logDensityMix(yF, w, muF, sigmaSqF, K, NF) + logDensityMix(yR, w, muR, sigmaSqR, K,NR) - penalty;
	
	
	//Rprintf("bic= %lf, length(a)=%d, lambda=%lf, dMu=%lf, penalty=%lf, logPF=%lf, logPR=%lf,  \n", bic, length(a), *lambda, *dMu, penalty,  logDensityMix(FF, ww, muF, sigmaSqFF, K, NF), logDensityMix(RR, ww, muR, sigmaSqRR, K,NR));
	
	if(length(a)>0){
		int J=length(a);
		double aNormF, aNormR, bNormF, bNormR, P0F=1, P0R=1, integrateF=0.0, integrateR=0.0; 
		gsl_matrix *cdfAF=gsl_matrix_calloc(K,J), *cdfBF=gsl_matrix_calloc(K,J);
		gsl_matrix *cdfAR=gsl_matrix_calloc(K,J),*cdfBR=gsl_matrix_calloc(K,J);
		gsl_matrix *H3F=gsl_matrix_calloc(K,J), *H3R=gsl_matrix_calloc(K,J);
		gsl_vector *PjF=gsl_vector_calloc(J), *PjR=gsl_vector_calloc(J);
		int *aP=INTEGER(a), *bP=INTEGER(b);
		
		
		for(i=0;i<J;i++)
		{
			for(j=0;j<K;j++)
			{
				aNormF=(aP[i]-muF[j])/sqrt(sigmaSqF[j]);
				aNormR=(aP[i]-muR[j])/sqrt(sigmaSqR[j]);
				bNormF=(bP[i]-muF[j])/sqrt(sigmaSqF[j]);
				bNormR=(bP[i]-muR[j])/sqrt(sigmaSqR[j]);
				
				gsl_matrix_set(cdfAF,j,i,gsl_cdf_tdist_P(aNormF,nu));
				gsl_matrix_set(cdfAR,j,i,gsl_cdf_tdist_P(aNormR,nu));
				gsl_matrix_set(cdfBF,j,i,gsl_cdf_tdist_P(bNormF,nu));
				gsl_matrix_set(cdfBR,j,i,gsl_cdf_tdist_P(bNormR,nu));			  
				
				gsl_matrix_set(H3F,j,i,gsl_matrix_get(cdfBF,j,i)-gsl_matrix_get(cdfAF,j,i));
				gsl_matrix_set(H3R,j,i,gsl_matrix_get(cdfBR,j,i)-gsl_matrix_get(cdfAR,j,i));
				
				gsl_vector_set(PjF,i, gsl_vector_get(PjF,i)+w[j]*gsl_matrix_get(H3F,j,i));
				gsl_vector_set(PjR,i, gsl_vector_get(PjR,i)+w[j]*gsl_matrix_get(H3R,j,i));
			}
		}
		
		for(i=0;i<J;i++)
		{
			P0F -= gsl_vector_get(PjF,i);
			P0R -= gsl_vector_get(PjR,i);
			integrateF += intplogp (aP[i], bP[i], w, muF, sigmaSqF, K);
			integrateR += intplogp (aP[i], bP[i], w, muR, sigmaSqR, K);
		}
		
		bic += NF/P0F*integrateF + NR/P0R*integrateR;
		
		// Free the memory 
		gsl_vector_free(PjF);  gsl_vector_free(PjR); gsl_matrix_free(H3F);  gsl_matrix_free(H3R);
		gsl_matrix_free(cdfAF); gsl_matrix_free(cdfBF); gsl_matrix_free(cdfAR); gsl_matrix_free(cdfBR);
	}
	
	SEXP ans;
	PROTECT(ans=allocVector(REALSXP, 1));
	REAL(ans)[0]=bic;
	UNPROTECT(1);
	return(ans);
	
}

SEXP initPara(SEXP F, SEXP R, SEXP kk) {
	int NF=length(F), NR=length(R), twoK, i=0, K=INTEGER(kk)[0];
	double *yF=REAL(F), *yR=REAL(R), weights=1.0/K,  varF, varR, muR[K], muF[K];

	gsl_sort(yF, 1, NF);
	gsl_sort(yR, 1, NR);
	varF = gsl_stats_variance(yF,1,NF)/K;
	varR = gsl_stats_variance(yR,1,NR)/K;	
	
	//Rprintf("K=%d, weights=%lf\n", K, weights);
	 
	int nProtected=0; //number of protected variables
	
	SEXP w, mu, delta, sigmaSqF, sigmaSqR; 
	PROTECT(w			= allocVector(REALSXP, K));	nProtected++; 
	PROTECT(mu			= allocVector(REALSXP, K));	nProtected++; 
	PROTECT(delta		= allocVector(REALSXP, K));	nProtected++; 
	PROTECT(sigmaSqF	= allocVector(REALSXP, K));	nProtected++; 
	PROTECT(sigmaSqR	= allocVector(REALSXP, K));	nProtected++; 
	
	twoK= K*2;
	for (i=0; i<K; i++) {
		muR[i] = gsl_stats_quantile_from_sorted_data(yR, 1, NR, (2.0*i+1.0)/twoK);
		muF[i] = gsl_stats_quantile_from_sorted_data(yF, 1, NF, (2.0*i+1.0)/twoK);
		REAL(mu)[i]			= (muF[i]+muR[i])/2;
		REAL(w)[i]			= weights;
		REAL(delta)[i]		= 150.0;
		REAL(sigmaSqF)[i]	= varF;
		REAL(sigmaSqR)[i]	= varR;
	}
	
	SEXP para;
	PROTECT(para = allocVector(VECSXP, 5)); nProtected++; 
	//PROTECT(para = NEW_LIST(5)); 
	SET_VECTOR_ELT(para, 0, w);
	SET_VECTOR_ELT(para, 1, mu);
	SET_VECTOR_ELT(para, 2, delta);
	SET_VECTOR_ELT(para, 3, sigmaSqF);
	SET_VECTOR_ELT(para, 4, sigmaSqR);


	SEXP names; 
	PROTECT(names = allocVector(STRSXP, 5)); nProtected++; 
	SET_STRING_ELT(names, 0, mkChar("w"));
	SET_STRING_ELT(names, 1, mkChar("mu"));
	SET_STRING_ELT(names, 2, mkChar("delta"));
	SET_STRING_ELT(names, 3, mkChar("sigmaSqF"));
	SET_STRING_ELT(names, 4, mkChar("sigmaSqR"));
	setAttrib(para, R_NamesSymbol, names);
	UNPROTECT(nProtected);
	return(para);
}

void printPara(SEXP para){
	double *w=REAL(VECTOR_ELT(para, 0)), *mu=REAL(VECTOR_ELT(para, 1)), *delta=REAL(VECTOR_ELT(para, 2)), *sigmaSqF=REAL(VECTOR_ELT(para, 3)), *sigmaSqR=REAL(VECTOR_ELT(para, 4));
	int i, K=length(VECTOR_ELT(para, 0));
	
	Rprintf("w=");
	for (i=0; i<K; i++) {
		Rprintf("%lf \t ",w[i]);
	}
	Rprintf("\n");
	
	Rprintf("mu=");
	for (i=0; i<K; i++) {
		Rprintf("%lf \t ",mu[i]);
	}
	Rprintf("\n");
	
	Rprintf("delta=");
	for (i=0; i<K; i++) {
		Rprintf("%lf \t ",delta[i]);
	}
	Rprintf("\n");	
	
	Rprintf("sigmaSqF=");
	for (i=0; i<K; i++) {
		Rprintf("%lf \t ",sigmaSqF[i]);
	}
	Rprintf("\n");
	
	Rprintf("sigmaSqR=");
	for (i=0; i<K; i++) {
		Rprintf("%lf \t ",sigmaSqR[i]);
	}
	Rprintf("\n");
}


int getInfMat(SEXP R, SEXP F, SEXP para, SEXP a, SEXP b, double rho, double xi, double alpha, double cst, double lambda, double nu, gsl_matrix* infMat, gsl_vector* se, gsl_vector* seF, gsl_vector* seR)
{
  int i=0,j=0, k=0, K=length(VECTOR_ELT(para, 0)), NF=length(F), NR=length(R);
  double *w=REAL(VECTOR_ELT(para, 0)), *mu=REAL(VECTOR_ELT(para, 1)), *delta=REAL(VECTOR_ELT(para, 2)), *sigmaSqF=REAL(VECTOR_ELT(para, 3)), *sigmaSqR=REAL(VECTOR_ELT(para, 4));
  double *yR=REAL(R), *yF=REAL(F), yNormF, yNormR; 
  gsl_vector *sumF=gsl_vector_calloc(NF), *sumR=gsl_vector_calloc(NR);
  gsl_vector *muF=gsl_vector_calloc(K), *muR=gsl_vector_calloc(K);
  gsl_matrix *rF=gsl_matrix_calloc(K,NF), *rR=gsl_matrix_calloc(K,NR);
  gsl_matrix *ruyF=gsl_matrix_calloc(K,NF), *ruyR=gsl_matrix_calloc(K,NR);
  gsl_vector *OneF=gsl_vector_calloc(NF), *OneR=gsl_vector_calloc(NR);
  gsl_matrix *qWoF=gsl_matrix_calloc(K,NF),*qWoR=gsl_matrix_calloc(K,NR);
  gsl_matrix *qMuoF=gsl_matrix_calloc(K,NF),*qMuoR=gsl_matrix_calloc(K,NR);
  gsl_matrix *qSigmaoF=gsl_matrix_calloc(K,NF),*qSigmaoR=gsl_matrix_calloc(K,NR);
  gsl_matrix *scoreObsF=gsl_matrix_calloc(5*K-1,NF),*scoreObsR=gsl_matrix_calloc(5*K-1,NR);
  gsl_vector *scoreObsBarF=gsl_vector_calloc(5*K-1), *scoreObsBarR=gsl_vector_calloc(5*K-1);
  gsl_matrix *infMatF=gsl_matrix_calloc(5*K-1,5*K-1), *infMatR=gsl_matrix_calloc(5*K-1,5*K-1), *infMatPri=gsl_matrix_calloc(5*K-1,5*K-1);
  gsl_vector_view scoreObsF_col, scoreObsR_col;
  gsl_vector *A=gsl_vector_calloc(5*K-1);
  gsl_matrix *DiagOne=gsl_matrix_calloc(5*K-1,5*K-1);
  int flag=0;
  
 /** Initialize the vector of Ones **/
 gsl_vector_set_all(OneF, 1.0);
 gsl_vector_set_all(OneR, 1.0);
 
 /** Turn off the error handler **/
 gsl_set_error_handler_off();
 
 for(j=0;j<K;j++)
 {
   gsl_vector_set(muF,j,mu[j]-delta[j]/2.0);
   gsl_vector_set(muR,j,mu[j]+delta[j]/2.0);
   
   for(i=0;i<NF;i++)
   {
     yNormF=(yF[i]-gsl_vector_get(muF,j))/sqrt(sigmaSqF[j]);
     gsl_matrix_set(rF,j,i,w[j]*gsl_ran_tdist_pdf(yNormF,nu)/sqrt(sigmaSqF[j]));
     /** I track the sum of the posterior prob so that I can normalize them **/
     gsl_vector_set(sumF,i, gsl_vector_get(sumF,i)+gsl_matrix_get(rF,j,i));
   }
   for(i=0;i<NR;i++)
   {
     yNormR=(yR[i]-gsl_vector_get(muR,j))/sqrt(sigmaSqR[j]);
     gsl_matrix_set(rR,j,i,w[j]*gsl_ran_tdist_pdf(yNormR,nu)/sqrt(sigmaSqR[j]));
     /** I track the sum of the posterior prob so that I can normalize them **/
     gsl_vector_set(sumR,i, gsl_vector_get(sumR,i)+gsl_matrix_get(rR,j,i));
   }
 }
 
 for(i=0;i<NF;i++)
 {
   for(j=0;j<K;j++)
   {
     yNormF=(yF[i]-gsl_vector_get(muF,j))/sqrt(sigmaSqF[j]);
     
     /** Normalize **/
     gsl_matrix_set(rF,j,i,gsl_matrix_get(rF,j,i)/gsl_vector_get(sumF,i));
     gsl_matrix_set(ruyF,j,i, gsl_matrix_get(rF,j,i)*(nu+1.)/(nu+yNormF*yNormF));
     //'q.w.o[k,n]=dQ_{k,obs}(y_n,theta)/d w_k', a matrix of 'K' by 'N'
     gsl_matrix_set(qWoF,j,i, gsl_matrix_get(rF,j,i)/w[j]);
     //'qMu.o[k,n]=dQ_{k,obs}(y_n,theta)/d mu_k', a matrix of 'K' by 'N'
     gsl_matrix_set(qMuoF,j,i, gsl_matrix_get(ruyF,j,i)*yNormF/sqrt(sigmaSqF[j]));
      //'q.sigma.o[k,n]=dQ_{k,obs}(y_n,theta)/d sigma.sq_k', a matrix of 'K' by 'N'
     gsl_matrix_set(qSigmaoF,j,i, sigmaSqF[j]/2.*(gsl_matrix_get(rF,j,i)-gsl_matrix_get(ruyF,j,i)*yNormF*yNormF));
   }
 }
 
 for(i=0;i<NR;i++)
 {
   for(j=0;j<K;j++)
   {
     yNormR=(yR[i]-gsl_vector_get(muR,j))/sqrt(sigmaSqR[j]);
     /** Normalize **/
     gsl_matrix_set(rR,j,i,gsl_matrix_get(rR,j,i)/gsl_vector_get(sumR,i));
     gsl_matrix_set(ruyR,j,i, gsl_matrix_get(rR,j,i)*(nu+1.)/(nu+yNormR*yNormR));
     //'q.w.o[k,n]=dQ_{k,obs}(y_n,theta)/d w_k', a matrix of 'K' by 'N'      
     gsl_matrix_set(qWoR,j,i, gsl_matrix_get(rR,j,i)/w[j]);
     //'qMu.o[k,n]=dQ_{k,obs}(y_n,theta)/d mu_k', a matrix of 'K' by 'N'
      gsl_matrix_set(qMuoR,j,i, gsl_matrix_get(ruyR,j,i)*yNormR/sqrt(sigmaSqR[j]));
      //'q.sigma.o[k,n]=dQ_{k,obs}(y_n,theta)/d sigma.sq_k', a matrix of 'K' by 'N'
      gsl_matrix_set(qSigmaoR,j,i, sigmaSqR[j]/2.*(gsl_matrix_get(rR,j,i)-gsl_matrix_get(ruyR,j,i)*yNormR*yNormR));
   }
 }

  //the score matrix for observed data, a matrix of '5*K-1' by 'N'
  for(i=0;i<NF;i++)
  {
    for(j=1;j<K;j++)
    {
      gsl_matrix_set(scoreObsF, j-1, i, gsl_matrix_get(qWoF,j,i));
    }
    for(j=0;j<K;j++)
    {
      gsl_matrix_set(scoreObsF, (K-1)+j, i, gsl_matrix_get(qMuoF,j,i));
      gsl_matrix_set(scoreObsF, 2*K-1+j, i, -gsl_matrix_get(qMuoF,j,i)/2.);
      gsl_matrix_set(scoreObsF, 3*K-1+j, i, gsl_matrix_get(qSigmaoF,j,i));
    }
  }
  for(i=0;i<NR;i++)
  {
    for(j=1;j<K;j++)
    {
      gsl_matrix_set(scoreObsR, j-1, i, gsl_matrix_get(qWoR,j,i));
    }
    for(j=0;j<K;j++)
    {
      gsl_matrix_set(scoreObsR, (K-1)+j, i, gsl_matrix_get(qMuoR,j,i));
      gsl_matrix_set(scoreObsR, 2*K-1+j, i, gsl_matrix_get(qMuoR,j,i)/2.);
      gsl_matrix_set(scoreObsR, 4*K-1+j, i, gsl_matrix_get(qSigmaoR,j,i));
    }
  }

    // Compute scoreObs*1 (sum over rows)
  gsl_blas_dgemv(CblasNoTrans, 1.0, scoreObsF, OneF, 0, scoreObsBarF);
  gsl_blas_dgemv(CblasNoTrans, 1.0, scoreObsR, OneR, 0, scoreObsBarR);
    // Scale the matrix to compute the mean
  gsl_vector_scale(scoreObsBarF, 1./NF);
  gsl_vector_scale(scoreObsBarR, 1./NR);
  
    //the information matrix, a matrix of '5*K-1' by '5*K-1'
  for(i=0;i<NF;i++)
  {
    scoreObsF_col=gsl_matrix_column(scoreObsF, i);
    gsl_blas_dsyr(CblasLower, 1.0, &scoreObsF_col.vector, infMatF);
  }
  for(i=0;i<NR;i++)
  {
    scoreObsR_col=gsl_matrix_column(scoreObsR, i);
    gsl_blas_dsyr(CblasLower, 1.0, &scoreObsR_col.vector, infMatR);
  }
  
  gsl_blas_dsyr(CblasLower, -1.0*NF, scoreObsBarF, infMatF);
  gsl_blas_dsyr(CblasLower, -1.0*NR, scoreObsBarR, infMatR);

  for(j=0;j<K;j++)
  {
    gsl_matrix_set(infMatPri,2*K-1+j,2*K-1+j,rho*(1./sigmaSqF[j]+1./sigmaSqR[j]));
    gsl_matrix_set(infMatPri,3*K-1+j,3*K-1+j,(alpha-1./2.)*sigmaSqF[j]*sigmaSqF[j]);
    gsl_matrix_set(infMatPri,4*K-1+j,4*K-1+j,(alpha-1./2.)*sigmaSqR[j]*sigmaSqR[j]);
  }

  for(j=0;j<K;j++)
  {
    gsl_matrix_set(infMatPri,2*K-1+j,3*K-1+j,(delta[j]-xi)*rho);
    gsl_matrix_set(infMatPri,2*K-1+j,4*K-1+j,(delta[j]-xi)*rho);
    gsl_matrix_set(infMatPri,3*K-1+j,2*K-1+j,(delta[j]-xi)*rho);
    gsl_matrix_set(infMatPri,4*K-1+j,2*K-1+j,(delta[j]-xi)*rho);

  }
    
    // Diagonal elements
  if(K>1)
  {
    gsl_matrix_set(infMatPri,K,K,-2*lambda);
    gsl_matrix_set(infMatPri,2*K-1,2*K-1,-2*lambda);
    
    for(j=1;j<(K-1);j++)
    {
      gsl_matrix_set(infMatPri,K+j,K+j,-4*lambda);
    }
      // Off diagonal ones
    for(j=0;j<(K-1);j++)
    {
      gsl_matrix_set(infMatPri,K+j+1,K+j,2*lambda);
      gsl_matrix_set(infMatPri,K+j,K+j+1,2*lambda);
    }
  }
    // I store everything in infMat
  gsl_matrix_add(infMat,infMatF);
  gsl_matrix_add(infMat,infMatR);
  gsl_matrix_add(infMat,infMatPri);

    //  gsl_matrix_fprintf (stdout, infMat, "%f");
  
  if(length(a)>0) // test if there is any missing reads
  {
	  int J=length(a);
	  double aNormF, aNormR, bNormF, bNormR, P0F=1, P0R=1, PhiF, PhiR, PhiFw, PhiRw, wPjF, wPjR; 
	  gsl_matrix *cdfAF=gsl_matrix_calloc(K,J), *cdfBF=gsl_matrix_calloc(K,J), *cdfAR=gsl_matrix_calloc(K,J),*cdfBR=gsl_matrix_calloc(K,J);
	  gsl_matrix *aSinF=gsl_matrix_calloc(K,J), *bSinF=gsl_matrix_calloc(K,J), *aSinR=gsl_matrix_calloc(K,J),*bSinR=gsl_matrix_calloc(K,J);
	  gsl_matrix *H3F=gsl_matrix_calloc(K,J), *H1F=gsl_matrix_calloc(K,J), *H2F=gsl_matrix_calloc(K,J),  *H0F=gsl_matrix_calloc(K,J);
	  gsl_matrix *H3R=gsl_matrix_calloc(K,J), *H1R=gsl_matrix_calloc(K,J), *H2R=gsl_matrix_calloc(K,J),  *H0R=gsl_matrix_calloc(K,J);
	  gsl_vector *PjF=gsl_vector_calloc(J), *PjR=gsl_vector_calloc(J);
	  //gsl_vector *OneJ=gsl_vector_calloc(J);
	  gsl_vector *NjF=gsl_vector_calloc(J), *NjR=gsl_vector_calloc(J);
	  gsl_matrix *qWmF=gsl_matrix_calloc(K,J),*qMumF=gsl_matrix_calloc(K,J),*qSigmamF=gsl_matrix_calloc(K,J) ;
	  gsl_matrix *qWmR=gsl_matrix_calloc(K,J),*qMumR=gsl_matrix_calloc(K,J),*qSigmamR=gsl_matrix_calloc(K,J) ;
	  gsl_matrix *scoreMissF=gsl_matrix_calloc(5*K-1,J),*scoreMissR=gsl_matrix_calloc(5*K-1,J);
	  gsl_vector *scoreMissBarF=gsl_vector_calloc(5*K-1), *scoreMissBarR=gsl_vector_calloc(5*K-1);
	  gsl_matrix *infMatMissF=gsl_matrix_calloc(5*K-1,5*K-1), *infMatMissR=gsl_matrix_calloc(5*K-1,5*K-1);
	  gsl_vector_view scoreMissF_col, scoreMissR_col;

	  //gsl_vector *rowSumH0F=gsl_vector_calloc(K), *rowSumH0R=gsl_vector_calloc(K);
	  //gsl_vector *rowSumH1F=gsl_vector_calloc(K), *rowSumH1R=gsl_vector_calloc(K);
  	  //gsl_vector *rowSumH2F=gsl_vector_calloc(K), *rowSumH2R=gsl_vector_calloc(K);
	  //gsl_vector *rowSumH3F=gsl_vector_calloc(K), *rowSumH3R=gsl_vector_calloc(K);
	  int *aP=INTEGER(a), *bP=INTEGER(b);
	  
	  /** Initialize the vector of Ones **/
	  //gsl_vector_set_all(OneJ, 1.0);
	  	  
	  
	  for(i=0;i<J;i++)
	  {
		  for(j=0;j<K;j++)
		  {
			  aNormF=(aP[i]-gsl_vector_get(muF,j))/sqrt(sigmaSqF[j]);
			  aNormR=(aP[i]-gsl_vector_get(muR,j))/sqrt(sigmaSqR[j]);
			  bNormF=(bP[i]-gsl_vector_get(muF,j))/sqrt(sigmaSqF[j]);
			  bNormR=(bP[i]-gsl_vector_get(muR,j))/sqrt(sigmaSqR[j]);
			  
			  gsl_matrix_set(cdfAF,j,i,gsl_cdf_tdist_P(aNormF,nu));
			  gsl_matrix_set(cdfAR,j,i,gsl_cdf_tdist_P(aNormR,nu));
			  gsl_matrix_set(cdfBF,j,i,gsl_cdf_tdist_P(bNormF,nu));
			  gsl_matrix_set(cdfBR,j,i,gsl_cdf_tdist_P(bNormR,nu));			  
			  
			  
			  gsl_matrix_set(aSinF,j,i,sin(atan(aNormF/2)));
			  gsl_matrix_set(aSinR,j,i,sin(atan(aNormR/2)));
			  gsl_matrix_set(bSinF,j,i,sin(atan(bNormF/2)));
			  gsl_matrix_set(bSinR,j,i,sin(atan(bNormR/2)));	
			  
			  
			  gsl_matrix_set(H3F,j,i,gsl_matrix_get(cdfBF,j,i)-gsl_matrix_get(cdfAF,j,i));
			  gsl_matrix_set(H3R,j,i,gsl_matrix_get(cdfBR,j,i)-gsl_matrix_get(cdfAR,j,i));
		  
			  gsl_matrix_set(H2F,j,i,cst*(fun2(gsl_matrix_get(bSinF,j,i)) - fun2(gsl_matrix_get(aSinF,j,i))));
			  gsl_matrix_set(H2R,j,i,cst*(fun2(gsl_matrix_get(bSinR,j,i)) - fun2(gsl_matrix_get(aSinR,j,i))));

			  gsl_matrix_set(H0F,j,i,cst*(fun0(gsl_matrix_get(bSinF,j,i)) - fun0(gsl_matrix_get(aSinF,j,i))));
			  gsl_matrix_set(H0R,j,i,cst*(fun0(gsl_matrix_get(bSinR,j,i)) - fun0(gsl_matrix_get(aSinR,j,i))));
 
			  gsl_matrix_set(H1F,j,i, (fun1(aNormF,4) - fun1(bNormF,4)) * cst/(nu+1));
			  gsl_matrix_set(H1R,j,i, (fun1(aNormR,4) - fun1(bNormR,4)) * cst/(nu+1));
			   
			  gsl_vector_set(PjF,i, gsl_vector_get(PjF,i)+w[j]*gsl_matrix_get(H3F,j,i));
			  gsl_vector_set(PjR,i, gsl_vector_get(PjR,i)+w[j]*gsl_matrix_get(H3R,j,i));
			  
		  }
	  }
	  			  
	  for(i=0;i<J;i++)
	  {
		  P0F -= gsl_vector_get(PjF,i);
		  P0R -= gsl_vector_get(PjR,i);
	  }
	  PhiF = NF/P0F;
	  PhiR = NR/P0R;
	  
	  for(i=0;i<J;i++)
	  {
		  gsl_vector_set(NjF,i, gsl_vector_get(PjF,i)*PhiF);
		  gsl_vector_set(NjR,i, gsl_vector_get(PjR,i)*PhiR);
	  }
	  
	  
	  //calculate the score matrix for missing data
	  for(i=0;i<J;i++)
	  {
		  for(j=0;j<K;j++)
		  {
			  wPjF = w[j] / gsl_vector_get(PjF,i);		  
			  //'q.w.m[k,j]=dQ_{k,mis}(y_n,theta)/d w_k', a matrix of 'K' by 'J'
			  gsl_matrix_set(qWmF,j,i, gsl_matrix_get(H3F,j,i)/gsl_vector_get(PjF,i));
			  //'qMu.m[k,j]=dQ_{k,mis}(y_n,theta)/d mu_k', a matrix of 'K' by 'J'
			  gsl_matrix_set(qMumF,j,i, wPjF*gsl_matrix_get(H1F,j,i)*2./sqrt(sigmaSqF[j]));
			  //'q.sigma.m[k,j]=dQ_{k,mis}(y_n,theta)/d sigma.sq_k', a matrix of 'K' by 'J'
			  gsl_matrix_set(qSigmamF,j,i, wPjF*sigmaSqR[j]/2.*(gsl_matrix_get(H0F,j,i)-nu*gsl_matrix_get(H2F,j,i)));
			  
			  wPjR = w[j] / gsl_vector_get(PjR,i);		  
			  //'q.w.m[k,j]=dQ_{k,mis}(y_n,theta)/d w_k', a matrix of 'K' by 'J'
			  gsl_matrix_set(qWmR,j,i, gsl_matrix_get(H3R,j,i)/gsl_vector_get(PjR,i));
			  //'qMu.m[k,j]=dQ_{k,mis}(y_n,theta)/d mu_k', a matrix of 'K' by 'J'
			  gsl_matrix_set(qMumR,j,i, wPjR*gsl_matrix_get(H1R,j,i)*2./sqrt(sigmaSqR[j]));
			  //'q.sigma.m[k,j]=dQ_{k,mis}(y_n,theta)/d sigma.sq_k', a matrix of 'K' by 'J'
			  gsl_matrix_set(qSigmamR,j,i, wPjR*sigmaSqR[j]/2.*(gsl_matrix_get(H0R,j,i)-nu*gsl_matrix_get(H2R,j,i)));
		  }
	  }
	  
	  	  
	  //the score matrix for miss data, a matrix of '5*K-1' by 'J'
	  for(i=0;i<J;i++)
	  {
		  for(j=1;j<K;j++)
		  {
			  gsl_matrix_set(scoreMissF, j-1, i, gsl_matrix_get(qWmF,j,i));
			  
			  gsl_matrix_set(scoreMissR, j-1, i, gsl_matrix_get(qWmR,j,i));
		  }
		  for(j=0;j<K;j++)
		  {
			  gsl_matrix_set(scoreMissF, (K-1)+j, i, gsl_matrix_get(qMumF,j,i));
			  gsl_matrix_set(scoreMissF, 2*K-1+j, i, -gsl_matrix_get(qMumF,j,i)/2.);
			  gsl_matrix_set(scoreMissF, 3*K-1+j, i, gsl_matrix_get(qSigmamF,j,i));
			  
			  gsl_matrix_set(scoreMissR, (K-1)+j, i, gsl_matrix_get(qMumR,j,i));
			  gsl_matrix_set(scoreMissR, 2*K-1+j, i, gsl_matrix_get(qMumR,j,i)/2.);
			  gsl_matrix_set(scoreMissR, 4*K-1+j, i, gsl_matrix_get(qSigmamR,j,i));
		  }
	  }
	  
	  // Compute scoreMiss*1 (sum over rows)
	  gsl_blas_dgemv(CblasNoTrans, 1.0/(1.0-P0F), scoreMissF, PjF, 0, scoreMissBarF);
	  gsl_blas_dgemv(CblasNoTrans, 1.0/(1.0-P0R), scoreMissR, PjR, 0, scoreMissBarR);
	  
	  // Scale the matrix to compute the mean
	  gsl_vector_scale(scoreMissBarF, 1./J);
	  gsl_vector_scale(scoreMissBarR, 1./J);
	  
	  //the information matrix, a matrix of '5*K-1' by '5*K-1'
	  for(i=0;i<J;i++)
	  {
		  scoreMissF_col=gsl_matrix_column(scoreMissF, i);
		  gsl_blas_dsyr(CblasLower, gsl_vector_get(NjF,i), &scoreMissF_col.vector, infMatMissF);

		  scoreMissR_col=gsl_matrix_column(scoreMissR, i);
		  gsl_blas_dsyr(CblasLower, gsl_vector_get(NjR,i), &scoreMissR_col.vector, infMatMissR);
	  }
	  
	  gsl_blas_dsyr(CblasLower, -1.0*NF*(1-P0F)/P0F, scoreMissBarF, infMatMissF);
	  gsl_blas_dsyr(CblasLower, -1.0*NR*(1-P0R)/P0R, scoreMissBarR, infMatMissR);
	  
	  // I store everything in infMat
	  gsl_matrix_add(infMat,infMatMissF);
	  gsl_matrix_add(infMat,infMatMissR);

	  
	  /** Free the memory **/ 
	  gsl_vector_free(PjF);gsl_vector_free(PjR); //gsl_vector_free(OneJ);
	  //gsl_vector_free(rowSumH0F);  gsl_vector_free(rowSumH1F); gsl_vector_free(rowSumH2F);    gsl_vector_free(rowSumH3F);
	  //gsl_vector_free(rowSumH0R);  gsl_vector_free(rowSumH1R); gsl_vector_free(rowSumH2R);    gsl_vector_free(rowSumH3R);
	  gsl_matrix_free(cdfAF); gsl_matrix_free(cdfBF); gsl_matrix_free(cdfAR); gsl_matrix_free(cdfBR);
	  gsl_matrix_free(aSinF); gsl_matrix_free(bSinF); gsl_matrix_free(aSinR); gsl_matrix_free(bSinR);
	  gsl_matrix_free(H0F); gsl_matrix_free(H1F);gsl_matrix_free(H2F);gsl_matrix_free(H3F);
	  gsl_matrix_free(H0R); gsl_matrix_free(H1R);gsl_matrix_free(H2R);gsl_matrix_free(H3R);
	  gsl_matrix_free(qWmF), gsl_matrix_free(qWmR);
	  gsl_matrix_free(qMumF), gsl_matrix_free(qMumR);
	  gsl_matrix_free(qSigmamF), gsl_matrix_free(qSigmamR);
	  gsl_matrix_free(scoreMissF),gsl_matrix_free(scoreMissR);
	  gsl_vector_free(scoreMissBarF), gsl_vector_free(scoreMissBarR);
	  gsl_matrix_free(infMatMissF), gsl_matrix_free(infMatMissR);	  
  }

  flag=gsl_linalg_cholesky_decomp(infMat);
    //There is an error so we simply return the flag
  if(flag==GSL_EDOM)
  {
      //printf("Cannot compute the cholesky decomposition\n");
      //Can't compute anything
    return(flag);
  }
  else
  {
    /* Compute the inverse of I^{-1/2} */
    /* Compute L'^{-1} */
    gsl_matrix_set_identity(DiagOne);
    flag=gsl_blas_dtrsm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, infMat, DiagOne);
    //There is an error so we simply return the flag
    if(flag!=0)
    {      
        //printf("Cannot solve with the cholesky decomposition\n");
        //Can't compute anything
      return(flag);
    }
    else
    {
      for(k=0;k<K;k++)
      {
        gsl_vector_set_zero(A);
        gsl_vector_set(A,K-1+k,1);
        /* Compute L'^{-1} A */
        flag=gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, DiagOne, A);
        if(flag!=0)
          return(flag);
        /* Compute the norm of that vector */
        gsl_vector_set(se,k,gsl_blas_dnrm2(A));
        // Make sure it's a number
        if(gsl_finite(gsl_vector_get(se,k))==0)
        {
          flag=1;
          return(flag);
        }
        gsl_vector_set_zero(A);
        gsl_vector_set(A,K-1+k,1);
        gsl_vector_set(A,K-1+K+k,-0.5);
        flag=gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, DiagOne, A);
        if(flag!=0)
          return(flag);
        gsl_vector_set(seF,k,gsl_blas_dnrm2(A));
        if(gsl_finite(gsl_vector_get(seF,k))==0)
        {
          flag=1;
          return(flag);
        }
        gsl_vector_set_zero(A);
        gsl_vector_set(A,K-1+k,1);
        gsl_vector_set(A,K-1+K+k,0.5);
        flag=gsl_blas_dtrmv (CblasUpper, CblasNoTrans, CblasNonUnit, DiagOne, A);
        if(flag!=0)
          return(flag);
        gsl_vector_set(seR,k,gsl_blas_dnrm2(A));
        if(gsl_finite(gsl_vector_get(seR,k))==0)
        {
          flag=1;
          return(flag);
        }
      }
    }    
  }
  
  // Copy infMat in DiagOne
  // So that I don't have to recompute the inverse of L'
  gsl_matrix_memcpy(infMat, DiagOne);
  
  gsl_vector_free(sumF), gsl_vector_free(sumR), gsl_vector_free(muF),gsl_vector_free(muR);
  gsl_matrix_free(rF), gsl_matrix_free(rR);
  gsl_matrix_free(ruyF), gsl_matrix_free(ruyR);
  gsl_vector_free(OneF), gsl_vector_free(OneR);
  gsl_matrix_free(qWoF), gsl_matrix_free(qWoR);
  gsl_matrix_free(qMuoF), gsl_matrix_free(qMuoR);
  gsl_matrix_free(qSigmaoF), gsl_matrix_free(qSigmaoR);
  gsl_matrix_free(scoreObsF),gsl_matrix_free(scoreObsR);
  gsl_vector_free(scoreObsBarF), gsl_vector_free(scoreObsBarR);
  gsl_matrix_free(infMatF), gsl_matrix_free(infMatR), gsl_matrix_free(infMatPri);
  gsl_vector_free(A);
  gsl_matrix_free(DiagOne);

  return(flag);
}


int mergePeak(SEXP para, gsl_matrix* infMat, gsl_vector* se, gsl_vector* seF, gsl_vector* seR, int *K, double nu, double nSe, double minSpacingPeaks, int dataType)
{
  int i=0,j=0,k=0,l=0,kMerge=0,flag=0;
  int K0=*K;
  gsl_matrix *Index=gsl_matrix_calloc(*K,*K);
  gsl_vector *A=gsl_vector_calloc(5**K-1),*B=gsl_vector_calloc(5**K-1),*C=gsl_vector_calloc(5**K-1);
  gsl_vector *OriginalW=gsl_vector_calloc(*K);
  bool loopflag=0;
  double DDelta, EisenbergerF, EisenbergerR, mycut=27/4;
  double minDiff=0,diff=0,SumCombW=0,SumW,TmpMu=0,TmpDelta=0,TmpSigmaSqF,TmpSigmaSqR;
  double *w=REAL(VECTOR_ELT(para, 0)), *mu=REAL(VECTOR_ELT(para, 1)), *delta=REAL(VECTOR_ELT(para, 2)), *sigmaSqF=REAL(VECTOR_ELT(para, 3)), *sigmaSqR=REAL(VECTOR_ELT(para, 4));

  gsl_matrix_set_identity(Index);
  
  for(k=0;k<*K;k++)
  {
    gsl_vector_set(OriginalW,k,w[k]);    
  }
  

  //calculate the loop flag	
  minDiff=minSpacingPeaks;
  if (dataType==1) 
  { //TF case, Find the maximum overlap for merging
	  for(k=0;k<*K-1;k++)
	  {    
		  diff=mu[k+1]-delta[k+1]/2.-nSe*gsl_vector_get(seF,k+1)-(mu[k]+delta[k]/2.+nSe*gsl_vector_get(seR,k));
		  if(diff<minDiff)
		  {
			  kMerge=k;
			  minDiff=diff;
		  }
	  }
  }	
  if (dataType==2) 
  {//histone case,  Find the centers of peaks with min distance, and Eisenberger distance less than cut off
	  for(k=0;k<*K-1;k++)
	  {		  
		  diff=mu[k+1]-mu[k];
		  DDelta=delta[k+1]-delta[k];
		  EisenbergerF = pow((diff-DDelta/2.),2.0) * (1./sigmaSqF[k+1]+ 1./sigmaSqF[k]);	
		  EisenbergerR = pow((diff+DDelta/2.),2.0) * (1./sigmaSqR[k+1]+ 1./sigmaSqR[k]);
		  if((diff<minDiff)&&(EisenbergerF<mycut)&&(EisenbergerR<mycut))
		  {
			  kMerge=k;
			  minDiff=diff;
		  }
	  }
  }
  
  loopflag=(minDiff<minSpacingPeaks);	

	
  while(loopflag)
  {
      // Compute the new parameters
    SumCombW=w[kMerge]+w[kMerge+1];
      // New mu
    TmpMu=(w[kMerge]*mu[kMerge]+w[kMerge+1]*mu[kMerge+1])/SumCombW;
      // New Delta
    TmpDelta=(w[kMerge]*delta[kMerge]+w[kMerge+1]*delta[kMerge+1])/SumCombW;
      // New SigmaSqF    
    TmpSigmaSqF=(nu-2.)/nu*(((nu/(nu-2.)*sigmaSqF[kMerge]+gsl_pow_2(mu[kMerge]-delta[kMerge]/2.))*w[kMerge]+(nu/(nu-2.)*sigmaSqF[kMerge+1]+gsl_pow_2(mu[kMerge+1]-delta[kMerge+1]/2.))*w[kMerge+1])/SumCombW-gsl_pow_2(TmpMu-TmpDelta/2.));
      // New SigmaSqR
    TmpSigmaSqR=(nu-2.)/nu*(((nu/(nu-2.)*sigmaSqR[kMerge]+gsl_pow_2(mu[kMerge]+delta[kMerge]/2.))*w[kMerge]+(nu/(nu-2.)*sigmaSqR[kMerge+1]+gsl_pow_2(mu[kMerge+1]+delta[kMerge+1]/2.))*w[kMerge+1])/SumCombW-gsl_pow_2(TmpMu+TmpDelta/2.));
      // Replace the current parameters by the new ones
    mu[kMerge]=TmpMu;
    delta[kMerge]=TmpDelta;
    w[kMerge]=SumCombW;
    sigmaSqF[kMerge]=TmpSigmaSqF;
    sigmaSqR[kMerge]=TmpSigmaSqR;
    
    
    SumW=0;
    for(k=0;k<K0;k++)
    {
      if((gsl_matrix_get(Index,kMerge,k)==1) | (gsl_matrix_get(Index,kMerge+1,k)==1))
      {
        gsl_matrix_set(Index,kMerge,k,1);
        SumW+=gsl_vector_get(OriginalW,k);
      }
    }

    for(k=0;k<K0;k++)
    {
      if(gsl_matrix_get(Index,kMerge,k)==1)
      {
        gsl_vector_set(A,K0-1+k,gsl_vector_get(OriginalW,k)/SumW);
        gsl_vector_set(B,K0-1+k,gsl_vector_get(OriginalW,k)/SumW);
        gsl_vector_set(B,K0-1+K0+k,-gsl_vector_get(OriginalW,k)/SumW/2.);
        gsl_vector_set(C,K0-1+k,gsl_vector_get(OriginalW,k)/SumW);
        gsl_vector_set(C,K0-1+K0+k,gsl_vector_get(OriginalW,k)/SumW/2.);
      }
    }
    

    flag=gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, infMat, A);
    if(flag!=0)
    {
      return(flag);
    }
    /* Compute the norm of that vector */
    gsl_vector_set(se,kMerge,gsl_blas_dnrm2(A));
    
    /* Reset the vector B */
    /* Computer T^-0.5 B */
    flag=gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, infMat, B);
    if(flag!=0)
    {
      return(flag);
    }
    
    /* Compute the norm of that vector */
    gsl_vector_set(seF,kMerge,gsl_blas_dnrm2(B));
    
    /* Reset the vector B */
    /* Computer T^-0.5 B */
    flag=gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, infMat, C);
    if(flag!=0)
    {
      return(flag);
    }
    
    /* Compute the norm of that vector */
    gsl_vector_set(seR,kMerge,gsl_blas_dnrm2(C));

    /* Reset the three vectors */
    gsl_vector_set_zero(A);
    gsl_vector_set_zero(B);
    gsl_vector_set_zero(C);

    /* Shift the components */
    for(k=(kMerge+1);k<*K-1;k++)
    {
      w[k]=w[k+1];
      mu[k]=mu[k+1];
      delta[k]=delta[k+1];
      sigmaSqF[k]=sigmaSqF[k+1];
      sigmaSqR[k]=sigmaSqR[k+1];
      
      gsl_vector_set(se,k,gsl_vector_get(se,k+1));
      gsl_vector_set(seF,k,gsl_vector_get(seF,k+1));
      gsl_vector_set(seR,k,gsl_vector_get(seR,k+1));         
      for(l=0;l<(*K);l++)
      {      
        gsl_matrix_set(Index,k,l,gsl_matrix_get(Index,k+1,l));
      }
    }

    /* Decrease the number of components by one */
    (*K)--;               

	  //calculate the loop flag	
	  minDiff=minSpacingPeaks;
	  if (dataType==1) 
	  { //TF case, Find the maximum overlap for merging
		  for(k=0;k<*K-1;k++)
		  {    
			  diff=mu[k+1]-delta[k+1]/2.-nSe*gsl_vector_get(seF,k+1)-(mu[k]+delta[k]/2.+nSe*gsl_vector_get(seR,k));
			  if(diff<minDiff)
			  {
				  kMerge=k;
				  minDiff=diff;
			  }
		  }
	  }	
	  if (dataType==2) 
	  {//histone case,  Find the centers of peaks with min distance, and Eisenberger distance less than cut off
		  for(k=0;k<*K-1;k++)
		  {		  
			  diff=mu[k+1]-mu[k];
			  DDelta=delta[k+1]-delta[k];
			  EisenbergerF = pow((diff-DDelta/2.),2.0) * (1./sigmaSqF[k+1]+ 1./sigmaSqF[k]);	
			  EisenbergerR = pow((diff+DDelta/2.),2.0) * (1./sigmaSqR[k+1]+ 1./sigmaSqR[k]);
			  if((diff<minDiff)&&(EisenbergerF<mycut)&&(EisenbergerR<mycut))
			  {
				  kMerge=k;
				  minDiff=diff;
			  }
		  }
	  }
	  loopflag=(minDiff<minSpacingPeaks);	  
   }

  gsl_matrix_free(Index);
  gsl_vector_free(A);
  gsl_vector_free(B);
  gsl_vector_free(C);
  gsl_vector_free(OriginalW);
    //Here i return 0, in the future we might want to use a flag to check errors
  return(flag);
}

