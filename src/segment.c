#include <R.h>
#include <Rmath.h>
#include <float.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>
//GSL stuff
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

SEXP getVector(SEXP list, SEXP ind);
SEXP getK(SEXP list);
SEXP getScore(SEXP list);
SEXP getChr(SEXP list);
SEXP getMap(SEXP list);
//PING functions
void wThreCountsPING(int *step, int *dataF, int *dataR, int *nReadsF, int *nReadsR, int *width, int *scoreF, int *scoreR);
void callRegionsLPING(int *center, int *nProbes, int *width, int *scoreF, int *scoreR, int *scoreRegionF, int *scoreRegionR, int *cutoff, int *StartRegion, int *EndRegion, int *nRegions, int maxStep, int kStep, int minL);
SEXP segRPING(SEXP chr, SEXP dataF, SEXP dataR, SEXP contF, SEXP contR, SEXP StartRegion, SEXP EndRegion, SEXP StartMap, SEXP EndMap, SEXP jitter, int nRegions);
//PICS functions
void wThreCounts(int *step, int *dataF, int *dataR, int *nReadsF, int *nReadsR, int *width, int *scoreF, int *scoreR);
void callRegions(int *center, int *lengthCenter, int *width, int *scoreF, int *scoreR, int *scoreRegionF, int *scoreRegionR, int *cutoff, int *StartRegion, int *EndRegion, int *nRegions);
void callRegionsL(int *center, int *nProbes, int *width, int *scoreF, int *scoreR, int *scoreRegionF, int *scoreRegionR, int *cutoff, int *StartRegion, int *EndRegion, int *nRegions, int maxStep, int kStep, int minL);
SEXP segReads(SEXP chr, SEXP dataF, SEXP dataR, SEXP contF, SEXP contR, SEXP StartMap, SEXP EndMap, SEXP jitter, SEXP width, SEXP cutoff, SEXP step, SEXP maxStep, SEXP minLength, SEXP pPackage);
SEXP segR(SEXP chr, SEXP dataF, SEXP dataR, SEXP contF, SEXP contR, SEXP StartRegion, SEXP EndRegion, SEXP StartMap, SEXP EndMap, SEXP jitter, int nRegions);

SEXP segReadsAll(SEXP data, SEXP dataC, SEXP StartMap, SEXP EndMap, SEXP jitter, SEXP paraSW, SEXP maxStep, SEXP minLength, SEXP pPackage);
SEXP getDensity(SEXP pics, SEXP strand, SEXP step, SEXP filter, SEXP sum, SEXP scale);
SEXP getDensityList(SEXP picsList, SEXP strand, SEXP step, SEXP filter, SEXP sum, SEXP scale);
SEXP getListElement(SEXP list, const char *str);
int testClass(SEXP list, int i);
int testObj(SEXP pPackage);

SEXP segReadsAll(SEXP data, SEXP dataC, SEXP StartMap, SEXP EndMap, SEXP jitter, SEXP paraSW, SEXP maxStep, SEXP minLength, SEXP pPackage)
{
  int i=0,nChr;
  SEXP d,dC;
  SEXP ans, chr;
  SEXP names;
  SEXP st,ed;
  SEXP contp,contm;
  
  //d=GET_SLOT(data,install("listData"));
  //dC=GET_SLOT(dataC,install("listData"));
  d=data;
  dC=dataC;
  nChr=length(d);
 
  //Rprintf("maxstep=%d\n", INTEGER_VALUE(maxStep));
  //Rprintf("step=%d\n", INTEGER_VALUE(VECTOR_ELT(paraSW,0)));
  //Rprintf("width=%d\n", INTEGER_VALUE(VECTOR_ELT(paraSW,1)));
  //Rprintf("minReads=%d\n", INTEGER_VALUE(VECTOR_ELT(paraSW,2)));
  

  PROTECT(names=getAttrib(d, R_NamesSymbol));
  PROTECT(ans=NEW_LIST(nChr));

  /*Rprintf("nChr=%d\n", nChr);*/
  for(i=0;i<nChr;i++)
  {
    chr=STRING_ELT(names, i);
    if(length(VECTOR_ELT(dC,i))>0)
    {
      contp=VECTOR_ELT(VECTOR_ELT(dC,i),0);
      contm=VECTOR_ELT(VECTOR_ELT(dC,i),1);
    }
    else
    {
      contp=R_NilValue;
      contm=R_NilValue;
    }
    if(length(StartMap)>0)
    {
      st=VECTOR_ELT(StartMap,i);
      ed=VECTOR_ELT(EndMap,i);
    }
    else
    {
      st=StartMap;
      ed=EndMap;
    }
	
    //Rprintf("process chr %s\n", mkChar(chr));
    SET_VECTOR_ELT(ans,i, segReads(chr, VECTOR_ELT(VECTOR_ELT(d,i),0), VECTOR_ELT(VECTOR_ELT(d,i),1), contp, contm, st, ed, jitter, VECTOR_ELT(paraSW,1), VECTOR_ELT(paraSW,2), VECTOR_ELT(paraSW,0),maxStep, minLength, pPackage));
    /*Rprintf("Finished chr %d \n",i);*/
  }
  UNPROTECT(2);
  return(ans);
}

SEXP segReads(SEXP chr, SEXP dataF, SEXP dataR, SEXP contF, SEXP contR, SEXP StartMap, SEXP EndMap, SEXP jitter, SEXP width, SEXP cutoff, SEXP step, SEXP maxStep, SEXP minLength, SEXP pPackage)
{
  //R//
  /*char * r_max = vmaxget(); //R//gets the state of the memory*/

  int *center;
  int *scoreF;
  int *scoreR;
  int *scoreRegionF;
  int *scoreRegionR;
  int nRegions;
  SEXP StartRegion;
  SEXP EndRegion;	

  int *dF=INTEGER(dataF), *dR=INTEGER(dataR);
  SEXP ans;
  /** Parameter used in the merging function **/
  int dMerge=2*INTEGER_VALUE(width);
  int i=0,nR=length(dataR),nF=length(dataF),count=0,lengthCenter=0;
  double m,M;

  int kStep= INTEGER_VALUE(width)/INTEGER_VALUE(step);
  //Rprintf("maxStep=%d,\t  kStep=%d\n", INTEGER_VALUE(maxStep) ,kStep);

  /** Sort the data **/
  R_isort(dF, nF);
  R_isort(dR, nR);
  if(length(contF)>0 & length(contR)>0)
  {
    R_isort(INTEGER(contF), length(contF));
    R_isort(INTEGER(contR), length(contR));
  }
  m=imin2(dF[0],dR[0]);
  M=imax2(dF[nF-1],dR[nR-1]);

  lengthCenter=(int)((M-m)/(INTEGER_VALUE(step)));
  //R//These 3 allocations could be Calloc
  /*center=(int*)Calloc(lengthCenter, int);*/
  /*scoreF=(int*)Calloc(lengthCenter, int);*/
  /*scoreR=(int*)Calloc(lengthCenter, int);*/
  center=(int*)R_alloc(lengthCenter, sizeof(int));
  scoreF=(int*)R_alloc(lengthCenter, sizeof(int));
  scoreR=(int*)R_alloc(lengthCenter, sizeof(int));

  /** Allocate memory for the start/end of each preprocessed region **/
  /** Because I do not know the size yet, I use the maximum possible size **/
  if(INTEGER_VALUE(width)<0)
  {
    Rprintf("width is negative (%d) and will cause memory allocation issues", INTEGER_VALUE(width));
  } 
  //R//These 2 allocations could be Calloc
  scoreRegionF=(int*)Calloc((M-m)/(2*INTEGER_VALUE(width)), int);
  scoreRegionR=(int*)Calloc((M-m)/(2*INTEGER_VALUE(width)), int);
  /*scoreRegionF=(int*)R_alloc((int)((M-m)/(2*INTEGER_VALUE(width))), sizeof(int));*/
  /*scoreRegionR=(int*)R_alloc((int)((M-m)/(2*INTEGER_VALUE(width))), sizeof(int));*/
  PROTECT(StartRegion=allocVector(INTSXP, lengthCenter));
  PROTECT(EndRegion=allocVector(INTSXP, lengthCenter));
	
	/* Create the vector of centers */
	for(i=0;i<lengthCenter;i++)
	{
		center[i]=m+i*INTEGER_VALUE(step);
	}

  if(testObj(pPackage))
  {
    wThreCounts(INTEGER(step), dF, dR, &nF, &nR, INTEGER(width), scoreF, scoreR);
  }
  else
  {
    wThreCountsPING(INTEGER(step), dF, dR, &nF, &nR, INTEGER(width), scoreF, scoreR);
  }
  if (INTEGER_VALUE(maxStep)>0) 
  {
	  //Rprintf("Call Long region\n");
    if(testObj(pPackage))
    {
      callRegionsL(center, &lengthCenter, &dMerge, scoreF, scoreR, scoreRegionF, scoreRegionR, INTEGER(cutoff), INTEGER(StartRegion), INTEGER(EndRegion), &nRegions, INTEGER_VALUE(maxStep), kStep, INTEGER_VALUE(minLength));
    }
    else
    {
      callRegionsLPING(center, &lengthCenter, &dMerge, scoreF, scoreR, scoreRegionF, scoreRegionR, INTEGER(cutoff), INTEGER(StartRegion), INTEGER(EndRegion), &nRegions, INTEGER_VALUE(maxStep), kStep, INTEGER_VALUE(minLength));
    }
  }else 
  {
	  	  //Rprintf("Call region\n");
	   callRegions(center, &lengthCenter, &dMerge, scoreF, scoreR, scoreRegionF, scoreRegionR, INTEGER(cutoff), INTEGER(StartRegion), INTEGER(EndRegion), &nRegions);
  }
	/*
	*/
	
  if(nRegions>0)
  {
    if(testObj(pPackage))
    {
      PROTECT(ans=segR(chr, dataF, dataR, contF, contR, StartRegion, EndRegion, StartMap, EndMap, jitter, nRegions));
    }
    else
    {
      PROTECT(ans=segRPING(chr, dataF, dataR, contF, contR, StartRegion, EndRegion, StartMap, EndMap, jitter, nRegions));
    }
  }
  else
  {  
    PROTECT(ans=R_NilValue);
  }
	
  UNPROTECT(3);
  /*vmaxset(r_max); //R//Returns to the saved memory state (or smth like that..)*/
  /*Free(center);*/
  /*Free(scoreF);*/
  /*Free(scoreR);*/
  /*Free(scoreRegionF);*/
  /*Free(scoreRegionR);*/
  return(ans);
}

void wThreCounts(int *step, int* dataF, int* dataR, int *nReadsF, int *nReadsR, int *width, int *scoreF, int *scoreR)
{
  int m=imin2(dataR[0],dataF[0]);
  int M=imax2(dataR[*nReadsR-1],dataF[*nReadsF-1]);
  int center=m;
  int nbR=0,nbF=0;
  int startF=0,startR=0;
  int nbCenters=0,i;

  while(center<M)
  {    
    /* Count the number of forward reads within width of center on the left side only */    
    nbF=0;    
    i=startF;
    while((i<*nReadsF) && (center-dataF[i])>*width)
    {
      i++;
    }
    startF=i;
    while((i<*nReadsF) && ((center-dataF[i])<=*width) && ((center-dataF[i])>=0))
    {      
      nbF++;
      i++;
    }

    /* Count the number of reverse reads within width of center on the right side only */        
    nbR=0;
    i=startR;
    while((i<*nReadsR) && (dataR[i]-center)<0)
    {
      i++;
    }
    startR=i;
    while((i<*nReadsR) && ((dataR[i]-center)<=*width) && ((dataR[i]-center)>=0))
    {      
      nbR++;
      i++;
    }

    scoreF[nbCenters]=nbF;
		scoreR[nbCenters]=nbR;
    nbCenters++;
    center+=*step;
  }  
}

void wThreCountsPING(int *step, int* dataF, int* dataR, int *nReadsF, int *nReadsR, int *width, int *scoreF, int *scoreR)
{
  int m=imin2(dataR[0],dataF[0]);
  int M=imax2(dataR[*nReadsR-1],dataF[*nReadsF-1]);
  int center=m;
  int nbR=0,nbF=0;
  int startF=0,startR=0;
  int nbCenters=0,i;
	int minDist=25; // the min-distance, when reads are too close to the center of window, they are not counted in the window score

  while(center<M)
  {    
    /* Count the number of forward reads within width of center on the left side only */    
    nbF=0;    
    i=startF;
    while((i<*nReadsF) && (center-dataF[i])>*width)
    {
      i++;
    }
    startF=i;
    while((i<*nReadsF) && ((center-dataF[i])<=*width) && ((center-dataF[i])>=minDist))
    {      
      nbF++;
      i++;
    }

    /* Count the number of reverse reads within width of center on the right side only */        
    nbR=0;
    i=startR;
    while((i<*nReadsR) && (dataR[i]-center)<minDist)
    {
      i++;
    }
    startR=i;
    while((i<*nReadsR) && ((dataR[i]-center)<=*width) && ((dataR[i]-center)>=minDist))
    {      
      nbR++;
      i++;
    }

    scoreF[nbCenters]=nbF;
		scoreR[nbCenters]=nbR;
    nbCenters++;
    center+=*step;
  }  
}

void callRegions(int *center, int *nProbes, int *width, int *scoreF, int *scoreR, int *scoreRegionF, int *scoreRegionR, int *cutoff, int *StartRegion, int *EndRegion, int *nRegions)
{
  int i=0,p=0;
  int w=0,ww=0,max=0,maxF=0,maxR=0;
  *nRegions=0;

  while(p<*nProbes)
  {
    if((scoreF[p]>=*cutoff) & (scoreR[p]>=*cutoff))
    {
     maxF=scoreF[p];
     maxR=scoreR[p];     
      (*nRegions)++;
      StartRegion[*nRegions-1]=center[p]-*width/2;
      w=p+1;
      max=p;
      ww=p;

      while((w<*nProbes) && ((center[w]-center[ww])<=*width))
      {
        if((scoreF[w]>=*cutoff) & (scoreR[w]>=*cutoff))
        {
         if(scoreF[w]>maxF)
         {
           maxF=scoreF[w];
         }
         
         if(scoreR[w]>maxR)
         {
           maxR=scoreR[w];
         }
          max=w;
          ww=max;
        }
        w++;
      }
     scoreRegionF[*nRegions-1]=maxF;
     scoreRegionR[*nRegions-1]=maxR;
     EndRegion[*nRegions-1]=center[max]+*width/2;
     p=w;
    }
    else
    {
      p++;
    }
  }
}

void callRegionsL(int *center, int *nProbes, int *width, int *scoreF, int *scoreR, int *scoreRegionF, int *scoreRegionR, int *cutoff, int *StartRegion, int *EndRegion, int *nRegions, int maxStep, int kStep, int minL)
{
	//Rprintf("\n Regions shorter than %i bps, will not be called. \n", minL);
	int p=0, w=0, ww=0, start=0, minScore=0,  min=0,  cutCenter=0, pp; 
	//p		:the index of current center (outer loop)
	//ww	:the index of temperory last center in the regions
	//w		:the index of current center (inner loop)
	//start	:the index of first center in current segment
	//minScore: save the minimium score of current region
	//min	:save the index of center correponding to the min score, used for cut window in ceter
	//cutCenter: a flag indicating if last region was cut inthe center of window
	//kStep	:INTEGER(width)/INTEGER(step), when cut in the center, I want the index of start window is cut point+kStep
	
	*nRegions=0;
	
	//Rprintf("nProbes=%i \n",*nProbes);
	while(p < *nProbes)
	{
		//Rprintf("Start loop: p=%i \n",p);
		if(((scoreF[p]>=*cutoff) && (scoreR[p]>=*cutoff))|| (cutCenter==1)) //current window has enough reads to be recorded
		{
			//maxF=scoreF[p];
			//maxR=scoreR[p];     
			(*nRegions)++;
			
			/*set start point of the region according to if last region was cut in the center*/
			if (cutCenter==0) 
			{
				StartRegion[*nRegions-1]=center[p]-*width/2;
				//update new minscore
				start=p;
				min=start;
 				minScore=imin2(scoreF[start], scoreR[start]);
				//Rprintf("\n Start a new segment \n");
			}else 
			{
				StartRegion[*nRegions-1]=EndRegion[*nRegions-2]+1;
				//update new mini score
				start=min+kStep;
				min=start;
				minScore=imin2(scoreF[start], scoreR[start]);
				for (pp=start; pp<=p; pp++) {
					if (scoreF[pp]<minScore)
					{
						minScore= scoreF[pp];
						min=pp;
					}
					if (scoreR[pp]<minScore) 
					{
						minScore= scoreR[pp];
						min=pp;
					}
					pp++;
					
				}
				//Rprintf("\n Continue from last cut window \n");
			}
							
			//Rprintf("StartRegion[%i]=%i \n", *nRegions-1,StartRegion[*nRegions-1]);
			//Rprintf("min=%i, \t start=%i, \t p=%i, \t kStep=%i \n", min, start, p, kStep);

			
			w=p+1;  
			ww=p;	

			
			while(((w-start)<=maxStep) && ((center[w]-center[ww])<=*width) && (w < *nProbes))
			{
				if((scoreF[w]>=*cutoff) & (scoreR[w]>=*cutoff))
				{
					ww=w;
					
					if (scoreF[w]<minScore)
					{
						minScore= scoreF[w];
						min=w;
					}
					if (scoreR[w]<minScore) 
					{
						minScore= scoreR[w];
						min=w;
					}

				}
				w++;
			}
			
			if (w == *nProbes)
			{
				EndRegion[*nRegions-1]=center[ww]+*width/2;
			}else if ((w-start)>maxStep) 
			{
				EndRegion[*nRegions-1]=center[min];
				cutCenter=1;
			}else {
				EndRegion[*nRegions-1]=center[ww]+*width/2;
				cutCenter=0;
			}
			//Rprintf("EndRegion[%i]=%i \n", *nRegions-1,EndRegion[*nRegions-1]);
			if ((EndRegion[*nRegions-1]-StartRegion[*nRegions-1]) < minL) 
			{
				//Rprintf("\n %i - %i < %i, go back one step \n", EndRegion[*nRegions-1], StartRegion[*nRegions-1], minL);
				(*nRegions)--;
			}


			p=w;
		} else // do nothing, look at next window
		{
			p++;
		}/*else if (cutCenter==1) //current window have too few reads, but last regions was cut
		  { 
		  (*nRegions)++;
		  StartRegion[*nRegions-1]=EndRegion[*nRegions-2]+1;
		  EndRegion[*nRegions-1]=center[p-1]+*width/2;
		  cutCenter=0;
		  Rprintf("StartRegion[%i]=%i \n", *nRegions-1,StartRegion[*nRegions-1]);
		  Rprintf("EndRegion[%i]=%i \n", *nRegions-1,EndRegion[*nRegions-1]);
		  p=p-1+kStep;
		  if ((EndRegion[*nRegions-1]-StartRegion[*nRegions-1]) < minL) 
		  {
		  Rprintf("\n %i - %i < %i, go back one step \n", EndRegion[*nRegions-1], StartRegion[*nRegions-1], minL);
		  (*nRegions)--;
		  }
		  Rprintf("\n End last cut window \n");
		  }*/
		//Rprintf("end loop: p=%i \n",p);
	}
}

void callRegionsLPING(int *center, int *nProbes, int *width, int *scoreF, int *scoreR, int *scoreRegionF, int *scoreRegionR, int *cutoff, int *StartRegion, int *EndRegion, int *nRegions, int maxStep, int kStep, int minL)
{
	//Rprintf("\n Regions shorter than %i bps, will not be called. \n", minL);
	int p=0, w=0, ww=0, start=0, minScore=0,  min=0,  cutCenter=0, pp; 
	//p		:the index of current center (outer loop)
	//ww	:the index of temperory last center in the regions
	//w		:the index of current center (inner loop)
	//start	:the index of first center in current segment
	//minScore: save the minimium score of current region
	//min	:save the index of center correponding to the min score, used for cut window in ceter
	//cutCenter: a flag indicating if last region was cut inthe center of window
	//kStep	:INTEGER(width)/INTEGER(step), when cut in the center, I want the index of start window is cut point+kStep
	
	*nRegions=0;
	
	//Rprintf("nProbes=%i \n",*nProbes);
	while(p < *nProbes)
	{
		//Rprintf("Start loop: p=%i \n",p);
		if(((scoreF[p]>=*cutoff) && (scoreR[p]>=*cutoff))|| (cutCenter==1)) //current window has enough reads to be recorded
		{
			//maxF=scoreF[p];
			//maxR=scoreR[p];     
			(*nRegions)++;
			
			/*set start point of the region according to if last region was cut in the center*/
			if (cutCenter==0) 
			{
				StartRegion[*nRegions-1]=center[p]-*width/2;
				//update new minscore
				start=p;
				min=start;
 				minScore=imin2(scoreF[start], scoreR[start]);
				//Rprintf("\n Start a new segment \n");
			}else 
			{
				StartRegion[*nRegions-1]=EndRegion[*nRegions-2]+1;
				//update new mini score
				start=min+kStep;
				min=start;
				minScore=imin2(scoreF[start], scoreR[start]);
				for (pp=start; pp<=p; pp++) {
					if (scoreF[pp]<minScore)
					{
						minScore= scoreF[pp];
						min=pp;
					}
					if (scoreR[pp]<minScore) 
					{
						minScore= scoreR[pp];
						min=pp;
					}
					pp++;
					
				}
				//Rprintf("\n Continue from last cut window \n");
			}
							
			//Rprintf("StartRegion[%i]=%i \n", *nRegions-1,StartRegion[*nRegions-1]);
			//Rprintf("min=%i, \t start=%i, \t p=%i, \t kStep=%i \n", min, start, p, kStep);

			
			w=p+1;  
			ww=p;	

			
			while(((w-start)<=maxStep) && ((center[w]-center[ww])<=*width) && (w < *nProbes))
			{
				if((scoreF[w]>=*cutoff) & (scoreR[w]>=*cutoff))
				{
					ww=w;
					
					if (scoreF[w]<minScore)
					{
						minScore= scoreF[w];
						min=w;
					}
					if (scoreR[w]<minScore) 
					{
						minScore= scoreR[w];
						min=w;
					}

				}
				w++;
			}
			
			if (w == *nProbes)
			{
				EndRegion[*nRegions-1]=center[ww]+*width/2;
			}else if ((ww-start)>=maxStep) 
			{
				EndRegion[*nRegions-1]=center[min];
				cutCenter=1;
			}else {
				EndRegion[*nRegions-1]=center[ww]+*width/2;
				cutCenter=0;
			}
			//Rprintf("EndRegion[%i]=%i \n", *nRegions-1,EndRegion[*nRegions-1]);
			if ((EndRegion[*nRegions-1]-StartRegion[*nRegions-1]) < minL) 
			{
				//Rprintf("\n %i - %i < %i, go back one step \n", EndRegion[*nRegions-1], StartRegion[*nRegions-1], minL);
				(*nRegions)--;
			}


			p=w;
		} else // do nothing, look at next window
		{
			p++;
		}/*else if (cutCenter==1) //current window have too few reads, but last regions was cut
		  { 
		  (*nRegions)++;
		  StartRegion[*nRegions-1]=EndRegion[*nRegions-2]+1;
		  EndRegion[*nRegions-1]=center[p-1]+*width/2;
		  cutCenter=0;
		  Rprintf("StartRegion[%i]=%i \n", *nRegions-1,StartRegion[*nRegions-1]);
		  Rprintf("EndRegion[%i]=%i \n", *nRegions-1,EndRegion[*nRegions-1]);
		  p=p-1+kStep;
		  if ((EndRegion[*nRegions-1]-StartRegion[*nRegions-1]) < minL) 
		  {
		  Rprintf("\n %i - %i < %i, go back one step \n", EndRegion[*nRegions-1], StartRegion[*nRegions-1], minL);
		  (*nRegions)--;
		  }
		  Rprintf("\n End last cut window \n");
		  }*/
		//Rprintf("end loop: p=%i \n",p);
	}
}

SEXP segR(SEXP chr, SEXP dataF, SEXP dataR, SEXP contF, SEXP contR, SEXP StartRegion, SEXP EndRegion, SEXP StartMap, SEXP EndMap, SEXP jitter, int nRegions)
{
	int nF = length(dataF), nR = length(dataR), ncF = length(contF), ncR = length(contR), nMap = length(StartMap);
	int i=0, j=0, pF=0, pR=0, pcF=0, pcR=0, pM=0;
	int indStart=0, indEnd=0, indStartM=0, indEndM=0, indStartF=0, indEndF=0,indStartR=0, indEndR=0, *rmap, nProtect;
	int minLoc, maxLoc, temp; //temparoryly save boundary each regions
	SEXP ans, map, seg, yF, yR, cF, cR, classDef;
	SEXP name;
	
	PROTECT(name=NEW_CHARACTER(1));
	SET_STRING_ELT(name,0,mkChar(CHAR(chr)));
	GetRNGstate();
	/** Define the list for the output **/
	PROTECT(ans = NEW_LIST(nRegions));
	for(i=0;i<nRegions;i++)
	{
		/** Initialize the protects **/
		nProtect=0;
		
		if (pM>0) pM--; //if PM>0 we want to move the counter one step back since last unmappable regions might overlap with two segments
		/*Process mappability profile*/
		if((pM<nMap) & (nMap>0))
		{
			while((pM<nMap) & (INTEGER(EndMap)[pM]<INTEGER(StartRegion)[i]))
			{
				pM++;
			}
			/** Makes sure the index does not go out of bound */
			indStartM=imin2(pM,nMap);
			minLoc=imax2(INTEGER(StartMap)[indStartM],INTEGER(StartRegion)[i]);
			
			/** Keep looking **/
			while((pM<nMap) & (INTEGER(StartMap)[pM]<INTEGER(EndRegion)[i]))
			{
				pM++;
			}
			indEndM=imin2(pM,nMap);
			maxLoc=imin2(INTEGER(EndMap)[indEndM],INTEGER(EndRegion)[i]);
			
			
			PROTECT(map = allocMatrix(INTSXP,indEndM-indStartM,2));
			nProtect++;
			rmap=INTEGER(map);
			
			for(j=indStartM;j<indEndM;j++)
			{
				rmap[j-indStartM+(indEndM-indStartM)*0]=imax2(INTEGER(StartMap)[j],INTEGER(StartRegion)[i]);
				rmap[j-indStartM+(indEndM-indStartM)*1]=imin2(INTEGER(EndMap)[j],INTEGER(EndRegion)[i]);
			}
			
			/*
			Rprintf("IndStartM=%i,\t IndEndM=%i, \t nmap=%i \n",indStartM,indEndM, nMap );
			
			Rprintf("StartMap: \t");
			for(j=indStartM;j<indEndM;j++)
			{
				Rprintf("%i\t",INTEGER(map)[j-indStartM+(indEndM-indStartM)*0]);
			}
			Rprintf("\n");
			Rprintf("EndMap: \t");
			for(j=indStartM;j<indEndM;j++)
			{
				Rprintf("%i\t",INTEGER(map)[j-indStartM+(indEndM-indStartM)*1]);
			}
			Rprintf("\n");
			 */
			
		}
		else
		{
			minLoc=INTEGER(StartRegion)[i];
			maxLoc=INTEGER(EndRegion)[i];
			PROTECT(map = allocMatrix(INTSXP,0,2));
			nProtect++;
		}
		
		//Rprintf("No Truncation: \t\t minLoc[%i]=%i,\t maxLoc[%i]=%i \n",i,INTEGER(StartRegion)[i],i,INTEGER(EndRegion)[i]);
		//Rprintf("Map Truncation: \t minLoc[%i]=%i,\t maxLoc[%i]=%i \n",i,minLoc,i,maxLoc);
		
		/* find the index of 1st yF bounded by regions start,
		 and define "minLoc=min(yF[1], StartMap[1])" */
		while((pF<nF) && (INTEGER(dataF)[pF]<INTEGER(StartRegion)[i]))
		{
			pF++;
		}
		indStartF=pF;
		temp=minLoc;
		if (temp>INTEGER(StartRegion)[i]) 
		{
			minLoc=imin2(temp, INTEGER(dataF)[indStartF]);
		}else {
			minLoc=imax2(temp, INTEGER(dataF)[indStartF]);
		}

		
	
		/* find the index of 1st yR bounded by minLoc */
		while((pR<nR) && (INTEGER(dataR)[pR]<minLoc))
		{
			pR++;
		}
		indStartR=pR;
		
		
		/* find the index of last yR bounded by regions ends,
		 and define "maxLoc=min(yF[max], EndMap[max])" */
		while((pR<nR) && (INTEGER(dataR)[pR]<=INTEGER(EndRegion)[i]))
		{
			pR++;
		}
		indEndR=imin2(pR-1,nR);
		indEndR=imax2(indEndR,indStartR);
		temp=maxLoc;
		if (temp<INTEGER(EndRegion)[i])
		{
			maxLoc=imax2(temp, INTEGER(dataR)[indEndR]);
		}else {
			maxLoc=imin2(temp, INTEGER(dataR)[indEndR]);
		}

		//Rprintf("Reads Truncation: \t minLoc[%i]=%i,\t minLoc[%i]=%i \n",i,minLoc,i,maxLoc);
		
		/* find the index of last yF bounded by maxLoc */
		while((pF<nF) && (INTEGER(dataF)[pF]<=maxLoc))
		{
			pF++;
		}
		indEndF=imin2(pF-1,nF);
		indEndF=imax2(indEndF,indStartF);
		
		
		/*
		 Rprintf("Start: yF[%i]=%i, \n", indStartF,  INTEGER(dataF)[indStartF]);
		 Rprintf("Start: yR[%i]=%i, \n", indStartR,  INTEGER(dataR)[indStartR]);
		 Rprintf("End: yF[%i]=%i,   \n", indEndF,    INTEGER(dataF)[indEndF]);
		 Rprintf("End: yR[%i]=%i,   \n", indEndR,    INTEGER(dataR)[indEndR]);		
		 */
		
		/** Split the data using the start/end index **/
		PROTECT(yF = allocVector(REALSXP,indEndF-indStartF+1));
		nProtect++;
		for(j=indStartF;j<=indEndF;j++)
		{
			REAL(yF)[j-indStartF]=INTEGER(dataF)[j]+rnorm(0,.1)*INTEGER(jitter)[0];
		}
		
		PROTECT(yR = allocVector(REALSXP,indEndR-indStartR+1));
		nProtect++;    
		for(j=indStartR;j<=indEndR;j++)
		{
			REAL(yR)[j-indStartR]=INTEGER(dataR)[j]+rnorm(0,.1)*INTEGER(jitter)[0];
		}
		
		
		/*Process control data*/
		if((ncF>0) & (ncR>0))
		{
			while((pcF<ncF) & (INTEGER(contF)[pcF]<minLoc))
			{
				pcF++;
			}
			indStart=pcF;
			
			while((pcF<ncF) & (INTEGER(contF)[pcF]<=maxLoc))
			{
				pcF++;
			}
			indEnd=imin2(pcF,ncF);
			
			PROTECT(cF = allocVector(REALSXP,indEnd-indStart));
			nProtect++;      
			for(j=indStart;j<indEnd;j++)
			{
				REAL(cF)[j-indStart]=INTEGER(contF)[j]+rnorm(0,.1)*INTEGER(jitter)[0];
			}
			
			while((pcR<ncR) & (INTEGER(contR)[pcR]<minLoc))
			{
				pcR++;
			}
			indStart=pcR;
			
			while((pcR<ncR) & (INTEGER(contR)[pcR]<=maxLoc))
			{
				pcR++;
			}
			indEnd=imin2(pcR,ncR);
			
			PROTECT(cR = allocVector(REALSXP,indEnd-indStart));
			nProtect++;      
			for(j=indStart;j<indEnd;j++)
			{
				REAL(cR)[j-indStart]=INTEGER(contR)[j]+rnorm(0,.1)*INTEGER(jitter)[0];
			}
		}
		else
		{
			cF = R_NilValue;
			cR = R_NilValue;
		}
		
		classDef=MAKE_CLASS("segReads");
		PROTECT(seg=NEW_OBJECT(classDef));
		nProtect++;    
		SET_SLOT(seg,mkChar("yF"),yF);
		SET_SLOT(seg,mkChar("yR"),yR);
		SET_SLOT(seg,mkChar("cF"),cF);
		SET_SLOT(seg,mkChar("cR"),cR);
		SET_SLOT(seg,mkChar("map"),map);
		SET_SLOT(seg,mkChar("chr"),name);
		SET_VECTOR_ELT(ans,i,seg);
		UNPROTECT(nProtect);
	}
	UNPROTECT(2);
	PutRNGstate();
	return(ans);
}

SEXP segRPING(SEXP chr, SEXP dataF, SEXP dataR, SEXP contF, SEXP contR, SEXP StartRegion, SEXP EndRegion, SEXP StartMap, SEXP EndMap, SEXP jitter, int nRegions)
{
	int nF = length(dataF), nR = length(dataR), ncF = length(contF), ncR = length(contR), nMap = length(StartMap);
	int i=0, j=0, pF=0, pR=0, pcF=0, pcR=0, pM=0;
	int indStart=0, indEnd=0, indStartM=0, indEndM=0, indStartF=0, indEndF=0,indStartR=0, indEndR=0, *rmap, nProtect;
	int minLocF, maxLocF, minLocR, maxLocR, temp; //temparoryly save boundary each regions
	SEXP ans, map, seg, yF, yR, cF, cR, classDef;
	SEXP name;

	int ext=50; // extend each segments by ext bps on each side, only forward/reverse reads exist in the extended regions on the left/right
	
	PROTECT(name=NEW_CHARACTER(1));
	SET_STRING_ELT(name,0,mkChar(CHAR(chr)));
	GetRNGstate();
	/** Define the list for the output **/
	PROTECT(ans = NEW_LIST(nRegions));
	for(i=0;i<nRegions;i++)
	{
		/** Initialize the protects **/
		nProtect=0;
		
		if (pM>0) pM--; //if PM>0 we want to move the counter one step back since last unmappable regions might overlap with two segments
		/*Process mappability profile*/
		if((pM<nMap) & (nMap>0))
		{
			while((pM<nMap) & (INTEGER(EndMap)[pM]<INTEGER(StartRegion)[i]))
			{
				pM++;
			}
			/** Makes sure the index does not go out of bound */
			indStartM=imin2(pM,nMap);
			minLocF=imax2(INTEGER(StartMap)[indStartM],INTEGER(StartRegion)[i]-ext);
			//minLocR=imax2(INTEGER(StartMap)[indStartM],INTEGER(StartRegion)[i]+ext);
			minLocR=INTEGER(StartRegion)[i]+ext;
			
			
			/** Keep looking **/
			while((pM<nMap) & (INTEGER(StartMap)[pM]<INTEGER(EndRegion)[i]))
			{
				pM++;
			}
			indEndM=imin2(pM,nMap);
			//maxLocF=imin2(INTEGER(EndMap)[indEndM],INTEGER(EndRegion)[i]-ext);
			maxLocR=imin2(INTEGER(EndMap)[indEndM],INTEGER(EndRegion)[i]+ext);
			maxLocF=INTEGER(EndRegion)[i]-ext;
			
			
			PROTECT(map = allocMatrix(INTSXP,indEndM-indStartM,2));
			nProtect++;
			rmap=INTEGER(map);
			
			for(j=indStartM;j<indEndM;j++)
			{
				rmap[j-indStartM+(indEndM-indStartM)*0]=imax2(INTEGER(StartMap)[j],INTEGER(StartRegion)[i]);
				rmap[j-indStartM+(indEndM-indStartM)*1]=imin2(INTEGER(EndMap)[j],INTEGER(EndRegion)[i]);
			}
			
			
		}
		else
		{
			minLocF=INTEGER(StartRegion)[i]-ext;
			minLocR=INTEGER(StartRegion)[i]+ext;
			maxLocF=INTEGER(EndRegion)[i]-ext;
			maxLocR=INTEGER(EndRegion)[i]+ext;
			PROTECT(map = allocMatrix(INTSXP,0,2));
			nProtect++;
		}
		
		//Rprintf("No Truncation: \t\t minLoc[%i]=%i,\t maxLoc[%i]=%i \n",i,INTEGER(StartRegion)[i],i,INTEGER(EndRegion)[i]);
		//Rprintf("Map Truncation: \t minLoc[%i]=%i,\t maxLoc[%i]=%i \n",i,minLoc,i,maxLoc);
		
		/* find the index of 1st yF bounded by regions start,
		 and define "minLoc=min(yF[1], StartMap[1])" */
		while((pF<nF) && (INTEGER(dataF)[pF]<((INTEGER(StartRegion)[i])-ext)))
		{
			pF++;
		}
		indStartF=pF;
		temp=minLocF;
		if (temp>(INTEGER(StartRegion)[i]-ext)) 
		{
			minLocF=imin2(temp, INTEGER(dataF)[indStartF]);
		}else {
			minLocF=imax2(temp, INTEGER(dataF)[indStartF]);
		}

		
	
		/* find the index of 1st yR bounded by minLoc */
		while((pR<nR) && (INTEGER(dataR)[pR]<minLocR))
		{
			pR++;
		}
		indStartR=pR;
		
		
		/* find the index of last yR bounded by regions ends,
		 and define "maxLoc=min(yF[max], EndMap[max])" */
		while((pR<nR) && (INTEGER(dataR)[pR]<=(INTEGER(EndRegion)[i]+ext)))
		{
			pR++;
		}
		indEndR=imin2(pR-1,nR);
		indEndR=imax2(indEndR,indStartR);
		temp=maxLocR;
		if (temp<INTEGER(EndRegion)[i])
		{
			maxLocR=imax2(temp, INTEGER(dataR)[indEndR]);
		}else {
			maxLocR=imin2(temp, INTEGER(dataR)[indEndR]);
		}

		//Rprintf("Reads Truncation: \t minLoc[%i]=%i,\t minLoc[%i]=%i \n",i,minLoc,i,maxLoc);
		
		/* find the index of last yF bounded by maxLoc */
		while((pF<nF) && (INTEGER(dataF)[pF]<=maxLocF))
		{
			pF++;
		}
		indEndF=imin2(pF-1,nF);
		indEndF=imax2(indEndF,indStartF);
		
		
		/*
		 Rprintf("Start: yF[%i]=%i, \n", indStartF,  INTEGER(dataF)[indStartF]);
		 Rprintf("Start: yR[%i]=%i, \n", indStartR,  INTEGER(dataR)[indStartR]);
		 Rprintf("End: yF[%i]=%i,   \n", indEndF,    INTEGER(dataF)[indEndF]);
		 Rprintf("End: yR[%i]=%i,   \n", indEndR,    INTEGER(dataR)[indEndR]);		
		 */
		
		/** Split the data using the start/end index **/
		PROTECT(yF = allocVector(REALSXP,indEndF-indStartF+1));
		nProtect++;
		for(j=indStartF;j<=indEndF;j++)
		{
			REAL(yF)[j-indStartF]=INTEGER(dataF)[j]+rnorm(0,2)*INTEGER(jitter)[0];
		}
		
		PROTECT(yR = allocVector(REALSXP,indEndR-indStartR+1));
		nProtect++;    
		for(j=indStartR;j<=indEndR;j++)
		{
			REAL(yR)[j-indStartR]=INTEGER(dataR)[j]+rnorm(0,2)*INTEGER(jitter)[0];
		}
		
		
		/*Process control data*/
		if((ncF>0) & (ncR>0))
		{
			while((pcF<ncF) & (INTEGER(contF)[pcF]<minLocF))
			{
				pcF++;
			}
			indStart=pcF;
			
			while((pcF<ncF) & (INTEGER(contF)[pcF]<=maxLocF))
			{
				pcF++;
			}
			indEnd=imin2(pcF,ncF);
			
			PROTECT(cF = allocVector(REALSXP,indEnd-indStart));
			nProtect++;      
			for(j=indStart;j<indEnd;j++)
			{
				REAL(cF)[j-indStart]=INTEGER(contF)[j]+rnorm(0,.1)*INTEGER(jitter)[0];
			}
			
			while((pcR<ncR) & (INTEGER(contR)[pcR]<minLocR))
			{
				pcR++;
			}
			indStart=pcR;
			
			while((pcR<ncR) & (INTEGER(contR)[pcR]<=maxLocR))
			{
				pcR++;
			}
			indEnd=imin2(pcR,ncR);
			
			PROTECT(cR = allocVector(REALSXP,indEnd-indStart));
			nProtect++;      
			for(j=indStart;j<indEnd;j++)
			{
				REAL(cR)[j-indStart]=INTEGER(contR)[j]+rnorm(0,.1)*INTEGER(jitter)[0];
			}
		}
		else
		{
			cF = R_NilValue;
			cR = R_NilValue;
		}
		
		classDef=MAKE_CLASS("segReads");
		PROTECT(seg=NEW_OBJECT(classDef));
		nProtect++;    
		SET_SLOT(seg,mkChar("yF"),yF);
		SET_SLOT(seg,mkChar("yR"),yR);
		SET_SLOT(seg,mkChar("cF"),cF);
		SET_SLOT(seg,mkChar("cR"),cR);
		SET_SLOT(seg,mkChar("map"),map);
		SET_SLOT(seg,mkChar("chr"),name);
		SET_VECTOR_ELT(ans,i,seg);
		UNPROTECT(nProtect);
	}
	UNPROTECT(2);
	PutRNGstate();
	return(ans);
}

SEXP getVector(SEXP list, SEXP ind)
{
  int i=0, j=0, K=0, n=length(list), nTotal=0, counter=0;
  SEXP ans;
  
  /** First we get the length of the object **/
  for(i=0;i<n;i++)
  {
    /** Make sure it's not an object of class picsError **/
    //if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
    if(testClass(list,i))
    {
      nTotal+=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
    }
  }

  /** Allocate the memory for the results **/
  PROTECT(ans = allocVector(REALSXP, nTotal));

  for(i=0;i<n;i++)
  {
    //if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
    if(testClass(list,i))
    {
      K=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
      for(j=0;j<K;j++)
      {      
        REAL(ans)[counter]=REAL(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), INTEGER(ind)[0]))[j];
        counter++;
      }
    }
  }
  UNPROTECT(1);
  return(ans);
}



SEXP getScore(SEXP list)
{
  int i=0, j=0, K=0, n=length(list), nTotal=0, counter=0;
  SEXP ans;
  
  /** First we get the length of the object **/
  for(i=0;i<n;i++)
  {
    /** Make sure it's not an object of class picsError **/
    //if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
    if(testClass(list,i))
    {
      nTotal+=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
    }
  }

  /** Allocate the memory for the results **/
  PROTECT(ans = allocVector(REALSXP, nTotal));

  for(i=0;i<n;i++)
  {
    //if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
    if(testClass(list,i))
    {
      K=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
      for(j=0;j<K;j++)
      {      
        REAL(ans)[counter]=REAL(GET_SLOT(VECTOR_ELT(list, i),install("score")))[j];
        counter++;
      }
    }
  }
  UNPROTECT(1);
  return(ans);
}

SEXP getScoreF(SEXP list)
{
	int i=0, j=0, K=0, n=length(list), nTotal=0, counter=0;
	SEXP ans;
	
	/** First we get the length of the object **/
	for(i=0;i<n;i++)
	{
		/** Make sure it's not an object of class picsError **/
		//if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
    		if(testClass(list,i))
		{
			nTotal+=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
		}
	}
	
	/** Allocate the memory for the results **/
	PROTECT(ans = allocVector(REALSXP, nTotal));
	
	for(i=0;i<n;i++)
	{
		//if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
    		if(testClass(list,i))
		{
			K=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
			for(j=0;j<K;j++)
			{      
				REAL(ans)[counter]=REAL(GET_SLOT(VECTOR_ELT(list, i),install("scoreF")))[j];
				counter++;
			}
		}
	}
	UNPROTECT(1);
	return(ans);
}

SEXP getScoreR(SEXP list)
{
	int i=0, j=0, K=0, n=length(list), nTotal=0, counter=0;
	SEXP ans;
	
	/** First we get the length of the object **/
	for(i=0;i<n;i++)
	{
		/** Make sure it's not an object of class picsError **/
		//if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
    		if(testClass(list,i))
		{
			nTotal+=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
		}
	}
	
	/** Allocate the memory for the results **/
	PROTECT(ans = allocVector(REALSXP, nTotal));
	
	for(i=0;i<n;i++)
	{
		//if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
    		if(testClass(list,i))
		{
			K=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
			for(j=0;j<K;j++)
			{      
				REAL(ans)[counter]=REAL(GET_SLOT(VECTOR_ELT(list, i),install("scoreR")))[j];
				counter++;
			}
		}
	}
	UNPROTECT(1);
	return(ans);
}

SEXP getChr(SEXP list)
{
  int i=0, j=0, K=0, n=length(list), nTotal=0, counter=0;
  SEXP ans;
  
  /** First we get the length of the object **/
  for(i=0;i<n;i++)
  {
    /** Make sure it's not an object of class picsError **/
    //if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
    if(testClass(list,i))
    {
      nTotal+=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
    }
  }

  /** Allocate the memory for the results **/
  PROTECT(ans = allocVector(STRSXP, nTotal));

  for(i=0;i<n;i++)
  {
    //if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
    if(testClass(list,i))
    {
      K=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
      for(j=0;j<K;j++)
      {      
        SET_STRING_ELT(ans,counter, STRING_ELT(GET_SLOT(VECTOR_ELT(list, i),install("chr")),0));
        counter++;
      }
    }
  }
  UNPROTECT(1);
  return(ans);
}

SEXP getMin(SEXP list)
{
	int i=0, j=0, K=0, n=length(list), nTotal=0, counter=0;
	SEXP ans;
	
	/** First we get the length of the object **/
	for(i=0;i<n;i++)
	{
		/** Make sure it's not an object of class picsError **/
		//if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
    		if(testClass(list,i))
		{
			nTotal+=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
		}
	}
	
	/** Allocate the memory for the results **/
	PROTECT(ans = allocVector(REALSXP, nTotal));
	
	for(i=0;i<n;i++)
	{
		//if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
    		if(testClass(list,i))
		{
			K=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
			for(j=0;j<K;j++)
			{      
				REAL(ans)[counter]= REAL(GET_SLOT(VECTOR_ELT(list, i),install("range")))[0];
				counter++;
			}
		}
	}
	UNPROTECT(1);
	return(ans);
}

SEXP getMax(SEXP list)
{
	int i=0, j=0, K=0, n=length(list), nTotal=0, counter=0;
	SEXP ans;
	
	/** First we get the length of the object **/
	for(i=0;i<n;i++)
	{
		/** Make sure it's not an object of class picsError **/
		//if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
    		if(testClass(list,i))
		{
			nTotal+=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
		}
	}
	
	/** Allocate the memory for the results **/
	PROTECT(ans = allocVector(REALSXP, nTotal));
	
	for(i=0;i<n;i++)
	{
		//if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
    		if(testClass(list,i))
		{
			K=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
			for(j=0;j<K;j++)
			{      
				REAL(ans)[counter]= REAL(GET_SLOT(VECTOR_ELT(list, i),install("range")))[1];
				counter++;
			}
		}
	}
	UNPROTECT(1);
	return(ans);
}


SEXP getMap(SEXP list)
{
  int i=0, j=0, n=length(list), nMap=0, lyF, lyR;
  double sumDiff=0, m=0, M=1, *yF, *yR;
  int *map;
  SEXP ans;

  PROTECT(ans = allocVector(REALSXP, n));

  /** First we get the length of the object **/
  for(i=0;i<n;i++)
  {
    /** Make sure it's not an object of class picsError **/
    nMap=INTEGER(GET_DIM(GET_SLOT(VECTOR_ELT(list, i),install("map"))))[0];
    /** Check that we have an interesection with a mappability interval **/
    sumDiff=0;

    if(nMap>0)
    {
      map=INTEGER(GET_SLOT(VECTOR_ELT(list, i),install("map")));
      for(j=0;j<nMap;j++)
      {
        sumDiff+=map[nMap+j]-map[j];
      }

      yF=REAL(GET_SLOT(VECTOR_ELT(list, i),install("yF")));
      yR=REAL(GET_SLOT(VECTOR_ELT(list, i),install("yR")));
      lyF=length(GET_SLOT(VECTOR_ELT(list, i),install("yF")));
      lyR=length(GET_SLOT(VECTOR_ELT(list, i),install("yR")));

      m=fmin2(fmin2(yF[0],yR[0]),map[0]);
      M=fmax2(fmax2(yF[lyF-1],yR[lyR-1]),map[2*nMap-1]);
    }
    REAL(ans)[i]=sumDiff/fmax2(M-m,1.0);
  }  
  UNPROTECT(1);
  return(ans);
}

SEXP getK(SEXP list)
{
  int i=0, n=length(list);
  SEXP ans;

  PROTECT(ans = allocVector(REALSXP, n));
  
  /** First we get the length of the object **/
  for(i=0;i<n;i++)
  {
    /** Make sure it's not an object of class picsError **/
    //if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
    if(testClass(list,i))
    {
      REAL(ans)[i]=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
    }
    else
    {
      REAL(ans)[i]=0;
    }
  }
  UNPROTECT(1);
  return(ans);
}

SEXP getSegL(SEXP list)
{
	int i=0, j=0, n=length(list), nF, nR, nProtect=0;
	SEXP ans, tempF, tempR, tempM, NNF, NNR,  LL, mm, MM, chr;
	double min, max;

	PROTECT(LL = allocVector(REALSXP, n)); nProtect++;
	PROTECT(mm = allocVector(REALSXP, n)); nProtect++;
	PROTECT(MM = allocVector(REALSXP, n)); nProtect++;
	PROTECT(NNF= allocVector(INTSXP, n)); nProtect++;
	PROTECT(NNR= allocVector(INTSXP, n)); nProtect++;
	PROTECT(chr= allocVector(STRSXP, n)); nProtect++;
	PROTECT(ans= NEW_LIST(6));			  nProtect++;

	/** First we get the length of the object **/
	for(i=0;i<n;i++)
	{
		tempF=GET_SLOT(VECTOR_ELT(list, i),install("yF"));
		tempR=GET_SLOT(VECTOR_ELT(list, i),install("yR"));
		nF= length(tempF);
		nR= length(tempR);
		tempM=GET_SLOT(VECTOR_ELT(list, i),install("map"));
		
		/*
		Rprintf("Map: \n");
		for (j=0; j<length(tempM); j++) {
			Rprintf("%i, \t",INTEGER(tempM)[j]);
		}
		Rprintf("\n");
		Rprintf("nMap=%i \n",length(tempM));
		*/ 
		 
		if (length(tempM)>0) {
			min=imin2(REAL(tempF)[0],INTEGER(tempM)[0]);
			max=imax2(REAL(tempR)[nR-1],INTEGER(tempM)[length(tempM)-1]);
		}else 
		{
			min=REAL(tempF)[0];
			max=REAL(tempR)[nR-1];
		}


		REAL(MM)[i]=max;
		REAL(mm)[i]=min;
		REAL(LL)[i]=max-min;
		INTEGER(NNF)[i]=nF;
		INTEGER(NNR)[i]=nR;
		SET_STRING_ELT(chr,i, STRING_ELT(GET_SLOT(VECTOR_ELT(list, i),install("chr")),0));
	}
	
	SET_VECTOR_ELT(ans,0,chr);
	SET_VECTOR_ELT(ans,1,NNF);
	SET_VECTOR_ELT(ans,2,NNR);
	SET_VECTOR_ELT(ans,3,LL);
	SET_VECTOR_ELT(ans,4,mm);
	SET_VECTOR_ELT(ans,5,MM);

	
	UNPROTECT(nProtect);
	return(ans);
}

SEXP getDensityList(SEXP picsList, SEXP strand, SEXP step, SEXP filter, SEXP sum, SEXP scale)
{
  int l=0,i=0;
  SEXP pics, List, ans, ansTmp;
  SEXP x, y, chr, chromosome,names;
  int nProtected=0,totalLength=0;
  double *range;

  PROTECT(List=GET_SLOT(picsList,install("List"))); nProtected++;
  PROTECT(ans=NEW_LIST(3)); nProtected++;

  for(l=0;l<length(List);l++)
  {
    pics=VECTOR_ELT(List,l);
    if((strcmp(CHAR(STRING_ELT(GET_CLASS(pics), 0)),"pics")==0) || (strcmp(CHAR(STRING_ELT(GET_CLASS(pics), 0)),"ping")==0))
    {
      range=REAL(GET_SLOT(pics,install("range")));
      totalLength+=(int)((range[1]-range[0])/REAL(step)[0]);
    }
  }

  PROTECT(x=NEW_NUMERIC(totalLength));nProtected++;
  PROTECT(y=NEW_NUMERIC(totalLength));nProtected++;
  PROTECT(chr=NEW_STRING(totalLength));nProtected++;
  
  totalLength=0;
  for(l=0;l<length(List);l++)
  { 
    pics=VECTOR_ELT(List,l);
    //if(strcmp(CHAR(STRING_ELT(GET_CLASS(pics), 0)),"pics")==0)
    if((strcmp(CHAR(STRING_ELT(GET_CLASS(pics), 0)),"pics")==0) || (strcmp(CHAR(STRING_ELT(GET_CLASS(pics), 0)),"ping")==0))
    {
    
      chromosome=GET_SLOT(pics,install("chr"));
      range=REAL(GET_SLOT(pics,install("range")));

      PROTECT(ansTmp=getDensity(pics, strand, step, filter, sum, scale));nProtected++;
      for(i=0;i<(int)((range[1]-range[0])/REAL(step)[0]);i++)
      {
        REAL(x)[i+totalLength]=REAL(VECTOR_ELT(ansTmp,0))[i];
        REAL(y)[i+totalLength]=REAL(VECTOR_ELT(ansTmp,1))[i];
        STRING_PTR(chr)[i+totalLength]=STRING_PTR(chromosome)[0];
      }    
      totalLength+=(int)((range[1]-range[0])/REAL(step)[0]);
      UNPROTECT(1);nProtected--;
    }
  }
    // I have added the option to interrupt R
    // R_CheckUserInterrupt();
  SET_VECTOR_ELT(ans,0,x);
  SET_VECTOR_ELT(ans,1,y);
  SET_VECTOR_ELT(ans,2,chr);
  PROTECT(names = allocVector(STRSXP, 3)); nProtected++; 
  SET_STRING_ELT(names, 0, mkChar("x"));
  SET_STRING_ELT(names, 1, mkChar("density"));
  SET_STRING_ELT(names, 2, mkChar("chr"));
  setAttrib(ans, R_NamesSymbol, names);  
  UNPROTECT(nProtected);
  return(ans);
}

SEXP getDensity(SEXP pics, SEXP strand, SEXP step, SEXP filter, SEXP sum, SEXP scale)
{
  int k=0,K=1,i=0;
  double *w, *mu, *delta, *sF, *sR, *range, *se, sumW=0;
  double *muFilter, *deltaFilter, *sFFilter, *sRFilter, *seFilter, *seFilterF, *seFilterR, *scoreFilter, *score;
  SEXP filter1;
  gsl_vector *One;  
  gsl_matrix *Density;
  int length, nProtected=0;
  SEXP ans, x, y, chr, names;
  gsl_vector_view view;
  chr=GET_SLOT(pics,install("chr"));
  
  if(strcmp(CHAR(STRING_ELT(GET_CLASS(pics), 0)),"picsError")==0)
  {
    return(R_NilValue);
  }
  else
  {
    K=length(VECTOR_ELT(GET_SLOT(pics,install("estimates")), 0));
    w=REAL(VECTOR_ELT(GET_SLOT(pics,install("estimates")), 0));
    mu=REAL(VECTOR_ELT(GET_SLOT(pics,install("estimates")), 1));
    delta=REAL(VECTOR_ELT(GET_SLOT(pics,install("estimates")),2));
    sF=REAL(VECTOR_ELT(GET_SLOT(pics,install("estimates")), 3));
    sR=REAL(VECTOR_ELT(GET_SLOT(pics,install("estimates")), 4));
    score=REAL(GET_SLOT(pics,install("score")));
    se=REAL(VECTOR_ELT(GET_SLOT(pics,install("estimates")), 5));

    deltaFilter=REAL(getListElement(filter, "delta"));
    sFFilter=REAL(getListElement(filter, "sigmaSqF"));
    sRFilter=REAL(getListElement(filter, "sigmaSqR"));
    seFilter=REAL(getListElement(filter, "se"));
    seFilterF=REAL(getListElement(filter, "seF"));
    seFilterR=REAL(getListElement(filter, "seR"));
    scoreFilter=REAL(getListElement(filter, "score"));
    //Rprintf("seFF=%f, seFR=%f\n", seFilterF, seFilterR);

    range=REAL(GET_SLOT(pics,install("range")));

    PROTECT(ans=NEW_LIST(2));nProtected++;

    // Compute the number of steps I want
    length=(int)((range[1]-range[0])/REAL(step)[0]);
    // Allocate the memory for x and y
    PROTECT(x=NEW_NUMERIC(length));nProtected++;
    PROTECT(y=NEW_NUMERIC(length));nProtected++;
    Density=gsl_matrix_calloc(length,K);
    One=gsl_vector_alloc(K);
    gsl_vector_set_all(One, 1.0);

    for(i=0;i<length;i++)
    {
      //Grid
      REAL(x)[i]=range[0]+i*REAL(step)[0];

      for(k=0;k<K;k++)
      {
      //Check if this is a valid binding event
        if((delta[k]>deltaFilter[0] & delta[k]<deltaFilter[1]) & 
          (sF[k]>sFFilter[0] & sF[k]<sFFilter[1]) & 
          (sR[k]>sRFilter[0] & sR[k]<sRFilter[1]) & 
          (se[k]>seFilter[0] & se[k]<seFilter[1]) &
          (se[k]>seFilterF[0] & se[k]<seFilterF[1]) &
          (se[k]>seFilterR[0] & se[k]<seFilterR[1]) &          
          (score[k]>scoreFilter[0] & score[k]<scoreFilter[1]))
        {
        //Keep track of the sum of the weights to renormalize the density
          sumW+=w[k];
          if(REAL(strand)[0]==1)
          {
            gsl_matrix_set(Density,i,k,w[k]*gsl_ran_tdist_pdf((REAL(x)[i]-mu[k]+delta[k]/2.0)/sqrt(sF[k]), 4.0)/sqrt(sF[k]));
          }
          else if(REAL(strand)[0]==-1)
          {
            gsl_matrix_set(Density,i,k,w[k]*gsl_ran_tdist_pdf((REAL(x)[i]-mu[k]-delta[k]/2.0)/sqrt(sR[k]), 4.0)/sqrt(sR[k]));
          }
          else if(REAL(strand)[0]==0)
          {
            //Here I compute the shift density
            gsl_matrix_set(Density,i,k,0.5*w[k]*gsl_ran_tdist_pdf((REAL(x)[i]-mu[k])/sqrt(sR[k]), 4.0)/sqrt(sR[k]));
            gsl_matrix_set(Density,i,k,gsl_matrix_get(Density,i,k)+0.5*w[k]*gsl_ran_tdist_pdf((REAL(x)[i]-mu[k])/sqrt(sF[k]), 4.0)/sqrt(sF[k]));
          }
          // Here I scale by the score
          if(LOGICAL(scale)[0])
          {
            gsl_matrix_set(Density,i,k,gsl_matrix_get(Density,i,k)*score[k]);
          }
        }
      }
    }
      // We compute the overall density
    if(LOGICAL(sum)[0])
    {
          //Sum up over all binding events
      view=gsl_vector_view_array(REAL(y),length);
      gsl_blas_dgemv(CblasNoTrans, 1.0, Density, One, 0, &view.vector);
          //Rescale
      if(sumW>0)
      {
        gsl_vector_scale(&view.vector, 1./sumW);
      }          
    }//Otherwise we take the most likely event at each position
    else
    {
      for(i=0;i<length;i++)
      {
        view=gsl_matrix_row(Density,i);
        REAL(y)[i]=gsl_vector_max(&view.vector);
      }
    }

    SET_VECTOR_ELT(ans,0,x);
    SET_VECTOR_ELT(ans,1,y);
    PROTECT(names = allocVector(STRSXP, 2)); nProtected++; 
    SET_STRING_ELT(names, 0, mkChar("x"));
    SET_STRING_ELT(names, 1, mkChar("density"));
    setAttrib(ans, R_NamesSymbol, names);  

    gsl_vector_free(One);
    gsl_matrix_free(Density);
    UNPROTECT(nProtected++);
    return(ans);
  }
}

SEXP getListElement(SEXP list, const char *str)
     {
       SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
       int i;
     
       for (i = 0; i < length(list); i++)
         if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
           elmt = VECTOR_ELT(list, i);
           break;
         }
       return elmt;
     }

//This function test if every object of the list is a "pics" or "ping"
int testClass(SEXP list, int i)
{
  int res=1;
  if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0 || strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list,     i)), 0)),"ping")==0)
  { res=1; }
  else
  { res=0; }
  return(res);
}

//This function finds out if the object has been called from 'PICS' or 'PING'
int testObj(SEXP pPackage)
{
  int res=1;
  if(strcmp(CHAR(STRING_ELT(pPackage,0)),"PICS")==0){
    res=1;
  } else{
    res=0;
  }
  return(res);
}
