#include <math.h>

#include "cfortran.h"
#include "minuit.h"
#include "mex.h"

int npar;
char gof_func[200];
double *p;
int    nfrhs;
mxArray *pflhs[1],*pfrhs[100];

void gof(int nvpar, double * grad, double * fval,
	 double par[], int iflag, void ( * A6 )() )
{
  int i;
  for(i=0;i<npar;i++)
     p[i]=par[i];
  mexCallMATLAB(1,pflhs,nfrhs,pfrhs,gof_func);

  *fval=*mxGetPr(pflhs[0]);
}
FCALLSCSUB6(gof,GOF,gof,
	    INT,PDOUBLE,PDOUBLE,DOUBLEV,INT,ROUTINE);

void mexFunction(int nlhs,
		 mxArray *plhs[],
		 int nrhs,
		 const mxArray *prhs[])
{
  double *inpar, *freepar, *lb, *ub;
  int i,j;
  int err; /* error return flag for Minuit calls */
  char cmd[100];

  if(nrhs<8)
    mexErrMsgTxt("Must have at least 8 input arguments");

  mxGetString(prhs[0],gof_func,200);

  npar = mxGetNumberOfElements(prhs[1]);
  inpar = mxGetPr(prhs[1]);
  freepar = mxGetPr(prhs[2]);
  lb=mxGetPr(prhs[3]);
  ub=mxGetPr(prhs[4]);

  /* Prepare gof function args */
  /* Array used to pass par to MATLAB */
  pfrhs[0]=mxDuplicateArray(prhs[1]); /* Duplicate to preserve shape */
  p=mxGetPr(pfrhs[0]);
  for(i=5;i<nrhs;i++)                           /* Copy optional args */
    pfrhs[i-4]=prhs[i];
  nfrhs=nrhs-4;

  MNINIT(5,6,7);
  sprintf(cmd,"set pri -1"); MNCOMD(NULL,cmd,err,0);

  for(i=0;i<npar;i++)  
  {
    char str[10];
    sprintf(str,"p%02d",i+1);

    /* Enter starting par values */
    if(inpar[i]==0)
      MNPARM(i+1,str,0,1e-6,lb[i],ub[i],err);
    else
      MNPARM(i+1,str,inpar[i],fabs(inpar[i]/10),lb[i],ub[i],err);

    /* fix par if requested */
    if(freepar[i]==0)
    {
      sprintf(str,"fix %d",i+1);
      MNCOMD(NULL,str,err,0);
    }
  }

  /* Do minimization */
  sprintf(cmd,"mini"); MNCOMD(gof_,cmd,err,0);

  /* Calc assym errors and show summary data
     (only parabolic errors are currently returned) */
  sprintf(cmd,"set pri 0"); MNCOMD(NULL,cmd,err,0);
  sprintf(cmd,"mino"); MNCOMD(gof_,cmd,err,0);

  /* Show the parameter correlation matrix */
  sprintf(cmd,"set pri -1"); MNCOMD(NULL,cmd,err,0);
  sprintf(cmd,"show corr"); MNCOMD(gof_,cmd,err,0);

  if(nlhs>0) /* Get par values */
  {
    double *outpar;
    plhs[0]=mxDuplicateArray(prhs[1]); /* Duplicate to preserve shape */
    outpar  = mxGetPr(plhs[0]);
    for(i=0;i<npar;i++)
    {
      char dstr[100];
      int  dint;
      double ddbl;
      MNPOUT(i+1,dstr,outpar[i],ddbl,ddbl,ddbl,dint);
    }
  }

  if(nlhs>1) /* Get par errors */
  {
    double *parerr;
    plhs[1]=mxDuplicateArray(prhs[1]); /* Duplicate to preserve shape */
    parerr  = mxGetPr(plhs[1]);
    for(i=0;i<npar;i++)
    {
      char dstr[100];
      int  dint;
      double ddbl;
      MNPOUT(i+1,dstr,ddbl,parerr[i],ddbl,ddbl,dint);
    }
  }

  {  /* Get gof and status */
    double gof,ddbl;
    int dint,stat;
    MNSTAT(gof,ddbl,ddbl,dint,dint,stat);

    if(nlhs>2)
    {
      plhs[2]=mxCreateDoubleMatrix(1,1,mxREAL);
      *mxGetPr(plhs[2])=gof;
    }
    if(nlhs>3)
    {
      plhs[3]=mxCreateDoubleMatrix(1,1,mxREAL);
      *mxGetPr(plhs[3])=stat;
    }
  }

  if(nlhs>4) /* Get covariance matrix */
  {  
    double *covar;
    int *ndims=mxGetDimensions(prhs[1]);
    plhs[4]=mxCreateDoubleMatrix(ndims[1],ndims[1],mxREAL);
    
    covar=mxGetPr(plhs[4]);
    MNEMAT(*covar,ndims[1]);
    
  }
}

int intrac_()
{
  return C2FLOGICAL(FALSE);
}
