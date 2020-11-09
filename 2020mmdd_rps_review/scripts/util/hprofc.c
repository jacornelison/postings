/*
   Doing profile histogram in matlab code is too slow
   Code the main loop in C but keep a .m file wrapper
   CLP 11/20/00
*/

#include <math.h>

#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[])
{
  double *px,*py;
  long   nx,ny;
  int    nbin;
  double low,high;
  double *ptic,*pmu,*psig,*pn;
  double bw;
  long   i;
  int    bin;

  if(nrhs!=5)
    mexErrMsgTxt("Must have 5 input arguments");
  if(nlhs!=4)
    mexErrMsgTxt("Must have 4 output arguments");

  px   =mxGetPr(prhs[0]);
  py   =mxGetPr(prhs[1]);

  nx=mxGetNumberOfElements(prhs[0]);
  ny=mxGetNumberOfElements(prhs[1]);

  if(nx!=ny)
    mexErrMsgTxt("number of elements in x and y must be equal");

  nbin=mxGetScalar(prhs[2]);
  low =mxGetScalar(prhs[3]);
  high=mxGetScalar(prhs[4]);

  if(nbin<1)
    mexErrMsgTxt("number of bins must be greater than zero");
  if(low>=high)
    mexErrMsgTxt("must have low<high");

  plhs[0]=mxCreateDoubleMatrix(1,nbin,mxREAL);
  plhs[1]=mxCreateDoubleMatrix(1,nbin,mxREAL);
  plhs[2]=mxCreateDoubleMatrix(1,nbin,mxREAL);
  plhs[3]=mxCreateDoubleMatrix(1,nbin,mxREAL);

  ptic=mxGetPr(plhs[0]);
  pmu =mxGetPr(plhs[1]);
  psig=mxGetPr(plhs[2]);
  pn  =mxGetPr(plhs[3]);

  bw=(high-low)/nbin;
  for(i=0;i<nbin;i++)
    ptic[i]=low+(i+0.5)*bw;

  /* Calc mean */
  for(i=0;i<nx;i++)
  {
    bin=floor((px[i]-low)/bw);
    if(bin>=0 && bin<nbin)
    {
      pn[bin]++;
      pmu[bin]+=py[i];
    }
  }
  for(bin=0;bin<nbin;bin++)
    if(pn[bin]>1)
      pmu[bin]/=pn[bin];

  /* Calc sigma */
  for(i=0;i<nx;i++)
  {
    bin=floor((px[i]-low)/bw);
    if(bin>=0 && bin<nbin)
      psig[bin]+=pow(py[i]-pmu[bin],2);
  }
  for(bin=0;bin<nbin;bin++)
    if(pn[bin]>1)
      psig[bin]=sqrt((1/(pn[bin]-1))*psig[bin]);
    else
      psig[bin]=0;
}
