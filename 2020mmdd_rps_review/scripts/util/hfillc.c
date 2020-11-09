/*
   Filling histograms in matlab code is too slow.
   Code the main loop in C but keep a .m file wrapper
   to allow easy filling in of omitted arguments etc. 
   CLP 12/1/00
*/

#include <math.h>

#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[])
{
  double *px,*py,*pw;
  long   nx,ny,nw;
  int    nxbin,nybin;
  double lx,hx,ly,hy;
  double *px_tic,*py_tic,*pn;
  double wx,wy;
  long   i;
  int    xbin,ybin;

  if(nrhs<4||nrhs>5)
    mexErrMsgTxt("Must have 4 or 5 input arguments");
  if(nlhs!=2)
    mexErrMsgTxt("Must have 2 output arguments");

  px   =mxGetPr(prhs[0]);
  if(nrhs==5)
    pw   =mxGetPr(prhs[4]);

  nx=mxGetNumberOfElements(prhs[0]);
  if(nrhs==5)
    nw=mxGetNumberOfElements(prhs[4]);

  if(nrhs==5 && (nx!=nw))
    mexErrMsgTxt("number of elements in x and w must be equal");

  nxbin=mxGetScalar(prhs[1]);
  lx   =mxGetScalar(prhs[2]);
  hx   =mxGetScalar(prhs[3]);

  if(nxbin<1)
    mexErrMsgTxt("number of bins must be greater than zero");
  if(lx>=hx)
    mexErrMsgTxt("must have low<high");

  plhs[0]=mxCreateDoubleMatrix(1,nxbin,mxREAL);
  plhs[1]=mxCreateDoubleMatrix(1,nxbin,mxREAL);

  px_tic=mxGetPr(plhs[0]);
  pn    =mxGetPr(plhs[1]);

  wx=(hx-lx)/nxbin;
  for(i=0;i<nxbin;i++)
    px_tic[i]=lx+(i+0.5)*wx;

  for(i=0;i<nx;i++)
  {
    xbin=floor((px[i]-lx)/wx);
    if(xbin>=0 && xbin<nxbin)
      if(nrhs==4)
	pn[xbin]++;
      else
	pn[xbin]+=pw[i];
  }
}
