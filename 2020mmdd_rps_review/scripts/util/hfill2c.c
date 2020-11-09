/*
   Make and fill a 2 dimensional histogram
   CLP 11/3/00
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
  int    nxbin,nybin,opt;
  double lx,hx,ly,hy;
  double *px_tic,*py_tic,*pn;
  double wx,wy;
  long   i;
  int    xbin,ybin;
  
  /* designed to be called from C-wrapper so we require all args present */
  if(nrhs<10||nrhs>10)
    mexErrMsgTxt("Must have 10 input arguments");
  if(nlhs<3||nlhs>3)
    mexErrMsgTxt("Must have 3 output arguments");

  px   =mxGetPr(prhs[0]);
  py   =mxGetPr(prhs[1]);
  pw   =mxGetPr(prhs[8]);

  nx=mxGetNumberOfElements(prhs[0]);
  ny=mxGetNumberOfElements(prhs[1]);
  nw=mxGetNumberOfElements(prhs[8]);

  nxbin=mxGetScalar(prhs[2]);
  lx   =mxGetScalar(prhs[3]);
  hx   =mxGetScalar(prhs[4]);

  nybin=mxGetScalar(prhs[5]);
  ly   =mxGetScalar(prhs[6]);
  hy   =mxGetScalar(prhs[7]);

  opt=mxGetScalar(prhs[9]);

  if(nxbin<1 || nybin<1)
    mexErrMsgTxt("number of x and y bins must be greater than zero");
  if(lx>=hx || ly>=hy)
    mexErrMsgTxt("must have lx<hx and ly<hy");

  plhs[0]=mxCreateDoubleMatrix(1,nxbin,mxREAL);
  plhs[1]=mxCreateDoubleMatrix(1,nybin,mxREAL);
  plhs[2]=mxCreateDoubleMatrix(nybin,nxbin,mxREAL);

  px_tic=mxGetPr(plhs[0]);
  py_tic=mxGetPr(plhs[1]);
  pn    =mxGetPr(plhs[2]);

  if(opt>0)
    for(i=0;i<nxbin*nybin;i++)
      pn[i]=sqrt(-1);

  wx=(hx-lx)/nxbin;
  wy=(hy-ly)/nybin;
  for(i=0;i<nxbin;i++)
    px_tic[i]=lx+(i+0.5)*wx;
  for(i=0;i<nybin;i++)
    py_tic[i]=ly+(i+0.5)*wy;

  for(i=0;i<nx;i++)
  {
    xbin=floor((px[i]-lx)/wx);
    ybin=floor((py[i]-ly)/wy);
    if(xbin>=0 && xbin<nxbin && ybin>=0 && ybin<nybin)
      if(opt==0)
	/* sum the bin contents */
	pn[xbin*nybin+ybin]+=pw[i];
      else
	/* find max/min in bin */
	switch(opt)
	{
	  case 1:
	    if(pn[xbin*nybin+ybin]<pw[i] | isnan(pn[xbin*nybin+ybin]))
	      pn[xbin*nybin+ybin]=pw[i];
	    break;
	  case 2:
	    if(pn[xbin*nybin+ybin]>pw[i] | isnan(pn[xbin*nybin+ybin]))
	      pn[xbin*nybin+ybin]=pw[i];
	    break;
	}
  }
}
