
#include <math.h>

#include "mex.h"
#include "matrix.h"

#include "point.h"

void mexFunction(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[])
{
  double *pra,*pdec,*pjd,*plat,*plon;
  long i,npts;
  double *paz,*pel;

  if(nrhs!=5)
    mexErrMsgTxt("Must have 5 input arguments");
  if(nlhs!=2)
    mexErrMsgTxt("Must have 2 output arguments");

  /* don't bother with further error checking here in C -
     do that in Matlab wrapper function */

  /* get pointers to input arrays */
  pra=mxGetPr(prhs[0]);
  pdec=mxGetPr(prhs[1]);
  pjd=mxGetPr(prhs[2]);
  plat=mxGetPr(prhs[3]);
  plon=mxGetPr(prhs[4]);
  npts=mxGetNumberOfElements(prhs[0]);

  /* create output arrays with same size */
  plhs[0]=mxCreateDoubleMatrix(npts,1,mxREAL);
  plhs[1]=mxCreateDoubleMatrix(npts,1,mxREAL);

  /* get pointers to real part of these */
  paz=mxGetPr(plhs[0]);
  pel=mxGetPr(plhs[1]);

  /* copy ra/dec into output arrays
     - radec2azel does conversion in place */
  for(i=0;i<npts;i++) {
    paz[i]=pra[i];
    pel[i]=pdec[i];
  }

  /* call the function - convert ra and dec in place to az el */
  radec2azel(paz,pel,pjd,plat,plon,npts);
}
