
#include <math.h>

#include "mex.h"
#include "matrix.h"

#include "point.h"

void mexFunction(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[])
{
  double *paz,*pel,*pjd,*plat,*plon;
  long i,npts;
  double *pra,*pdec;

  if(nrhs!=5)
    mexErrMsgTxt("Must have 5 input arguments");
  if(nlhs!=2)
    mexErrMsgTxt("Must have 2 output arguments");

  /* don't bother with further error checking here in C -
     do that in Matlab wrapper function */

  /* get pointers to input arrays */
  paz=mxGetPr(prhs[0]);
  pel=mxGetPr(prhs[1]);
  pjd=mxGetPr(prhs[2]);
  plat=mxGetPr(prhs[3]);
  plon=mxGetPr(prhs[4]);
  npts=mxGetNumberOfElements(prhs[0]);

  /* create output arrays with same size */
  plhs[0]=mxCreateDoubleMatrix(npts,1,mxREAL);
  plhs[1]=mxCreateDoubleMatrix(npts,1,mxREAL);

  /* get pointers to real part of these */
  pra=mxGetPr(plhs[0]);
  pdec=mxGetPr(plhs[1]);

  /* copy az/el into output arrays
     - azel2radec does conversion in place */
  for(i=0;i<npts;i++) {
    pra[i]=paz[i];
    pdec[i]=pel[i];
  }

  /* call the function - convert az and el in place to ra dec */
  azel2radec(pra,pdec,pjd,plat,plon,npts);
}
