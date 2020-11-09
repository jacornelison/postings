
#include <math.h>

#include "mex.h"
#include "matrix.h"

#include "point.h"

void mexFunction(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[])
{
  double *pjd;
  long i,npts;
  double *pra,*pdec,*pelon;

  if(nrhs!=1)
    mexErrMsgTxt("Must have 1 input arguments");
  if(nlhs!=3)
    mexErrMsgTxt("Must have 3 output arguments");

  /* don't bother with further error checking here in C -
     do that in Matlab wrapper function */

  /* get pointers to input arrays */
  pjd=mxGetPr(prhs[0]);
  npts=mxGetNumberOfElements(prhs[0]);

  /* create output arrays with same size */
  plhs[0]=mxCreateDoubleMatrix(npts,1,mxREAL);
  plhs[1]=mxCreateDoubleMatrix(npts,1,mxREAL);
  plhs[2]=mxCreateDoubleMatrix(npts,1,mxREAL);

  /* get pointers to real part of these */
  pra=mxGetPr(plhs[0]);
  pdec=mxGetPr(plhs[1]);
  pelon=mxGetPr(plhs[2]);

  /* call the function */
  for (i=0; i<npts; i++)
    sunpos(pjd[i],pra+i,pdec+i,pelon+i);
}
