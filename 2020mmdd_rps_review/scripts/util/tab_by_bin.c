#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#ifndef HAVE_OCTAVE
#  include "matrix.h"
#endif

void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    double * tab_out;
    double * ibin;
    double * xdat;
    mxArray * mx_tab_out;
    int nchan, nsamp, nbins;
    double xx;
    int nn;
    int num_nans = 0;

    int i, j;

    if (nrhs < 3)
        mexErrMsgTxt ("Usage: tab = tab_by_bins (x, i, n)");

    ibin = mxGetData (prhs[1]);
    xdat = mxGetData (prhs[0]);
    nbins = mxGetScalar (prhs[2]);
    nchan = mxGetM (prhs[0]);
    nsamp = mxGetN (prhs[0]);

    printf ("%d channels, %d samples, %d bins.\n", nchan, nsamp, nbins);
    mx_tab_out = mxCreateDoubleMatrix (nchan*3, nbins, mxREAL);
    tab_out = mxGetData (mx_tab_out);
    for (i=0; i<(nchan*3)*nbins; i++)
      tab_out[i] = 0;

    for (i=0; i<nchan; i++)
    {
      for (j=0; j<nsamp; j++)
      {
        xx = xdat[i+j*nchan];
	nn = ibin[j];
        if ((nn < 0) || (nn > nbins))
        {
          mxDestroyArray (mx_tab_out);
          printf ("Bin %d out of bounds! (max is %d.)\n", nn, nbins);
          mexErrMsgTxt ("Exiting tab_by_bin.\n");
        }
        if ((xx == xx) && (nn != 0)) 
        {
          tab_out[(3*i+0) + (nn-1) * (3*nchan)] += 1;
          tab_out[(3*i+1) + (nn-1) * (3*nchan)] += xx;
          tab_out[(3*i+2) + (nn-1) * (3*nchan)] += xx * xx;
        }
	if (xx != xx)
          num_nans++;
      }
    }

    printf ("Found %d NaNs.\n", num_nans);
    plhs[0] = mx_tab_out;

    return;
}

