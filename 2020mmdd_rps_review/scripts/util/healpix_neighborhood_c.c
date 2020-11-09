/*
 *  Helper function for healpix neighborhoods of finite radius.
 *  Given a table of nearest-neighbor pixels, desired radius r,
 *  and pixel center locations, find all pixels within neigh-
 *  borhoods of size r.
 *  RWO 12/08/11
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#ifndef HAVE_OCTAVE
#  include "matrix.h"
#endif
#include <math.h>

/*
 * Structure for keeping lists of pixels --
 * needs to be easy to expand if needed.
 */
struct pixlist
{
  int * val;
  int n;
  int nused;
};

/* initialize pixlist structure */
void init_pixlist (struct pixlist * pl, int n)
{
  pl->n = n;
  pl->nused = 0;
  pl->val = (int *)malloc (n * sizeof(int));
}

/* ... and here's the easy way to expand it. */
void expand_pixlist (struct pixlist * pl, int n)
{
  printf("Expanding pixlist from %d to %d.\n", pl->n, n);
  pl->val = (int *)realloc (pl->val, n * sizeof(int));
  pl->n = n;
}

/* free allocated memory in pixlist when done. */
void free_pixlist (struct pixlist * pl)
{
  free(pl->val);
  pl->n = 0;
  pl->nused = 0;
}

/*
 * Add a pixel to the list, if it's not already present.
 * If it is already present, do nothing.
 * Expand the list if it's already full.
 */
void add_to_pixlist (struct pixlist * pl, int v)
{
  /* printf("Adding %d to pixlist.\n", v); */
  int i;
  for (i=0; i<pl->nused; i++)
  {
    if (pl->val[i] == v)
      return;
  }
  if (pl->nused == pl->n)
    expand_pixlist (pl, pl->n + 1);
  pl->val[pl->nused] = v;
  pl->nused++;
}

/*
 * Simple bubble sort.  Not the world's
 * fastest sort algorithm, but close enough.
 */
void sort_pixlist (struct pixlist * pl)
{
  int i;
  int done = 0;
  int tmp;

  while (done == 0)
  {
    done = 1;
    for (i=1; i<pl->nused; i++)
    {
      if (pl->val[i] < pl->val[i-1])
      {
        tmp = pl->val[i];
        pl->val[i] = pl->val[i-1];
        pl->val[i-1] = tmp;
        done = 0;
      }
    }
  }
}

void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    int npix;            /* # pixels in map */
    int maxn0;           /* Max # nearest neighbors (should be 8) */
    int maxnout;         /* Max # of pixels in neighborhood (depends on r) */

    double * n0;         /* Input array: nearest neighbor table */
    double * pxvec;      /* Input array: npix by 3 list of pixel center unit vectors */
    double * neighb_out; /* Output array: neighborhood table */
    mxArray * mx_neighb_out; 

    struct pixlist pl;   /* Temporary neighborhood list */
    int i, j, k, l, m;   /* Loop variables */

    double r;            /* Requested neighborhood size */
    double dpcut;        /* Corresponding dot-product threshold */

    if (nrhs < 3)
        mexErrMsgTxt ("Usage: neighb = hpix_neighborhood (n0, pxvec, r)");

    /* Get inputs */
    npix = mxGetM (prhs[0]);
    maxn0 = mxGetN (prhs[0]);
    n0 = mxGetPr (prhs[0]);
    if ((mxGetM (prhs[1]) != npix) || (mxGetN (prhs[1]) != 3))
      mexErrMsgTxt ("Input pxvec must be npix by 3");
    pxvec = mxGetPr (prhs[1]);
    if ((mxGetM (prhs[2]) != 1) || (mxGetN (prhs[2]) != 1))
      mexErrMsgTxt ("Input r must be a numeric scalar");
    r = mxGetScalar (prhs[2]);

    /* Allocate output structure */
    maxnout = maxn0;
    mx_neighb_out = mxCreateDoubleMatrix(npix,maxnout,mxREAL);
    neighb_out = mxGetPr(mx_neighb_out);

    /* Initialize temporary pixel list */
    init_pixlist(&pl, maxnout);

    /* Calculate cutoff dot-product corresponding to r */
    dpcut = cos (r * 3.141592653589793238 / 180.0);

    /* Loop over pixels */
    for (j=0; j<npix; j++)
    {
        int last_nused = 0;

        pl.nused = maxn0;
        /* copy nearest neighbors into temp array */
        for (i=0; i<maxn0; i++)
        {
          double tmp = *(((double *)n0) + j + (i*npix));
          /* replace any NaNs with -1 */
          if (!(tmp>=0))
            tmp = -1;
          pl.val[i] = (int)tmp;
        }

        /* repeat expansion up to 100x */
        /* but stop when we're no longer adding pixels */
        for (i=0; i<100; i++)
        {
          int cur_nused = pl.nused;
          /* Loop over pixels already in neighborhood */
          for (k=last_nused; k<cur_nused; k++)
          {
            int tmp = pl.val[k];
            if (!(tmp >= 0) || !(tmp < npix))
              continue;
            /* And over all of their neighbors */
            for (l=0; l<maxn0; l++)
            {
              /* Is it a valid pixel number? */
              int tmp2 = (int)*(((double *)n0) + tmp + (l*npix));
              if ((tmp2 >= 0) && (tmp2 < npix))
              {
                /* Is it within r?  Calculate dot product */
                double dp = 0;
                for (m=0; m<3; m++)
                  dp += *(((double *)pxvec) + j + (m * npix)) * *(((double *)pxvec) + tmp2 + (m * npix));

                /* If so, add it to the neighborhood list. */
                if (dp >= dpcut)
                  add_to_pixlist(&pl, tmp2);
              }
            }
          }
          /* If the last round didn't add any new neighbors, we're done. */
          if (cur_nused == pl.nused)
            break;
          last_nused = cur_nused;
        }
        /* Sort for good measure. */
        sort_pixlist(&pl);

        /* Put temporary neighborhood list into output array. */

        /* First, expand output array if needed */
        if (pl.nused > maxnout)
        {
          printf("Expanding neighb_out from %d to %d.\n", maxnout, pl.nused);
          neighb_out = mxRealloc(neighb_out, npix * pl.nused * sizeof(double));
          /* NaN out newly created entries for previous pixels */
          for (i=0; i<j; i++)
            for (k=maxnout; k<pl.nused; k++)
              *(((double *)neighb_out) + i + (k*npix)) = 0.0 / 0.0;
          maxnout = pl.nused;
          mxSetN(mx_neighb_out,maxnout);
          mxSetPr(mx_neighb_out,neighb_out);
        }

        /* And copy in pixel list */
        for (i=0; i<maxnout; i++)
        {
          if (i < pl.nused)
            *(((double *)neighb_out) + j + (i*npix)) = pl.val[i];
          else
            *(((double *)neighb_out) + j + (i*npix)) = 0.0 / 0.0;
        }
    }

    free_pixlist(&pl);
    plhs[0] = mx_neighb_out;

    return;
}

