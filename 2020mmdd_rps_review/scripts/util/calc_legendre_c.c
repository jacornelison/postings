/*
 * Matlab-equivalent call sequence is
 *
 *   [P,F10,F12,F22] = calc_legendre(z, lmax, verbose);
 *
 * INPUTS
 *
 *   z          Arguments to Legendre polynomials
 *   lmax       Maximum ell to compute
 *   verbose    non-zero values cause status messages to be printed to stdout
 *
 * OUTPUTS
 *   All outputs are 2D length(z) x (lmax+1) matrices.
 *
 *   P          m=0 Legendre polynomial values
 *   F10        Polarization weighting function (Eqn 16)
 *   F12        Polarization weighting function (Eqn 17)
 *   F22        Polarization weighting function (Eqn 18)
 *
 * 2014-12-12 JBW
 *   Rewritten version of JET's calc_legendre.m for speed.
 *
 *   Implemented changes:
 *     - Cache friendliness: accessing in column-major order in accord with
 *       Matlab's matrix layout.
 *     - No use of temporaries. (Every element is computed from a known kernel
 *       function with no need for temporary values. It's unclear whether
 *       Matlab produces temporary-free code or not, but the speed difference
 *       suggests it doesn't.)
 *
 *   Futher ideas:
 *     - Try adding parallel computations (OpenMP support should be trivial
 *       assuming Matlab isn't allergic to it.)
 *     - Smaller computational blocks to increase cache coherency.
 *
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include "mex.h"
#include "matrix.h"

/* Let the calc_legendre function to be compiled for either C or Fortran
 * storage order. (Matlab uses Fortran --- i.e. column-major --- order, so
 * that is the default if not specified)
 *
 * LEGENDRE_MAJOR_ORDER == 0 indicates C
 * LEGENDRE_MAJOR_ORDER == 1 indicates FORTRAN
 */
#ifndef LEGENDRE_MAJOR_ORDER
#    define LEGENDRE_MAJOR_ORDER 1
#endif

void calc_legendre(int32_t verbose, const double* z, size_t  N, int32_t lmax,
                   double* P, double* P2, double* F10, double* F12, double* F22);

/*
 * [P,F10,F12,F22] = calc_legendre_c(z, lmax, verbose);
 */
void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
    /* Inputs */
    const double* z       = NULL;
    int32_t       lmax    = 0;
    int32_t       verbose = 0;
    /* Outputs */
    double* P    = NULL;
    double* F10  = NULL;
    double* F12  = NULL;
    double* F22  = NULL;

    /* Internal */
    size_t   N    = 0;
    mxArray* mxP2 = NULL;
    double*  P2   = NULL;

    /***** Verify inputs *****/

    if (nrhs != 3)
    {
        mexErrMsgIdAndTxt("matrix:calc_legendre_c:nrhs",
            "Three input arguments required.");
    }
    if (nlhs != 4)
    {
        mexErrMsgIdAndTxt("matrix:calc_legendre_c:nlhs",
            "Four output arguments required.");
    }

    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
    {
        mexErrMsgIdAndTxt("matrix:calc_legendre_c:notDouble",
            "Input z must be of type double.");
    }
    if (mxIsSparse(prhs[0]))
    {
        mexErrMsgIdAndTxt("matrix:calc_legendre_c:isSparse",
            "Input z is sparse. Must be dense.");
    }

    if (!mxIsInt32(prhs[1]) || mxGetNumberOfElements(prhs[1])!=1)
    {
        mexErrMsgIdAndTxt("matrix:calc_legendre_c:notInt32",
            "Input lmax must be scalar of type int32");
    }

    if (!mxIsInt32(prhs[2]) || mxGetNumberOfElements(prhs[2])!=1)
    {
        mexErrMsgIdAndTxt("matrix:calc_legendre_c:notInt32",
            "Input verbose must be scalar of type int32");
    }

    /* Collect data into more convenient forms */
    z = (double*)mxGetData(prhs[0]);
    N = mxGetNumberOfElements(prhs[0]);
    lmax = *((int32_t*)mxGetData(prhs[1]));
    verbose = *((int32_t*)mxGetData(prhs[2]));

    if (lmax < 4)
    {
        mexWarnMsgIdAndTxt("matrix:calc_legendre_c:lmaxMin",
            "lmax < 4. Output will be for lmax = 4.");
        lmax = 4;
    }

    /***** Allocate appropriate sized working matrices *****/

    /* Allocate backing storage */
	mxP2    = mxCreateDoubleMatrix(N, lmax+1, mxREAL);
    plhs[0] = mxCreateDoubleMatrix(N, lmax+1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(N, lmax+1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(N, lmax+1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(N, lmax+1, mxREAL);

    /* Get C pointers to arrays */
    P2  = mxGetData(mxP2);
    P   = mxGetData(plhs[0]);
    F10 = mxGetData(plhs[1]);
    F12 = mxGetData(plhs[2]);
    F22 = mxGetData(plhs[3]);

    calc_legendre(verbose, z, N, lmax, P, P2, F10, F12, F22);

    /* Free the array which is no longer needed */
    mxDestroyArray(mxP2);
    P2 = NULL;
}

/**
 * Calculates Legendre polynomials and three weighting functions in
 * accordance with the recursion relations given in Appendix A of
 * astro-ph/0012120v3, "How to measure CMB polarization power spectra without
 * losing information" (2001) M. Tegmark & A. Oliveira-Costa.
 *
 * INPUTS
 *
 *   verbose    non-zero values cause status messages to be printed to stdout
 *   z          Arguments to Legendre polynomials. Array of length N.
 *   N          Length of z.
 *   lmax       Maximum ell to compute.
 *
 * OUTPUTS
 *   All outputs are 2D matrices and must be pre-allocated arrays (i.e. no
 *   memory is allocated within the function). The form of the matrices
 *   depends on the choice of LEGENDRE_MAJOR_ORDER at compile time:
 *
 *   - LEGENDRE_MAJOR_ORDER == 0  =>  (lmax + 1) by length(z)
 *   - LEGENDRE_MAJOR_ORDER == 1  =>  length(z) by (lmax + 1)
 *
 *   P          m=0 Legendre polynomial values
 *   P2         m=2 Legendre polynomial values
 *   F10        Polarization weighting function (Eqn 16)
 *   F12        Polarization weighting function (Eqn 17)
 *   F22        Polarization weighting function (Eqn 18)
 **/
void calc_legendre(int32_t verbose, const double* z, size_t N, int32_t lmax,
                   double* P, double* P2, double* F10, double* F12, double* F22)
{
    int32_t ll;
    size_t  ii;
    char*   zpone = NULL;
    char*   zmone = NULL;

    #if LEGENDRE_MAJOR_ORDER==1
    #    define IDX(r,c) ((r) + (N)*(c))
    #else
    #    define IDX(r,c) ((r)*(N) + (c))
    #endif

    /***** m=0 Legendre functions (Eqn A19) *****/

    if (verbose)
        printf("Computing m=0 Legendre functions...\n");

    /* Initial conditions */

    /* P_0(z) = 1 */
    for (ii=0; ii<N; ++ii)
        P[IDX(ii,0)] = 1;

    /* P_1(z) = z */
    for (ii=0; ii<N; ++ii)
        P[IDX(ii,1)] = z[ii];

    /* Recursion relation fill ll >= 2 */
    for (ll=2; ll<lmax+1; ++ll)
    {
        int32_t llm1 = 2*ll - 1; /* two-ell minus 1 */
        int32_t lm1 = ll - 1;    /* ell minus 1 */
        int32_t lm2 = ll - 2;    /* ell minus 2 */
        double  ool = 1.0 / ll;  /* one over ell */

        for (ii=0; ii<N; ++ii)
            P[IDX(ii,ll)] = ool * (llm1 * z[ii] * P[IDX(ii,lm1)]
                - lm1 * P[IDX(ii,lm2)]);
    }

    /***** m=2 Legendre functions (Eqn A20) *****/

    if (verbose)
        printf("Computing m=2 Legendre functions...\n");

    /* Initial conditions */

    /* P^2_2(z) = 3(1-z^2) */
    for (ii=0; ii<N; ++ii)
        P2[IDX(ii,2)] = 3 * (1 - z[ii]*z[ii]);

    /* P^2_3(z) = 5z P^2_2(z) */
    for (ii=0; ii<N; ++ii)
        P2[IDX(ii,3)] = 15 * z[ii] * (1 - z[ii]*z[ii]);

    /* Recursion relation fill ll >= 4 */
    for (ll=4; ll<lmax+1; ++ll)
    {
        int32_t llm1 = 2*ll - 1; /* two-ell minus 1 */
        int32_t lm1 = ll - 1;    /* ell minus 1 */
        int32_t lp1 = ll + 1;    /* ell plus 1 */
        int32_t lm2 = ll - 2;    /* ell minus 2 */
        double  oolm2 = 1.0/(ll-2); /* one over (ell minus 2) */

        for (ii=0; ii<N; ++ii)
            P2[IDX(ii,ll)] = oolm2 * (llm1 * z[ii] * P2[IDX(ii,lm1)]
                - lp1 * P2[IDX(ii,lm2)]);
    }

    /*
     * The (1-z^2)~=0 terms may cause answers to blow up. Checking for when
     * z ~= +-1 happens often, so pre-check now and use this to decide in
     * the tighter loops.
     */
    zpone = (char*)mxMalloc(sizeof(char) * N);
    zmone = (char*)mxMalloc(sizeof(char) * N);
    for (ii=0; ii<N; ++ii)
    {
        zpone[ii] = 0;
        zmone[ii] = 0;
        if (fabs(z[ii]-1.0) < 1e-13)
            zpone[ii] = 1;
        else if (fabs(z[ii]+1.0) < 1e-13)
            zmone[ii] = 1;
    }

    /***** F10(z) polarization weight (Eqn A16) *****/

    if (verbose)
        printf("Computing F10 T-Pol weighting functions...\n");

    /* Set all ll==0 and ll==1 cases to zero since den == Inf */
    for (ii=0; ii<N*2; ++ii)
        F10[ii] = 0.0;

    for (ll=2; ll<lmax+1; ++ll)
    {
        int32_t llm1 = ll*(ll-1); /* ell times (ell minus 1) */
        int32_t lm1 = ll - 1;     /* ell minus 1 */
        double  den = 2.0 / sqrt((double)lm1 * ll * (ll+1) * (ll+2));
        for (ii=0; ii<N; ++ii)
        {
            double omzz;

            if (zpone[ii] || zmone[ii])
            {
                F10[IDX(ii,ll)] = 0;
                continue;
            }
            omzz = 1.0 / (1 - z[ii]*z[ii]);
            F10[IDX(ii,ll)] = den * (ll*z[ii]*omzz * P[IDX(ii,lm1)]
                - (ll*omzz + 0.5*llm1) * P[IDX(ii,ll)]);
        }
    }

    /***** F12(z) polarization weight (Eqn A17) *****/

    if (verbose)
        printf("Computing F12 weighting functions...\n");

    /* Set all ll==0 and ll==1 cases to zero since den == Inf */
    for (ii=0; ii<N*2; ++ii)
        F12[ii] = 0.0;

    for (ll=2; ll<lmax+1; ++ll)
    {
        int32_t llm1 = ll*(ll-1); /* ell times (ell minus 1) */
        int32_t lm1 = ll - 1;     /* ell minus 1 */
        int32_t lm4 = ll - 4;     /* ell minus 4 */
        int32_t lp2 = ll + 2;     /* ell plus 2 */
        double  den = 2.0 / ((double)lm1 * ll * (ll+1) * lp2);
        for (ii=0; ii<N; ++ii)
        {
            double omzz;

            if (zpone[ii])
            {
                F12[IDX(ii,ll)] = 0.5;
                continue;
            }
            if (zmone[ii])
            {
                if (ll % 2 == 0)
                    F12[IDX(ii,ll)] = 0.5;
                else
                    F12[IDX(ii,ll)] = -0.5;
                continue;
            }

            omzz = 1.0 / (1 - z[ii]*z[ii]);
            F12[IDX(ii,ll)] = den * (lp2*z[ii]*omzz * P2[IDX(ii,lm1)]
                - (lm4*omzz + 0.5*llm1) * P2[IDX(ii,ll)]);
        }
    }

    /***** F22(z) polarization weight (Eqn A18) *****/

    if (verbose)
        printf("Computing F22 weighting functions...\n");

    /* Set all ll==0 and ll==1 cases to zero since den == Inf */
    for (ii=0; ii<N*2; ++ii)
        F22[ii] = 0.0;

    for (ll=2; ll<lmax+1; ++ll)
    {
        int32_t lm1 = ll - 1;     /* ell minus 1 */
        int32_t lp2 = ll + 2;     /* ell plus 2 */
        double  den = 4.0 / ((double)lm1 * ll * (ll+1) * lp2);
        for (ii=0; ii<N; ++ii)
        {
            double omzz;

            if (zpone[ii])
            {
                F22[IDX(ii,ll)] = -0.5;
                continue;
            }
            if (zmone[ii])
            {
                if (ll % 2 == 0)
                    F22[IDX(ii,ll)] = 0.5;
                else
                    F22[IDX(ii,ll)] = -0.5;
                continue;
            }

            omzz = 1.0 / (1 - z[ii]*z[ii]);
            F22[IDX(ii,ll)] = den * omzz * (lp2 * P2[IDX(ii,lm1)]
                - lm1*z[ii] * P2[IDX(ii,ll)]);
        }
    }

    /* Free the plus/minus flag memory */
    mxFree(zpone);
    mxFree(zmone);

    #undef IDX
}
