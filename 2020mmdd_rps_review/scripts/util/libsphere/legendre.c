#include "libsphere.h"
#include "memory.h"
#include "debug.h"
#include <math.h>

#pragma GCC visibility push(hidden)

MKALLOC(float64, double);
MKALLOC(coeff, legendre_coeff);

double initcond(legendre_norm norm)
{
    switch (norm) {
        case LEGENDRE_NORM_UNIT:
            return 1.0;

        case LEGENDRE_NORM_SPHERE:
            return sqrt(0.25 * M_1_PI);

        default:
            return 0.0;
    }
}

double coeff_alpha(legendre_norm norm, int l, int m)
{
    double coeff = 0.0;
    double ld = l;
    double md = m;
    double fac1, fac2;

    switch (norm) {
        case LEGENDRE_NORM_UNIT:
            coeff = (2*ld - 1.0) / (ld - md);
            break;

        case LEGENDRE_NORM_SPHERE:
            fac1 = (2*ld + 1.0) / ((2*ld - 3.0) * (ld*ld - md*md));
            fac2 = 4 * ((ld-1.0)*(ld-1.0)) - 1.0;
            coeff = sqrt(fac1 * fac2);
            break;
    }

    return coeff;
}

double coeff_beta(legendre_norm norm, int l, int m)
{
    double coeff = 0.0;
    double ld = l;
    double md = m;
    double fac1, fac2;

    switch (norm) {
        case LEGENDRE_NORM_UNIT:
            coeff = (ld + md - 1.0) / (ld - md);
            break;

        case LEGENDRE_NORM_SPHERE:
            fac1 = (2*ld + 1.0) / ((2*ld - 3.0) * (ld*ld - md*md));
            fac2 = (ld-1.0)*(ld-1.0) - md*md;
            coeff = sqrt(fac1 * fac2);
            break;
    }

    return coeff;
}

double coeff_mu(legendre_norm norm, int m)
{
    double coeff = 0.0;
    double md = m;

    switch (norm) {
        case LEGENDRE_NORM_UNIT:
            coeff = 2*md - 1.0;
            break;

        case LEGENDRE_NORM_SPHERE:
            coeff = sqrt(1.0 + 1.0/(2*md));
            break;
    }

    return coeff;
}

double coeff_nu(legendre_norm norm, int m)
{
    double coeff = 0.0;
    double md = m;

    switch (norm) {
        case LEGENDRE_NORM_UNIT:
            coeff = 2*md + 1.0;
            break;

        case LEGENDRE_NORM_SPHERE:
            coeff = sqrt(2*md + 3.0);
            break;
    }

    return coeff;
}

double raise_lm_1term(legendre_coeff* coeff, int m, double y, double plm)
{
    double mu = coeff->mu[m + 1];
    return -mu * y * plm;
}

double raise_l_1term(legendre_coeff* coeff, int m, double x, double plm1m)
{
    double nu = coeff->nu[m];
    return nu * x * plm1m;
}

double raise_l_2term(legendre_coeff* coeff, int l, int m, double x,
                     double plm1m, double plm2m)
{
    int stride = coeff->lmax + 1;
    int elem = (l + 1) + stride*m;
    double alpha = coeff->alpha[elem];
    double beta  = coeff->beta[elem];
    return alpha * x * plm1m - beta * plm2m;
}

#pragma GCC visibility pop

legendre_coeff* legendre_gentab(legendre_norm norm, int lmax, int mmax)
{
    dbglog("called with norm = %i, lmax = %i, mmax = %i\n",
            (int)norm, lmax, mmax);

    legendre_coeff* coeff = coeff_malloc(1, "legendre_gentab: coeff");

    coeff->lmax = lmax;
    coeff->mmax = mmax;
    coeff->P00  = initcond(norm);

    size_t n = (mmax + 1) * (lmax + 1);
    size_t m = lmax + 1;
    coeff->alpha = float64_malloc(n, "legendre_gentab: coeff->alpha");
    coeff->beta  = float64_malloc(n, "legendre_gentab: coeff->beta");
    coeff->mu    = float64_malloc(m, "legendre_gentab: coeff->mu");
    coeff->nu    = float64_malloc(m, "legendre_gentab: coeff->nu");

    for (int mm = 0; mm <= mmax; ++mm)
    {
        coeff->mu[mm] = mm == 0 ? 0.0 : coeff_mu(norm, mm);
        coeff->nu[mm] = coeff_nu(norm, mm);

        for (int ll = mm+1; ll <= lmax; ++ll)
        {
            coeff->alpha[ll + m*mm] = coeff_alpha(norm, ll, mm);
            coeff->beta[ ll + m*mm] = coeff_beta( norm, ll, mm);
        }
    }

    return coeff;
}

void legendre_free(legendre_coeff** coeff)
{
    dbglog("called with coeff = %p, *coeff = %p\n",
            (void*)(coeff), (void*)(*coeff));

    if (!coeff)
        return;
    legendre_coeff* _coeff = *coeff;
    float64_free(&(_coeff->nu));
    float64_free(&(_coeff->mu));
    float64_free(&(_coeff->beta));
    float64_free(&(_coeff->alpha));
    coeff_free(coeff);
}

void legendrev(legendre_coeff* coeff, double* v, int lmax, int m, double x)
{
    dbglog("called with coeff = %p, v = %p, lmax = %i, m = %i, x = %lf\n",
            (void*)coeff, (void*)v, lmax, m, x);

    double pl = coeff->P00;
    double plp1 = 0.0;
    double plm1 = 0.0;

    v[0] = pl;
    if (lmax == 0) {
        return;
    }

    double y = sqrt((1.0 + x) * (1.0 - x));
    // Iterate along the main diagonal until reaching the target m
    for (int n = 0; n < m; ++n)
    {
        plp1 = raise_lm_1term(coeff, n, y, pl);
        pl = plp1;
        v[n] = 0.0;
    }
    v[m] = pl;

    // First step is to boost one in l to P_{m+1}^m using a single-term
    // recurrence
    plp1 = raise_l_1term(coeff, m, x, pl);
    v[m+1] = plp1;
    plm1 = pl;
    pl = plp1;

    // Finish by iterating to P_lmax^m using two-term recurrence
    for (int l = m+1; l < lmax; ++l)
    {
        plp1 = raise_l_2term(coeff, l, m, x, pl, plm1);
        v[l+1] = plp1;
        plm1 = pl;
        pl = plp1;
    }
}
