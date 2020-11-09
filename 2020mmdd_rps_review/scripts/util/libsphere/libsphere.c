#include "libsphere.h"
#include <math.h>

static inline int approxeq(double a, double b) {
    double m, eps;
    a = fabs(a);
    b = fabs(b);
    m = a > b ? a : b;
    eps = nextafter(m, HUGE_VAL) - m;
    // There are 53 significand bits in double floats, so allowing only the
    // lowest ~1/4 of digits to differ happens at a scale which is
    // ~ eps * 2^floor(53/4). This ends up being 8192.
    return fabs(a - b) < 8192*eps;
}

// Same as what approxeq does, but completely evaluated for b == 1.0
#define ONE_8192ULPS    1.81898940354586e-12
static inline int approxone(double a) {
    return fabs(a - 1.0) < ONE_8192ULPS;
}

static inline double clamp(double x, double l, double h) {
    x = x < l ? l : x;
    x = x > h ? h : x;
    return x;
}

/*
 * Calculates the angular distance (in radians) between two points on a sphere,
 * where the points are given as pairs of unit vectors.
 *
 * ARGUMENTS
 *   len    Length and dimension-1-size of `dist` and `rvec`, respectively.
 *   dist   Output vector of pixel distances.
 *   rvec   Nx3 (column-major) matrix of unit vectors given as (x, y, z)
 *          entries.
 *   r0     Index of particular entry in `rvec` which is source point;
 *          distances of all other pixels are with respect to this point
 *          (including itself).
 */
TARGET_CLONES
int distance(ssize_t len, double* dist, double* rvec, ssize_t r0)
{
    if (r0 > len)
        return -1;

    // Input vector
    double* x = rvec;
    double* y = x + len;
    double* z = y + len;

    double x0 = x[r0];
    double y0 = y[r0];
    double z0 = z[r0];

    for (ssize_t ii = 0; ii < len; ++ii) {
        dist[ii] = acos(clamp(x0*x[ii] + y0*y[ii] + z0*z[ii], -1.0, 1.0));
    }

    return 0;
}


/*
 * Calculates cosine of the angular distance between two points on a sphere,
 * where the points are given as pairs of unit vectors.
 *
 * ARGUMENTS
 *   len    Length and dimension-1-size of `dist` and `rvec`, respectively.
 *   dist   Output vector of pixel distances.
 *   rvec   Nx3 (column-major) matrix of unit vectors given as (x, y, z)
 *          entries.
 *   r0     Index of particular entry in `rvec` which is source point;
 *          distances of all other pixels are with respect to this point
 *          (including itself).
 */
TARGET_CLONES
int cosdistance(ssize_t len, double* dist, double* rvec, ssize_t r0)
{
    if (r0 > len)
        return -1;

    // Input vector
    double* x = rvec;
    double* y = x + len;
    double* z = y + len;

    double x0 = x[r0];
    double y0 = y[r0];
    double z0 = z[r0];

    for (ssize_t ii = 0; ii < len; ++ii) {
        dist[ii] = clamp(x0*x[ii] + y0*y[ii] + z0*z[ii], -1.0, 1.0);
    }

    return 0;
}

TARGET_CLONES
static inline void _bearing2_kernel(
        double x0, double y0, double z0,
        double x1, double y1, double z1,
        double* c, double* s)
{
    // Check if abs-value of dot product is nearly 1, meaning the two
    // vectors are essentially parallel.
    if (approxone(fabs(x0*x1 + y0*y1 + z0*z1))) {
        *c = 1.0;
        *s = 0.0;
        return;
    }
    // Also handle if r0 is parallel to zhat specially.
    if (approxone(fabs(z0))) {
        *c = copysign(1.0, -z0);
        *s = 0.0;
        return;
    }

    // Compute normalized auxiliary vector a1 which is cross product of
    // r0 and r1.
    double a1x = y0*z1 - z0*y1;
    double a1y = z0*x1 - x0*z1;
    double a1z = x0*y1 - y0*x1;
    double a1n = sqrt(a1x*a1x + a1y*a1y + a1z*a1z);
    a1x /= a1n;
    a1y /= a1n;
    a1z /= a1n;

    // Compute another auxiliary vector perpendicular to the Pole, r0 x zhat.
    // This simplifies because we know zhat == (0, 0, 1).
    double a2x = y0;
    double a2y = -x0;
    double a2n = sqrt(a2x*a2x + a2y*a2y);
    a2x /= a2n;
    a2y /= a2n;

    // numerator and denominator in the sense of atan2 arguments.

    // (a1 × a2) ⋅ r0, written out full
    double num = x0 * (/* 0 */ - a1z*a2y)
               + y0 * (a1z*a2x /* - 0 */)
               + z0 * (a1x*a2y - a1y*a2x);
    num = clamp(num, -1.0, 1.0);

    double den = a1x*a2x + a1y*a2y;
    den = clamp(den, -1.0, 1.0);

    *c = den;
    *s = num;
}

/* Calculates the bearing angle between the merdian at rvec[r0] and the
 * great circle connecting it to all other points in rvec (including itself).
 *
 * ARGUMENTS
 *   len    Length and dimension-1-size of `ang` and `rvec`, respectively.
 *   ang    Output vector of bearing angles.
 *   rvec   Nx3 (column-major) matrix of unit vectors given as (x, y, z)
 *          entries.
 *   r0     Index of particular entry in `rvec` which is source point;
 *          angles to all other pixels are with respect to this point
 *          (including itself).
 */
TARGET_CLONES
int bearing(ssize_t len, double* ang, double* rvec, ssize_t r0)
{
    if (r0 > len)
        return -1;

    // Input vector
    double* x = rvec;
    double* y = x + len;
    double* z = y + len;

    double x0 = x[r0];
    double y0 = y[r0];
    double z0 = z[r0];

    for (ssize_t ii = 0; ii < len; ++ii) {
        double c, s;
        _bearing2_kernel(x0, y0, z0, x[ii], y[ii], z[ii], &c, &s);
        ang[ii] = atan2(s, c);
    }

    return 0;
}

/* Calculates the polarization angle components between rvec[r0] and all
 * pairs in rvec (including itself). The polarization angle is twice the
 * bearing angle. `c01` and `s01` are the cosines/sines of the angle at
 * rvec[r0], and `c10` and `s10` are the cosines/sines at the other point.
 *
 * ARGUMENTS
 *   len    Length and dimension-1-size of `c01`/`s01`/`c10`/`s10` and
 *          `rvec`, respectively.
 *   c01    Output vector of polarization cosines at rvec[r0].
 *   s01    Output vector of polarization sines at rvec[r0].
 *   c10    Output vector of polarization cosines at the rvec[r0] complement.
 *   s10    Output vector of polarization sines at the rvec[r0] complement.
 *   rvec   Nx3 (column-major) matrix of unit vectors given as (x, y, z)
 *          entries.
 *   r0     Index of particular entry in `rvec` which is source point;
 *          angles to all other pixels are with respect to this point
 *          (including itself).
 */
TARGET_CLONES
int polbearings(ssize_t len, double* c01, double* s01,
        double* c10, double* s10, double* rvec, ssize_t r0)
{
    if (r0 > len)
        return -1;

    // Input vector
    double* x = rvec;
    double* y = x + len;
    double* z = y + len;

    double x0 = x[r0];
    double y0 = y[r0];
    double z0 = z[r0];

    for (ssize_t ii = 0; ii < len; ++ii) {
        double c, s;

        _bearing2_kernel(x0, y0, z0, x[ii], y[ii], z[ii], &c, &s);
        // Flip signs of both to move from quandrants 3 and 4 back to 1 and 2
        // iff sine is negative.
        if (s < 0) {
            c = -c;
            s = -s;
        }
        c01[ii] = c*c - s*s; // cos(2a) = cos(a)^2 - sin(a)^2
        s01[ii] = 2 * c * s; // sin(2a) = 2 sin(a) cos(a)

        _bearing2_kernel(x[ii], y[ii], z[ii], x0, y0, z0, &c, &s);
        if (s < 0) {
            c = -c;
            s = -s;
        }
        c10[ii] = c*c - s*s;
        s10[ii] = 2 * c * s;
    }

    return 0;
}
