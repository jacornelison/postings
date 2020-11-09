#ifndef SPHERELIB_H
#define SPHERELIB_H

#ifdef __cplusplus
extern "C" {
#endif

/* To disable some optimization tricks, uncomment the following line or
 * pass -DNOTRICKS during compilation.
 */
//#define NOTRICKS

#include <unistd.h>

int distance(ssize_t len, double* dist, double* rvec, ssize_t r0);
int cosdistance(ssize_t len, double* dist, double* rvec, ssize_t r0);
int bearing(ssize_t len, double* ang, double* rvec, ssize_t r0);
int polbearings(ssize_t len, double* c01, double* s01,
        double* c10, double* s10, double* rvec, ssize_t r0);

typedef enum {
    LEGENDRE_NORM_UNIT      = 1,
    LEGENDRE_NORM_SPHERE    = 2
} legendre_norm;

typedef struct {
    int     lmax;
    int     mmax;
    double  P00;
    double* alpha;
    double* beta;
    double* mu;
    double* nu;
} legendre_coeff;

legendre_coeff* legendre_gentab(legendre_norm norm, int lmax, int mmax);
void legendre_free(legendre_coeff** coeff);
void legendrev(legendre_coeff* coeff, double* v, int lmax, int m, double x);

#ifndef NOTRICKS
    // 2018-05-19: holybicep01 has sandybridge cores
    //             Odyssey3 compute nodes provide broadwell cores
    #define TARGET_CLONES __attribute__((\
            target_clones("arch=sandybridge",\
                          "arch=broadwell",\
                          "default")\
            ))
#else
    #define TARGET_CLONES
#endif

#ifdef __cplusplus
}
#endif

#endif /* SPHERELIB_H */

