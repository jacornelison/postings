#ifndef GENLEGENDRE_H
#define GENLEGENDRE_H

#ifdef __cplusplus
extern "C" {
#endif

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

legendre_coeff* legendre_gentab(legendre_norm norm, int lmax);
void legendre_free(legendre_coeff* coeff);
void legendre_vector(legendre_coeff* coeff, double* v, int lmax, int m, double x);

#ifdef __cplusplus
}
#endif

#endif /* GENLEGENDRE_H */
