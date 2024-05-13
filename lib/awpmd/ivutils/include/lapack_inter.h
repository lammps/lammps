// Interface for LAPACK function

#ifndef LAPACK_INTER_H
#define LAPACK_INTER_H

#include <complex>
typedef int lapack_int;
typedef complex<float> lapack_complex_float;
typedef complex<double> lapack_complex_double;

#define DGETRF dgetrf_
#define DGETRS dgetrs_
#define DGETRI dgetri_
#define ZPPTRF zpptrf_
#define ZPPTRI zpptri_

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
void dgetrf_(const lapack_int *m, const lapack_int *n, double *a,
             const lapack_int *lda, lapack_int *ipiv, lapack_int *info);
void dgetrs_(const char *trans, const lapack_int *n, const lapack_int *nrhs,
             const double *a, const lapack_int *lda, const lapack_int *ipiv,
             double *b, const lapack_int *ldb, lapack_int *info);
void dgetri_(const lapack_int *n, double *a, const lapack_int *lda,
             const lapack_int *ipiv, double *work, const lapack_int *lwork,
             lapack_int *info);
void zpptrf_(const char *uplo, const lapack_int *n, lapack_complex_double *ap,
             lapack_int *info);
void zpptri_(const char *uplo, const lapack_int *n, lapack_complex_double *ap,
             lapack_int *info);
#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* lapack_intER_H */
