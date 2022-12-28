#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dgesv_(integer *n, integer *nrhs, doublereal *a, integer *lda, integer *ipiv, doublereal *b,
           integer *ldb, integer *info)
{
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    extern int dgetrf_(integer *, integer *, doublereal *, integer *, integer *, integer *),
        xerbla_(char *, integer *, ftnlen),
        dgetrs_(char *, integer *, integer *, doublereal *, integer *, integer *, doublereal *,
                integer *, integer *, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    *info = 0;
    if (*n < 0) {
        *info = -1;
    } else if (*nrhs < 0) {
        *info = -2;
    } else if (*lda < max(1, *n)) {
        *info = -4;
    } else if (*ldb < max(1, *n)) {
        *info = -7;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DGESV ", &i__1, (ftnlen)6);
        return 0;
    }
    dgetrf_(n, n, &a[a_offset], lda, &ipiv[1], info);
    if (*info == 0) {
        dgetrs_((char *)"No transpose", n, nrhs, &a[a_offset], lda, &ipiv[1], &b[b_offset], ldb, info,
                (ftnlen)12);
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
