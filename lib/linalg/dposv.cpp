#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dposv_(char *uplo, integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *b,
           integer *ldb, integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int xerbla_(char *, integer *, ftnlen),
        dpotrf_(char *, integer *, doublereal *, integer *, integer *, ftnlen),
        dpotrs_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *,
                integer *, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    *info = 0;
    if (!lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1) && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*nrhs < 0) {
        *info = -3;
    } else if (*lda < max(1, *n)) {
        *info = -5;
    } else if (*ldb < max(1, *n)) {
        *info = -7;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DPOSV ", &i__1, (ftnlen)6);
        return 0;
    }
    dpotrf_(uplo, n, &a[a_offset], lda, info, (ftnlen)1);
    if (*info == 0) {
        dpotrs_(uplo, n, nrhs, &a[a_offset], lda, &b[b_offset], ldb, info, (ftnlen)1);
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
