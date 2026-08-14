#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c_n1 = -1;
int dsysv_(char *uplo, integer *n, integer *nrhs, doublereal *a, integer *lda, integer *ipiv,
           doublereal *b, integer *ldb, doublereal *work, integer *lwork, integer *info,
           ftnlen uplo_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int xerbla_(char *, integer *, ftnlen),
        dsytrf_(char *, integer *, doublereal *, integer *, integer *, doublereal *, integer *,
                integer *, ftnlen);
    integer lwkopt;
    logical lquery;
    extern int dsytrs_(char *, integer *, integer *, doublereal *, integer *, integer *,
                       doublereal *, integer *, integer *, ftnlen),
        dsytrs2_(char *, integer *, integer *, doublereal *, integer *, integer *, doublereal *,
                 integer *, doublereal *, integer *, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --work;
    *info = 0;
    lquery = *lwork == -1;
    if (!lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1) && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*nrhs < 0) {
        *info = -3;
    } else if (*lda < max(1, *n)) {
        *info = -5;
    } else if (*ldb < max(1, *n)) {
        *info = -8;
    } else if (*lwork < 1 && !lquery) {
        *info = -10;
    }
    if (*info == 0) {
        if (*n == 0) {
            lwkopt = 1;
        } else {
            dsytrf_(uplo, n, &a[a_offset], lda, &ipiv[1], &work[1], &c_n1, info, (ftnlen)1);
            lwkopt = (integer)work[1];
        }
        work[1] = (doublereal)lwkopt;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DSYSV ", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    dsytrf_(uplo, n, &a[a_offset], lda, &ipiv[1], &work[1], lwork, info, (ftnlen)1);
    if (*info == 0) {
        if (*lwork < *n) {
            dsytrs_(uplo, n, nrhs, &a[a_offset], lda, &ipiv[1], &b[b_offset], ldb, info, (ftnlen)1);
        } else {
            dsytrs2_(uplo, n, nrhs, &a[a_offset], lda, &ipiv[1], &b[b_offset], ldb, &work[1], info,
                     (ftnlen)1);
        }
    }
    work[1] = (doublereal)lwkopt;
    return 0;
}
#ifdef __cplusplus
}
#endif
