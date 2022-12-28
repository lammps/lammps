#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b6 = -1.;
static integer c__1 = 1;
static doublereal c_b27 = 1.;
int dsygs2_(integer *itype, char *uplo, integer *n, doublereal *a, integer *lda, doublereal *b,
            integer *ldb, integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    doublereal d__1;
    integer k;
    doublereal ct, akk, bkk;
    extern int dsyr2_(char *, integer *, doublereal *, doublereal *, integer *, doublereal *,
                      integer *, doublereal *, integer *, ftnlen),
        dscal_(integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int daxpy_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *);
    logical upper;
    extern int dtrmv_(char *, char *, char *, integer *, doublereal *, integer *, doublereal *,
                      integer *, ftnlen, ftnlen, ftnlen),
        dtrsv_(char *, char *, char *, integer *, doublereal *, integer *, doublereal *, integer *,
               ftnlen, ftnlen, ftnlen),
        xerbla_(char *, integer *, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    *info = 0;
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    if (*itype < 1 || *itype > 3) {
        *info = -1;
    } else if (!upper && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -2;
    } else if (*n < 0) {
        *info = -3;
    } else if (*lda < max(1, *n)) {
        *info = -5;
    } else if (*ldb < max(1, *n)) {
        *info = -7;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DSYGS2", &i__1, (ftnlen)6);
        return 0;
    }
    if (*itype == 1) {
        if (upper) {
            i__1 = *n;
            for (k = 1; k <= i__1; ++k) {
                akk = a[k + k * a_dim1];
                bkk = b[k + k * b_dim1];
                d__1 = bkk;
                akk /= d__1 * d__1;
                a[k + k * a_dim1] = akk;
                if (k < *n) {
                    i__2 = *n - k;
                    d__1 = 1. / bkk;
                    dscal_(&i__2, &d__1, &a[k + (k + 1) * a_dim1], lda);
                    ct = akk * -.5;
                    i__2 = *n - k;
                    daxpy_(&i__2, &ct, &b[k + (k + 1) * b_dim1], ldb, &a[k + (k + 1) * a_dim1],
                           lda);
                    i__2 = *n - k;
                    dsyr2_(uplo, &i__2, &c_b6, &a[k + (k + 1) * a_dim1], lda,
                           &b[k + (k + 1) * b_dim1], ldb, &a[k + 1 + (k + 1) * a_dim1], lda,
                           (ftnlen)1);
                    i__2 = *n - k;
                    daxpy_(&i__2, &ct, &b[k + (k + 1) * b_dim1], ldb, &a[k + (k + 1) * a_dim1],
                           lda);
                    i__2 = *n - k;
                    dtrsv_(uplo, (char *)"Transpose", (char *)"Non-unit", &i__2, &b[k + 1 + (k + 1) * b_dim1], ldb,
                           &a[k + (k + 1) * a_dim1], lda, (ftnlen)1, (ftnlen)9, (ftnlen)8);
                }
            }
        } else {
            i__1 = *n;
            for (k = 1; k <= i__1; ++k) {
                akk = a[k + k * a_dim1];
                bkk = b[k + k * b_dim1];
                d__1 = bkk;
                akk /= d__1 * d__1;
                a[k + k * a_dim1] = akk;
                if (k < *n) {
                    i__2 = *n - k;
                    d__1 = 1. / bkk;
                    dscal_(&i__2, &d__1, &a[k + 1 + k * a_dim1], &c__1);
                    ct = akk * -.5;
                    i__2 = *n - k;
                    daxpy_(&i__2, &ct, &b[k + 1 + k * b_dim1], &c__1, &a[k + 1 + k * a_dim1],
                           &c__1);
                    i__2 = *n - k;
                    dsyr2_(uplo, &i__2, &c_b6, &a[k + 1 + k * a_dim1], &c__1,
                           &b[k + 1 + k * b_dim1], &c__1, &a[k + 1 + (k + 1) * a_dim1], lda,
                           (ftnlen)1);
                    i__2 = *n - k;
                    daxpy_(&i__2, &ct, &b[k + 1 + k * b_dim1], &c__1, &a[k + 1 + k * a_dim1],
                           &c__1);
                    i__2 = *n - k;
                    dtrsv_(uplo, (char *)"No transpose", (char *)"Non-unit", &i__2, &b[k + 1 + (k + 1) * b_dim1],
                           ldb, &a[k + 1 + k * a_dim1], &c__1, (ftnlen)1, (ftnlen)12, (ftnlen)8);
                }
            }
        }
    } else {
        if (upper) {
            i__1 = *n;
            for (k = 1; k <= i__1; ++k) {
                akk = a[k + k * a_dim1];
                bkk = b[k + k * b_dim1];
                i__2 = k - 1;
                dtrmv_(uplo, (char *)"No transpose", (char *)"Non-unit", &i__2, &b[b_offset], ldb,
                       &a[k * a_dim1 + 1], &c__1, (ftnlen)1, (ftnlen)12, (ftnlen)8);
                ct = akk * .5;
                i__2 = k - 1;
                daxpy_(&i__2, &ct, &b[k * b_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
                i__2 = k - 1;
                dsyr2_(uplo, &i__2, &c_b27, &a[k * a_dim1 + 1], &c__1, &b[k * b_dim1 + 1], &c__1,
                       &a[a_offset], lda, (ftnlen)1);
                i__2 = k - 1;
                daxpy_(&i__2, &ct, &b[k * b_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
                i__2 = k - 1;
                dscal_(&i__2, &bkk, &a[k * a_dim1 + 1], &c__1);
                d__1 = bkk;
                a[k + k * a_dim1] = akk * (d__1 * d__1);
            }
        } else {
            i__1 = *n;
            for (k = 1; k <= i__1; ++k) {
                akk = a[k + k * a_dim1];
                bkk = b[k + k * b_dim1];
                i__2 = k - 1;
                dtrmv_(uplo, (char *)"Transpose", (char *)"Non-unit", &i__2, &b[b_offset], ldb, &a[k + a_dim1], lda,
                       (ftnlen)1, (ftnlen)9, (ftnlen)8);
                ct = akk * .5;
                i__2 = k - 1;
                daxpy_(&i__2, &ct, &b[k + b_dim1], ldb, &a[k + a_dim1], lda);
                i__2 = k - 1;
                dsyr2_(uplo, &i__2, &c_b27, &a[k + a_dim1], lda, &b[k + b_dim1], ldb, &a[a_offset],
                       lda, (ftnlen)1);
                i__2 = k - 1;
                daxpy_(&i__2, &ct, &b[k + b_dim1], ldb, &a[k + a_dim1], lda);
                i__2 = k - 1;
                dscal_(&i__2, &bkk, &a[k + a_dim1], lda);
                d__1 = bkk;
                a[k + k * a_dim1] = akk * (d__1 * d__1);
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
