#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublecomplex c_b1 = {1., 0.};
static integer c__1 = 1;
int zhegs2_(integer *itype, char *uplo, integer *n, doublecomplex *a, integer *lda,
            doublecomplex *b, integer *ldb, integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    doublereal d__1, d__2;
    doublecomplex z__1;
    integer k;
    doublecomplex ct;
    doublereal akk, bkk;
    extern int zher2_(char *, integer *, doublecomplex *, doublecomplex *, integer *,
                      doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    logical upper;
    extern int zaxpy_(integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *,
                      integer *),
        ztrmv_(char *, char *, char *, integer *, doublecomplex *, integer *, doublecomplex *,
               integer *, ftnlen, ftnlen, ftnlen),
        ztrsv_(char *, char *, char *, integer *, doublecomplex *, integer *, doublecomplex *,
               integer *, ftnlen, ftnlen, ftnlen),
        xerbla_(char *, integer *, ftnlen),
        zdscal_(integer *, doublereal *, doublecomplex *, integer *),
        zlacgv_(integer *, doublecomplex *, integer *);
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
        xerbla_((char *)"ZHEGS2", &i__1, (ftnlen)6);
        return 0;
    }
    if (*itype == 1) {
        if (upper) {
            i__1 = *n;
            for (k = 1; k <= i__1; ++k) {
                i__2 = k + k * a_dim1;
                akk = a[i__2].r;
                i__2 = k + k * b_dim1;
                bkk = b[i__2].r;
                d__1 = bkk;
                akk /= d__1 * d__1;
                i__2 = k + k * a_dim1;
                a[i__2].r = akk, a[i__2].i = 0.;
                if (k < *n) {
                    i__2 = *n - k;
                    d__1 = 1. / bkk;
                    zdscal_(&i__2, &d__1, &a[k + (k + 1) * a_dim1], lda);
                    d__1 = akk * -.5;
                    ct.r = d__1, ct.i = 0.;
                    i__2 = *n - k;
                    zlacgv_(&i__2, &a[k + (k + 1) * a_dim1], lda);
                    i__2 = *n - k;
                    zlacgv_(&i__2, &b[k + (k + 1) * b_dim1], ldb);
                    i__2 = *n - k;
                    zaxpy_(&i__2, &ct, &b[k + (k + 1) * b_dim1], ldb, &a[k + (k + 1) * a_dim1],
                           lda);
                    i__2 = *n - k;
                    z__1.r = -1., z__1.i = -0.;
                    zher2_(uplo, &i__2, &z__1, &a[k + (k + 1) * a_dim1], lda,
                           &b[k + (k + 1) * b_dim1], ldb, &a[k + 1 + (k + 1) * a_dim1], lda,
                           (ftnlen)1);
                    i__2 = *n - k;
                    zaxpy_(&i__2, &ct, &b[k + (k + 1) * b_dim1], ldb, &a[k + (k + 1) * a_dim1],
                           lda);
                    i__2 = *n - k;
                    zlacgv_(&i__2, &b[k + (k + 1) * b_dim1], ldb);
                    i__2 = *n - k;
                    ztrsv_(uplo, (char *)"C", (char *)"N", &i__2, &b[k + 1 + (k + 1) * b_dim1], ldb,
                           &a[k + (k + 1) * a_dim1], lda, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                    i__2 = *n - k;
                    zlacgv_(&i__2, &a[k + (k + 1) * a_dim1], lda);
                }
            }
        } else {
            i__1 = *n;
            for (k = 1; k <= i__1; ++k) {
                i__2 = k + k * a_dim1;
                akk = a[i__2].r;
                i__2 = k + k * b_dim1;
                bkk = b[i__2].r;
                d__1 = bkk;
                akk /= d__1 * d__1;
                i__2 = k + k * a_dim1;
                a[i__2].r = akk, a[i__2].i = 0.;
                if (k < *n) {
                    i__2 = *n - k;
                    d__1 = 1. / bkk;
                    zdscal_(&i__2, &d__1, &a[k + 1 + k * a_dim1], &c__1);
                    d__1 = akk * -.5;
                    ct.r = d__1, ct.i = 0.;
                    i__2 = *n - k;
                    zaxpy_(&i__2, &ct, &b[k + 1 + k * b_dim1], &c__1, &a[k + 1 + k * a_dim1],
                           &c__1);
                    i__2 = *n - k;
                    z__1.r = -1., z__1.i = -0.;
                    zher2_(uplo, &i__2, &z__1, &a[k + 1 + k * a_dim1], &c__1,
                           &b[k + 1 + k * b_dim1], &c__1, &a[k + 1 + (k + 1) * a_dim1], lda,
                           (ftnlen)1);
                    i__2 = *n - k;
                    zaxpy_(&i__2, &ct, &b[k + 1 + k * b_dim1], &c__1, &a[k + 1 + k * a_dim1],
                           &c__1);
                    i__2 = *n - k;
                    ztrsv_(uplo, (char *)"N", (char *)"N", &i__2, &b[k + 1 + (k + 1) * b_dim1], ldb,
                           &a[k + 1 + k * a_dim1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                }
            }
        }
    } else {
        if (upper) {
            i__1 = *n;
            for (k = 1; k <= i__1; ++k) {
                i__2 = k + k * a_dim1;
                akk = a[i__2].r;
                i__2 = k + k * b_dim1;
                bkk = b[i__2].r;
                i__2 = k - 1;
                ztrmv_(uplo, (char *)"N", (char *)"N", &i__2, &b[b_offset], ldb, &a[k * a_dim1 + 1], &c__1,
                       (ftnlen)1, (ftnlen)1, (ftnlen)1);
                d__1 = akk * .5;
                ct.r = d__1, ct.i = 0.;
                i__2 = k - 1;
                zaxpy_(&i__2, &ct, &b[k * b_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
                i__2 = k - 1;
                zher2_(uplo, &i__2, &c_b1, &a[k * a_dim1 + 1], &c__1, &b[k * b_dim1 + 1], &c__1,
                       &a[a_offset], lda, (ftnlen)1);
                i__2 = k - 1;
                zaxpy_(&i__2, &ct, &b[k * b_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
                i__2 = k - 1;
                zdscal_(&i__2, &bkk, &a[k * a_dim1 + 1], &c__1);
                i__2 = k + k * a_dim1;
                d__2 = bkk;
                d__1 = akk * (d__2 * d__2);
                a[i__2].r = d__1, a[i__2].i = 0.;
            }
        } else {
            i__1 = *n;
            for (k = 1; k <= i__1; ++k) {
                i__2 = k + k * a_dim1;
                akk = a[i__2].r;
                i__2 = k + k * b_dim1;
                bkk = b[i__2].r;
                i__2 = k - 1;
                zlacgv_(&i__2, &a[k + a_dim1], lda);
                i__2 = k - 1;
                ztrmv_(uplo, (char *)"C", (char *)"N", &i__2, &b[b_offset], ldb, &a[k + a_dim1], lda, (ftnlen)1,
                       (ftnlen)1, (ftnlen)1);
                d__1 = akk * .5;
                ct.r = d__1, ct.i = 0.;
                i__2 = k - 1;
                zlacgv_(&i__2, &b[k + b_dim1], ldb);
                i__2 = k - 1;
                zaxpy_(&i__2, &ct, &b[k + b_dim1], ldb, &a[k + a_dim1], lda);
                i__2 = k - 1;
                zher2_(uplo, &i__2, &c_b1, &a[k + a_dim1], lda, &b[k + b_dim1], ldb, &a[a_offset],
                       lda, (ftnlen)1);
                i__2 = k - 1;
                zaxpy_(&i__2, &ct, &b[k + b_dim1], ldb, &a[k + a_dim1], lda);
                i__2 = k - 1;
                zlacgv_(&i__2, &b[k + b_dim1], ldb);
                i__2 = k - 1;
                zdscal_(&i__2, &bkk, &a[k + a_dim1], lda);
                i__2 = k - 1;
                zlacgv_(&i__2, &a[k + a_dim1], lda);
                i__2 = k + k * a_dim1;
                d__2 = bkk;
                d__1 = akk * (d__2 * d__2);
                a[i__2].r = d__1, a[i__2].i = 0.;
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
