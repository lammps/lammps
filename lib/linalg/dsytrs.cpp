#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b7 = -1.;
static integer c__1 = 1;
static doublereal c_b19 = 1.;
int dsytrs_(char *uplo, integer *n, integer *nrhs, doublereal *a, integer *lda, integer *ipiv,
            doublereal *b, integer *ldb, integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    doublereal d__1;
    integer j, k;
    doublereal ak, bk;
    integer kp;
    doublereal akm1, bkm1;
    extern int dger_(integer *, integer *, doublereal *, doublereal *, integer *, doublereal *,
                     integer *, doublereal *, integer *);
    doublereal akm1k;
    extern int dscal_(integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    doublereal denom;
    extern int dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *,
                      doublereal *, integer *, doublereal *, doublereal *, integer *, ftnlen),
        dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
    logical upper;
    extern int xerbla_(char *, integer *, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    *info = 0;
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    if (!upper && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*nrhs < 0) {
        *info = -3;
    } else if (*lda < max(1, *n)) {
        *info = -5;
    } else if (*ldb < max(1, *n)) {
        *info = -8;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DSYTRS", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0 || *nrhs == 0) {
        return 0;
    }
    if (upper) {
        k = *n;
    L10:
        if (k < 1) {
            goto L30;
        }
        if (ipiv[k] > 0) {
            kp = ipiv[k];
            if (kp != k) {
                dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            i__1 = k - 1;
            dger_(&i__1, nrhs, &c_b7, &a[k * a_dim1 + 1], &c__1, &b[k + b_dim1], ldb,
                  &b[b_dim1 + 1], ldb);
            d__1 = 1. / a[k + k * a_dim1];
            dscal_(nrhs, &d__1, &b[k + b_dim1], ldb);
            --k;
        } else {
            kp = -ipiv[k];
            if (kp != k - 1) {
                dswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            i__1 = k - 2;
            dger_(&i__1, nrhs, &c_b7, &a[k * a_dim1 + 1], &c__1, &b[k + b_dim1], ldb,
                  &b[b_dim1 + 1], ldb);
            i__1 = k - 2;
            dger_(&i__1, nrhs, &c_b7, &a[(k - 1) * a_dim1 + 1], &c__1, &b[k - 1 + b_dim1], ldb,
                  &b[b_dim1 + 1], ldb);
            akm1k = a[k - 1 + k * a_dim1];
            akm1 = a[k - 1 + (k - 1) * a_dim1] / akm1k;
            ak = a[k + k * a_dim1] / akm1k;
            denom = akm1 * ak - 1.;
            i__1 = *nrhs;
            for (j = 1; j <= i__1; ++j) {
                bkm1 = b[k - 1 + j * b_dim1] / akm1k;
                bk = b[k + j * b_dim1] / akm1k;
                b[k - 1 + j * b_dim1] = (ak * bkm1 - bk) / denom;
                b[k + j * b_dim1] = (akm1 * bk - bkm1) / denom;
            }
            k += -2;
        }
        goto L10;
    L30:
        k = 1;
    L40:
        if (k > *n) {
            goto L50;
        }
        if (ipiv[k] > 0) {
            i__1 = k - 1;
            dgemv_((char *)"T", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &a[k * a_dim1 + 1], &c__1, &c_b19,
                   &b[k + b_dim1], ldb, (ftnlen)1);
            kp = ipiv[k];
            if (kp != k) {
                dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            ++k;
        } else {
            i__1 = k - 1;
            dgemv_((char *)"T", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &a[k * a_dim1 + 1], &c__1, &c_b19,
                   &b[k + b_dim1], ldb, (ftnlen)1);
            i__1 = k - 1;
            dgemv_((char *)"T", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &a[(k + 1) * a_dim1 + 1], &c__1,
                   &c_b19, &b[k + 1 + b_dim1], ldb, (ftnlen)1);
            kp = -ipiv[k];
            if (kp != k) {
                dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            k += 2;
        }
        goto L40;
    L50:;
    } else {
        k = 1;
    L60:
        if (k > *n) {
            goto L80;
        }
        if (ipiv[k] > 0) {
            kp = ipiv[k];
            if (kp != k) {
                dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            if (k < *n) {
                i__1 = *n - k;
                dger_(&i__1, nrhs, &c_b7, &a[k + 1 + k * a_dim1], &c__1, &b[k + b_dim1], ldb,
                      &b[k + 1 + b_dim1], ldb);
            }
            d__1 = 1. / a[k + k * a_dim1];
            dscal_(nrhs, &d__1, &b[k + b_dim1], ldb);
            ++k;
        } else {
            kp = -ipiv[k];
            if (kp != k + 1) {
                dswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            if (k < *n - 1) {
                i__1 = *n - k - 1;
                dger_(&i__1, nrhs, &c_b7, &a[k + 2 + k * a_dim1], &c__1, &b[k + b_dim1], ldb,
                      &b[k + 2 + b_dim1], ldb);
                i__1 = *n - k - 1;
                dger_(&i__1, nrhs, &c_b7, &a[k + 2 + (k + 1) * a_dim1], &c__1, &b[k + 1 + b_dim1],
                      ldb, &b[k + 2 + b_dim1], ldb);
            }
            akm1k = a[k + 1 + k * a_dim1];
            akm1 = a[k + k * a_dim1] / akm1k;
            ak = a[k + 1 + (k + 1) * a_dim1] / akm1k;
            denom = akm1 * ak - 1.;
            i__1 = *nrhs;
            for (j = 1; j <= i__1; ++j) {
                bkm1 = b[k + j * b_dim1] / akm1k;
                bk = b[k + 1 + j * b_dim1] / akm1k;
                b[k + j * b_dim1] = (ak * bkm1 - bk) / denom;
                b[k + 1 + j * b_dim1] = (akm1 * bk - bkm1) / denom;
            }
            k += 2;
        }
        goto L60;
    L80:
        k = *n;
    L90:
        if (k < 1) {
            goto L100;
        }
        if (ipiv[k] > 0) {
            if (k < *n) {
                i__1 = *n - k;
                dgemv_((char *)"T", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], ldb, &a[k + 1 + k * a_dim1],
                       &c__1, &c_b19, &b[k + b_dim1], ldb, (ftnlen)1);
            }
            kp = ipiv[k];
            if (kp != k) {
                dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            --k;
        } else {
            if (k < *n) {
                i__1 = *n - k;
                dgemv_((char *)"T", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], ldb, &a[k + 1 + k * a_dim1],
                       &c__1, &c_b19, &b[k + b_dim1], ldb, (ftnlen)1);
                i__1 = *n - k;
                dgemv_((char *)"T", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], ldb,
                       &a[k + 1 + (k - 1) * a_dim1], &c__1, &c_b19, &b[k - 1 + b_dim1], ldb,
                       (ftnlen)1);
            }
            kp = -ipiv[k];
            if (kp != k) {
                dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
            }
            k += -2;
        }
        goto L90;
    L100:;
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
