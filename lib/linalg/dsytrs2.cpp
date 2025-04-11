#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b10 = 1.;
int dsytrs2_(char *uplo, integer *n, integer *nrhs, doublereal *a, integer *lda, integer *ipiv,
             doublereal *b, integer *ldb, doublereal *work, integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    doublereal d__1;
    integer i__, j, k;
    doublereal ak, bk;
    integer kp;
    doublereal akm1, bkm1, akm1k;
    extern int dscal_(integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    doublereal denom;
    integer iinfo;
    extern int dswap_(integer *, doublereal *, integer *, doublereal *, integer *),
        dtrsm_(char *, char *, char *, char *, integer *, integer *, doublereal *, doublereal *,
               integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    logical upper;
    extern int xerbla_(char *, integer *, ftnlen),
        dsyconv_(char *, char *, integer *, doublereal *, integer *, integer *, doublereal *,
                 integer *, ftnlen, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --work;
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
        xerbla_((char *)"DSYTRS2", &i__1, (ftnlen)7);
        return 0;
    }
    if (*n == 0 || *nrhs == 0) {
        return 0;
    }
    dsyconv_(uplo, (char *)"C", n, &a[a_offset], lda, &ipiv[1], &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
    if (upper) {
        k = *n;
        while (k >= 1) {
            if (ipiv[k] > 0) {
                kp = ipiv[k];
                if (kp != k) {
                    dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                --k;
            } else {
                kp = -ipiv[k];
                if (kp == -ipiv[k - 1]) {
                    dswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                k += -2;
            }
        }
        dtrsm_((char *)"L", (char *)"U", (char *)"N", (char *)"U", n, nrhs, &c_b10, &a[a_offset], lda, &b[b_offset], ldb, (ftnlen)1,
               (ftnlen)1, (ftnlen)1, (ftnlen)1);
        i__ = *n;
        while (i__ >= 1) {
            if (ipiv[i__] > 0) {
                d__1 = 1. / a[i__ + i__ * a_dim1];
                dscal_(nrhs, &d__1, &b[i__ + b_dim1], ldb);
            } else if (i__ > 1) {
                if (ipiv[i__ - 1] == ipiv[i__]) {
                    akm1k = work[i__];
                    akm1 = a[i__ - 1 + (i__ - 1) * a_dim1] / akm1k;
                    ak = a[i__ + i__ * a_dim1] / akm1k;
                    denom = akm1 * ak - 1.;
                    i__1 = *nrhs;
                    for (j = 1; j <= i__1; ++j) {
                        bkm1 = b[i__ - 1 + j * b_dim1] / akm1k;
                        bk = b[i__ + j * b_dim1] / akm1k;
                        b[i__ - 1 + j * b_dim1] = (ak * bkm1 - bk) / denom;
                        b[i__ + j * b_dim1] = (akm1 * bk - bkm1) / denom;
                    }
                    --i__;
                }
            }
            --i__;
        }
        dtrsm_((char *)"L", (char *)"U", (char *)"T", (char *)"U", n, nrhs, &c_b10, &a[a_offset], lda, &b[b_offset], ldb, (ftnlen)1,
               (ftnlen)1, (ftnlen)1, (ftnlen)1);
        k = 1;
        while (k <= *n) {
            if (ipiv[k] > 0) {
                kp = ipiv[k];
                if (kp != k) {
                    dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                ++k;
            } else {
                kp = -ipiv[k];
                if (k < *n && kp == -ipiv[k + 1]) {
                    dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                k += 2;
            }
        }
    } else {
        k = 1;
        while (k <= *n) {
            if (ipiv[k] > 0) {
                kp = ipiv[k];
                if (kp != k) {
                    dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                ++k;
            } else {
                kp = -ipiv[k + 1];
                if (kp == -ipiv[k]) {
                    dswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                k += 2;
            }
        }
        dtrsm_((char *)"L", (char *)"L", (char *)"N", (char *)"U", n, nrhs, &c_b10, &a[a_offset], lda, &b[b_offset], ldb, (ftnlen)1,
               (ftnlen)1, (ftnlen)1, (ftnlen)1);
        i__ = 1;
        while (i__ <= *n) {
            if (ipiv[i__] > 0) {
                d__1 = 1. / a[i__ + i__ * a_dim1];
                dscal_(nrhs, &d__1, &b[i__ + b_dim1], ldb);
            } else {
                akm1k = work[i__];
                akm1 = a[i__ + i__ * a_dim1] / akm1k;
                ak = a[i__ + 1 + (i__ + 1) * a_dim1] / akm1k;
                denom = akm1 * ak - 1.;
                i__1 = *nrhs;
                for (j = 1; j <= i__1; ++j) {
                    bkm1 = b[i__ + j * b_dim1] / akm1k;
                    bk = b[i__ + 1 + j * b_dim1] / akm1k;
                    b[i__ + j * b_dim1] = (ak * bkm1 - bk) / denom;
                    b[i__ + 1 + j * b_dim1] = (akm1 * bk - bkm1) / denom;
                }
                ++i__;
            }
            ++i__;
        }
        dtrsm_((char *)"L", (char *)"L", (char *)"T", (char *)"U", n, nrhs, &c_b10, &a[a_offset], lda, &b[b_offset], ldb, (ftnlen)1,
               (ftnlen)1, (ftnlen)1, (ftnlen)1);
        k = *n;
        while (k >= 1) {
            if (ipiv[k] > 0) {
                kp = ipiv[k];
                if (kp != k) {
                    dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                --k;
            } else {
                kp = -ipiv[k];
                if (k > 1 && kp == -ipiv[k - 1]) {
                    dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
                }
                k += -2;
            }
        }
    }
    dsyconv_(uplo, (char *)"R", n, &a[a_offset], lda, &ipiv[1], &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
    return 0;
}
#ifdef __cplusplus
}
#endif
