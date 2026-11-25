#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static doublereal c_b11 = -1.;
static doublereal c_b13 = 0.;
int dsytri_(char *uplo, integer *n, doublereal *a, integer *lda, integer *ipiv, doublereal *work,
            integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, i__1;
    doublereal d__1;
    doublereal d__;
    integer k;
    doublereal t, ak;
    integer kp;
    doublereal akp1;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, integer *);
    doublereal temp, akkp1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *),
        dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer kstep;
    logical upper;
    extern int dsymv_(char *, integer *, doublereal *, doublereal *, integer *, doublereal *,
                      integer *, doublereal *, doublereal *, integer *, ftnlen),
        xerbla_(char *, integer *, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    --work;
    *info = 0;
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    if (!upper && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*lda < max(1, *n)) {
        *info = -4;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DSYTRI", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    if (upper) {
        for (*info = *n; *info >= 1; --(*info)) {
            if (ipiv[*info] > 0 && a[*info + *info * a_dim1] == 0.) {
                return 0;
            }
        }
    } else {
        i__1 = *n;
        for (*info = 1; *info <= i__1; ++(*info)) {
            if (ipiv[*info] > 0 && a[*info + *info * a_dim1] == 0.) {
                return 0;
            }
        }
    }
    *info = 0;
    if (upper) {
        k = 1;
    L30:
        if (k > *n) {
            goto L40;
        }
        if (ipiv[k] > 0) {
            a[k + k * a_dim1] = 1. / a[k + k * a_dim1];
            if (k > 1) {
                i__1 = k - 1;
                dcopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
                i__1 = k - 1;
                dsymv_(uplo, &i__1, &c_b11, &a[a_offset], lda, &work[1], &c__1, &c_b13,
                       &a[k * a_dim1 + 1], &c__1, (ftnlen)1);
                i__1 = k - 1;
                a[k + k * a_dim1] -= ddot_(&i__1, &work[1], &c__1, &a[k * a_dim1 + 1], &c__1);
            }
            kstep = 1;
        } else {
            t = (d__1 = a[k + (k + 1) * a_dim1], abs(d__1));
            ak = a[k + k * a_dim1] / t;
            akp1 = a[k + 1 + (k + 1) * a_dim1] / t;
            akkp1 = a[k + (k + 1) * a_dim1] / t;
            d__ = t * (ak * akp1 - 1.);
            a[k + k * a_dim1] = akp1 / d__;
            a[k + 1 + (k + 1) * a_dim1] = ak / d__;
            a[k + (k + 1) * a_dim1] = -akkp1 / d__;
            if (k > 1) {
                i__1 = k - 1;
                dcopy_(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
                i__1 = k - 1;
                dsymv_(uplo, &i__1, &c_b11, &a[a_offset], lda, &work[1], &c__1, &c_b13,
                       &a[k * a_dim1 + 1], &c__1, (ftnlen)1);
                i__1 = k - 1;
                a[k + k * a_dim1] -= ddot_(&i__1, &work[1], &c__1, &a[k * a_dim1 + 1], &c__1);
                i__1 = k - 1;
                a[k + (k + 1) * a_dim1] -=
                    ddot_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[(k + 1) * a_dim1 + 1], &c__1);
                i__1 = k - 1;
                dcopy_(&i__1, &a[(k + 1) * a_dim1 + 1], &c__1, &work[1], &c__1);
                i__1 = k - 1;
                dsymv_(uplo, &i__1, &c_b11, &a[a_offset], lda, &work[1], &c__1, &c_b13,
                       &a[(k + 1) * a_dim1 + 1], &c__1, (ftnlen)1);
                i__1 = k - 1;
                a[k + 1 + (k + 1) * a_dim1] -=
                    ddot_(&i__1, &work[1], &c__1, &a[(k + 1) * a_dim1 + 1], &c__1);
            }
            kstep = 2;
        }
        kp = (i__1 = ipiv[k], abs(i__1));
        if (kp != k) {
            i__1 = kp - 1;
            dswap_(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &c__1);
            i__1 = k - kp - 1;
            dswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + (kp + 1) * a_dim1], lda);
            temp = a[k + k * a_dim1];
            a[k + k * a_dim1] = a[kp + kp * a_dim1];
            a[kp + kp * a_dim1] = temp;
            if (kstep == 2) {
                temp = a[k + (k + 1) * a_dim1];
                a[k + (k + 1) * a_dim1] = a[kp + (k + 1) * a_dim1];
                a[kp + (k + 1) * a_dim1] = temp;
            }
        }
        k += kstep;
        goto L30;
    L40:;
    } else {
        k = *n;
    L50:
        if (k < 1) {
            goto L60;
        }
        if (ipiv[k] > 0) {
            a[k + k * a_dim1] = 1. / a[k + k * a_dim1];
            if (k < *n) {
                i__1 = *n - k;
                dcopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
                i__1 = *n - k;
                dsymv_(uplo, &i__1, &c_b11, &a[k + 1 + (k + 1) * a_dim1], lda, &work[1], &c__1,
                       &c_b13, &a[k + 1 + k * a_dim1], &c__1, (ftnlen)1);
                i__1 = *n - k;
                a[k + k * a_dim1] -= ddot_(&i__1, &work[1], &c__1, &a[k + 1 + k * a_dim1], &c__1);
            }
            kstep = 1;
        } else {
            t = (d__1 = a[k + (k - 1) * a_dim1], abs(d__1));
            ak = a[k - 1 + (k - 1) * a_dim1] / t;
            akp1 = a[k + k * a_dim1] / t;
            akkp1 = a[k + (k - 1) * a_dim1] / t;
            d__ = t * (ak * akp1 - 1.);
            a[k - 1 + (k - 1) * a_dim1] = akp1 / d__;
            a[k + k * a_dim1] = ak / d__;
            a[k + (k - 1) * a_dim1] = -akkp1 / d__;
            if (k < *n) {
                i__1 = *n - k;
                dcopy_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
                i__1 = *n - k;
                dsymv_(uplo, &i__1, &c_b11, &a[k + 1 + (k + 1) * a_dim1], lda, &work[1], &c__1,
                       &c_b13, &a[k + 1 + k * a_dim1], &c__1, (ftnlen)1);
                i__1 = *n - k;
                a[k + k * a_dim1] -= ddot_(&i__1, &work[1], &c__1, &a[k + 1 + k * a_dim1], &c__1);
                i__1 = *n - k;
                a[k + (k - 1) * a_dim1] -= ddot_(&i__1, &a[k + 1 + k * a_dim1], &c__1,
                                                 &a[k + 1 + (k - 1) * a_dim1], &c__1);
                i__1 = *n - k;
                dcopy_(&i__1, &a[k + 1 + (k - 1) * a_dim1], &c__1, &work[1], &c__1);
                i__1 = *n - k;
                dsymv_(uplo, &i__1, &c_b11, &a[k + 1 + (k + 1) * a_dim1], lda, &work[1], &c__1,
                       &c_b13, &a[k + 1 + (k - 1) * a_dim1], &c__1, (ftnlen)1);
                i__1 = *n - k;
                a[k - 1 + (k - 1) * a_dim1] -=
                    ddot_(&i__1, &work[1], &c__1, &a[k + 1 + (k - 1) * a_dim1], &c__1);
            }
            kstep = 2;
        }
        kp = (i__1 = ipiv[k], abs(i__1));
        if (kp != k) {
            if (kp < *n) {
                i__1 = *n - kp;
                dswap_(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + kp * a_dim1], &c__1);
            }
            i__1 = kp - k - 1;
            dswap_(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[kp + (k + 1) * a_dim1], lda);
            temp = a[k + k * a_dim1];
            a[k + k * a_dim1] = a[kp + kp * a_dim1];
            a[kp + kp * a_dim1] = temp;
            if (kstep == 2) {
                temp = a[k + (k - 1) * a_dim1];
                a[k + (k - 1) * a_dim1] = a[kp + (k - 1) * a_dim1];
                a[kp + (k - 1) * a_dim1] = temp;
            }
        }
        k -= kstep;
        goto L50;
    L60:;
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
