#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
int dsytf2_(char *uplo, integer *n, doublereal *a, integer *lda, integer *ipiv, integer *info,
            ftnlen uplo_len)
{
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;
    double sqrt(doublereal);
    integer i__, j, k;
    doublereal t, r1, d11, d12, d21, d22;
    integer kk, kp;
    doublereal wk, wkm1, wkp1;
    integer imax, jmax;
    extern int dsyr_(char *, integer *, doublereal *, doublereal *, integer *, doublereal *,
                     integer *, ftnlen);
    doublereal alpha;
    extern int dscal_(integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer kstep;
    logical upper;
    doublereal absakk;
    extern integer idamax_(integer *, doublereal *, integer *);
    extern logical disnan_(doublereal *);
    extern int xerbla_(char *, integer *, ftnlen);
    doublereal colmax, rowmax;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
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
        xerbla_((char *)"DSYTF2", &i__1, (ftnlen)6);
        return 0;
    }
    alpha = (sqrt(17.) + 1.) / 8.;
    if (upper) {
        k = *n;
    L10:
        if (k < 1) {
            goto L70;
        }
        kstep = 1;
        absakk = (d__1 = a[k + k * a_dim1], abs(d__1));
        if (k > 1) {
            i__1 = k - 1;
            imax = idamax_(&i__1, &a[k * a_dim1 + 1], &c__1);
            colmax = (d__1 = a[imax + k * a_dim1], abs(d__1));
        } else {
            colmax = 0.;
        }
        if (max(absakk, colmax) == 0. || disnan_(&absakk)) {
            if (*info == 0) {
                *info = k;
            }
            kp = k;
        } else {
            if (absakk >= alpha * colmax) {
                kp = k;
            } else {
                i__1 = k - imax;
                jmax = imax + idamax_(&i__1, &a[imax + (imax + 1) * a_dim1], lda);
                rowmax = (d__1 = a[imax + jmax * a_dim1], abs(d__1));
                if (imax > 1) {
                    i__1 = imax - 1;
                    jmax = idamax_(&i__1, &a[imax * a_dim1 + 1], &c__1);
                    d__2 = rowmax, d__3 = (d__1 = a[jmax + imax * a_dim1], abs(d__1));
                    rowmax = max(d__2, d__3);
                }
                if (absakk >= alpha * colmax * (colmax / rowmax)) {
                    kp = k;
                } else if ((d__1 = a[imax + imax * a_dim1], abs(d__1)) >= alpha * rowmax) {
                    kp = imax;
                } else {
                    kp = imax;
                    kstep = 2;
                }
            }
            kk = k - kstep + 1;
            if (kp != kk) {
                i__1 = kp - 1;
                dswap_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &c__1);
                i__1 = kk - kp - 1;
                dswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (kp + 1) * a_dim1], lda);
                t = a[kk + kk * a_dim1];
                a[kk + kk * a_dim1] = a[kp + kp * a_dim1];
                a[kp + kp * a_dim1] = t;
                if (kstep == 2) {
                    t = a[k - 1 + k * a_dim1];
                    a[k - 1 + k * a_dim1] = a[kp + k * a_dim1];
                    a[kp + k * a_dim1] = t;
                }
            }
            if (kstep == 1) {
                r1 = 1. / a[k + k * a_dim1];
                i__1 = k - 1;
                d__1 = -r1;
                dsyr_(uplo, &i__1, &d__1, &a[k * a_dim1 + 1], &c__1, &a[a_offset], lda, (ftnlen)1);
                i__1 = k - 1;
                dscal_(&i__1, &r1, &a[k * a_dim1 + 1], &c__1);
            } else {
                if (k > 2) {
                    d12 = a[k - 1 + k * a_dim1];
                    d22 = a[k - 1 + (k - 1) * a_dim1] / d12;
                    d11 = a[k + k * a_dim1] / d12;
                    t = 1. / (d11 * d22 - 1.);
                    d12 = t / d12;
                    for (j = k - 2; j >= 1; --j) {
                        wkm1 = d12 * (d11 * a[j + (k - 1) * a_dim1] - a[j + k * a_dim1]);
                        wk = d12 * (d22 * a[j + k * a_dim1] - a[j + (k - 1) * a_dim1]);
                        for (i__ = j; i__ >= 1; --i__) {
                            a[i__ + j * a_dim1] = a[i__ + j * a_dim1] - a[i__ + k * a_dim1] * wk -
                                                  a[i__ + (k - 1) * a_dim1] * wkm1;
                        }
                        a[j + k * a_dim1] = wk;
                        a[j + (k - 1) * a_dim1] = wkm1;
                    }
                }
            }
        }
        if (kstep == 1) {
            ipiv[k] = kp;
        } else {
            ipiv[k] = -kp;
            ipiv[k - 1] = -kp;
        }
        k -= kstep;
        goto L10;
    } else {
        k = 1;
    L40:
        if (k > *n) {
            goto L70;
        }
        kstep = 1;
        absakk = (d__1 = a[k + k * a_dim1], abs(d__1));
        if (k < *n) {
            i__1 = *n - k;
            imax = k + idamax_(&i__1, &a[k + 1 + k * a_dim1], &c__1);
            colmax = (d__1 = a[imax + k * a_dim1], abs(d__1));
        } else {
            colmax = 0.;
        }
        if (max(absakk, colmax) == 0. || disnan_(&absakk)) {
            if (*info == 0) {
                *info = k;
            }
            kp = k;
        } else {
            if (absakk >= alpha * colmax) {
                kp = k;
            } else {
                i__1 = imax - k;
                jmax = k - 1 + idamax_(&i__1, &a[imax + k * a_dim1], lda);
                rowmax = (d__1 = a[imax + jmax * a_dim1], abs(d__1));
                if (imax < *n) {
                    i__1 = *n - imax;
                    jmax = imax + idamax_(&i__1, &a[imax + 1 + imax * a_dim1], &c__1);
                    d__2 = rowmax, d__3 = (d__1 = a[jmax + imax * a_dim1], abs(d__1));
                    rowmax = max(d__2, d__3);
                }
                if (absakk >= alpha * colmax * (colmax / rowmax)) {
                    kp = k;
                } else if ((d__1 = a[imax + imax * a_dim1], abs(d__1)) >= alpha * rowmax) {
                    kp = imax;
                } else {
                    kp = imax;
                    kstep = 2;
                }
            }
            kk = k + kstep - 1;
            if (kp != kk) {
                if (kp < *n) {
                    i__1 = *n - kp;
                    dswap_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 + kp * a_dim1], &c__1);
                }
                i__1 = kp - kk - 1;
                dswap_(&i__1, &a[kk + 1 + kk * a_dim1], &c__1, &a[kp + (kk + 1) * a_dim1], lda);
                t = a[kk + kk * a_dim1];
                a[kk + kk * a_dim1] = a[kp + kp * a_dim1];
                a[kp + kp * a_dim1] = t;
                if (kstep == 2) {
                    t = a[k + 1 + k * a_dim1];
                    a[k + 1 + k * a_dim1] = a[kp + k * a_dim1];
                    a[kp + k * a_dim1] = t;
                }
            }
            if (kstep == 1) {
                if (k < *n) {
                    d11 = 1. / a[k + k * a_dim1];
                    i__1 = *n - k;
                    d__1 = -d11;
                    dsyr_(uplo, &i__1, &d__1, &a[k + 1 + k * a_dim1], &c__1,
                          &a[k + 1 + (k + 1) * a_dim1], lda, (ftnlen)1);
                    i__1 = *n - k;
                    dscal_(&i__1, &d11, &a[k + 1 + k * a_dim1], &c__1);
                }
            } else {
                if (k < *n - 1) {
                    d21 = a[k + 1 + k * a_dim1];
                    d11 = a[k + 1 + (k + 1) * a_dim1] / d21;
                    d22 = a[k + k * a_dim1] / d21;
                    t = 1. / (d11 * d22 - 1.);
                    d21 = t / d21;
                    i__1 = *n;
                    for (j = k + 2; j <= i__1; ++j) {
                        wk = d21 * (d11 * a[j + k * a_dim1] - a[j + (k + 1) * a_dim1]);
                        wkp1 = d21 * (d22 * a[j + (k + 1) * a_dim1] - a[j + k * a_dim1]);
                        i__2 = *n;
                        for (i__ = j; i__ <= i__2; ++i__) {
                            a[i__ + j * a_dim1] = a[i__ + j * a_dim1] - a[i__ + k * a_dim1] * wk -
                                                  a[i__ + (k + 1) * a_dim1] * wkp1;
                        }
                        a[j + k * a_dim1] = wk;
                        a[j + (k + 1) * a_dim1] = wkp1;
                    }
                }
            }
        }
        if (kstep == 1) {
            ipiv[k] = kp;
        } else {
            ipiv[k] = -kp;
            ipiv[k + 1] = -kp;
        }
        k += kstep;
        goto L40;
    }
L70:
    return 0;
}
#ifdef __cplusplus
}
#endif
