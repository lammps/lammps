#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlagts_(integer *job, integer *n, doublereal *a, doublereal *b, doublereal *c__,
            doublereal *d__, integer *in, doublereal *y, doublereal *tol, integer *info)
{
    integer i__1;
    doublereal d__1, d__2, d__3, d__4, d__5;
    double d_lmp_sign(doublereal *, doublereal *);
    integer k;
    doublereal ak, eps, temp, pert, absak, sfmin;
    extern doublereal dlamch_(char *, ftnlen);
    extern int xerbla_(char *, integer *, ftnlen);
    doublereal bignum;
    --y;
    --in;
    --d__;
    --c__;
    --b;
    --a;
    *info = 0;
    if (abs(*job) > 2 || *job == 0) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLAGTS", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    eps = dlamch_((char *)"Epsilon", (ftnlen)7);
    sfmin = dlamch_((char *)"Safe minimum", (ftnlen)12);
    bignum = 1. / sfmin;
    if (*job < 0) {
        if (*tol <= 0.) {
            *tol = abs(a[1]);
            if (*n > 1) {
                d__1 = *tol, d__2 = abs(a[2]), d__1 = max(d__1, d__2), d__2 = abs(b[1]);
                *tol = max(d__1, d__2);
            }
            i__1 = *n;
            for (k = 3; k <= i__1; ++k) {
                d__4 = *tol, d__5 = (d__1 = a[k], abs(d__1)), d__4 = max(d__4, d__5),
                d__5 = (d__2 = b[k - 1], abs(d__2)), d__4 = max(d__4, d__5),
                d__5 = (d__3 = d__[k - 2], abs(d__3));
                *tol = max(d__4, d__5);
            }
            *tol *= eps;
            if (*tol == 0.) {
                *tol = eps;
            }
        }
    }
    if (abs(*job) == 1) {
        i__1 = *n;
        for (k = 2; k <= i__1; ++k) {
            if (in[k - 1] == 0) {
                y[k] -= c__[k - 1] * y[k - 1];
            } else {
                temp = y[k - 1];
                y[k - 1] = y[k];
                y[k] = temp - c__[k - 1] * y[k];
            }
        }
        if (*job == 1) {
            for (k = *n; k >= 1; --k) {
                if (k <= *n - 2) {
                    temp = y[k] - b[k] * y[k + 1] - d__[k] * y[k + 2];
                } else if (k == *n - 1) {
                    temp = y[k] - b[k] * y[k + 1];
                } else {
                    temp = y[k];
                }
                ak = a[k];
                absak = abs(ak);
                if (absak < 1.) {
                    if (absak < sfmin) {
                        if (absak == 0. || abs(temp) * sfmin > absak) {
                            *info = k;
                            return 0;
                        } else {
                            temp *= bignum;
                            ak *= bignum;
                        }
                    } else if (abs(temp) > absak * bignum) {
                        *info = k;
                        return 0;
                    }
                }
                y[k] = temp / ak;
            }
        } else {
            for (k = *n; k >= 1; --k) {
                if (k <= *n - 2) {
                    temp = y[k] - b[k] * y[k + 1] - d__[k] * y[k + 2];
                } else if (k == *n - 1) {
                    temp = y[k] - b[k] * y[k + 1];
                } else {
                    temp = y[k];
                }
                ak = a[k];
                pert = d_lmp_sign(tol, &ak);
            L40:
                absak = abs(ak);
                if (absak < 1.) {
                    if (absak < sfmin) {
                        if (absak == 0. || abs(temp) * sfmin > absak) {
                            ak += pert;
                            pert *= 2;
                            goto L40;
                        } else {
                            temp *= bignum;
                            ak *= bignum;
                        }
                    } else if (abs(temp) > absak * bignum) {
                        ak += pert;
                        pert *= 2;
                        goto L40;
                    }
                }
                y[k] = temp / ak;
            }
        }
    } else {
        if (*job == 2) {
            i__1 = *n;
            for (k = 1; k <= i__1; ++k) {
                if (k >= 3) {
                    temp = y[k] - b[k - 1] * y[k - 1] - d__[k - 2] * y[k - 2];
                } else if (k == 2) {
                    temp = y[k] - b[k - 1] * y[k - 1];
                } else {
                    temp = y[k];
                }
                ak = a[k];
                absak = abs(ak);
                if (absak < 1.) {
                    if (absak < sfmin) {
                        if (absak == 0. || abs(temp) * sfmin > absak) {
                            *info = k;
                            return 0;
                        } else {
                            temp *= bignum;
                            ak *= bignum;
                        }
                    } else if (abs(temp) > absak * bignum) {
                        *info = k;
                        return 0;
                    }
                }
                y[k] = temp / ak;
            }
        } else {
            i__1 = *n;
            for (k = 1; k <= i__1; ++k) {
                if (k >= 3) {
                    temp = y[k] - b[k - 1] * y[k - 1] - d__[k - 2] * y[k - 2];
                } else if (k == 2) {
                    temp = y[k] - b[k - 1] * y[k - 1];
                } else {
                    temp = y[k];
                }
                ak = a[k];
                pert = d_lmp_sign(tol, &ak);
            L70:
                absak = abs(ak);
                if (absak < 1.) {
                    if (absak < sfmin) {
                        if (absak == 0. || abs(temp) * sfmin > absak) {
                            ak += pert;
                            pert *= 2;
                            goto L70;
                        } else {
                            temp *= bignum;
                            ak *= bignum;
                        }
                    } else if (abs(temp) > absak * bignum) {
                        ak += pert;
                        pert *= 2;
                        goto L70;
                    }
                }
                y[k] = temp / ak;
            }
        }
        for (k = *n; k >= 2; --k) {
            if (in[k - 1] == 0) {
                y[k - 1] -= c__[k - 1] * y[k];
            } else {
                temp = y[k - 1];
                y[k - 1] = y[k];
                y[k] = temp - c__[k - 1] * y[k];
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
