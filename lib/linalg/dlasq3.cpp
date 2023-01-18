#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlasq3_(integer *i0, integer *n0, doublereal *z__, integer *pp, doublereal *dmin__,
            doublereal *sigma, doublereal *desig, doublereal *qmax, integer *nfail, integer *iter,
            integer *ndiv, logical *ieee, integer *ttype, doublereal *dmin1, doublereal *dmin2,
            doublereal *dn, doublereal *dn1, doublereal *dn2, doublereal *g, doublereal *tau)
{
    integer i__1;
    doublereal d__1, d__2;
    double sqrt(doublereal);
    doublereal s, t;
    integer j4, nn;
    doublereal eps, tol;
    integer n0in, ipn4;
    doublereal tol2, temp;
    extern int dlasq4_(integer *, integer *, doublereal *, integer *, integer *, doublereal *,
                       doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                       doublereal *, integer *, doublereal *),
        dlasq5_(integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                logical *, doublereal *),
        dlasq6_(integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                doublereal *, doublereal *, doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern logical disnan_(doublereal *);
    --z__;
    n0in = *n0;
    eps = dlamch_((char *)"Precision", (ftnlen)9);
    tol = eps * 100.;
    d__1 = tol;
    tol2 = d__1 * d__1;
L10:
    if (*n0 < *i0) {
        return 0;
    }
    if (*n0 == *i0) {
        goto L20;
    }
    nn = (*n0 << 2) + *pp;
    if (*n0 == *i0 + 1) {
        goto L40;
    }
    if (z__[nn - 5] > tol2 * (*sigma + z__[nn - 3]) &&
        z__[nn - (*pp << 1) - 4] > tol2 * z__[nn - 7]) {
        goto L30;
    }
L20:
    z__[(*n0 << 2) - 3] = z__[(*n0 << 2) + *pp - 3] + *sigma;
    --(*n0);
    goto L10;
L30:
    if (z__[nn - 9] > tol2 * *sigma && z__[nn - (*pp << 1) - 8] > tol2 * z__[nn - 11]) {
        goto L50;
    }
L40:
    if (z__[nn - 3] > z__[nn - 7]) {
        s = z__[nn - 3];
        z__[nn - 3] = z__[nn - 7];
        z__[nn - 7] = s;
    }
    t = (z__[nn - 7] - z__[nn - 3] + z__[nn - 5]) * .5;
    if (z__[nn - 5] > z__[nn - 3] * tol2 && t != 0.) {
        s = z__[nn - 3] * (z__[nn - 5] / t);
        if (s <= t) {
            s = z__[nn - 3] * (z__[nn - 5] / (t * (sqrt(s / t + 1.) + 1.)));
        } else {
            s = z__[nn - 3] * (z__[nn - 5] / (t + sqrt(t) * sqrt(t + s)));
        }
        t = z__[nn - 7] + (s + z__[nn - 5]);
        z__[nn - 3] *= z__[nn - 7] / t;
        z__[nn - 7] = t;
    }
    z__[(*n0 << 2) - 7] = z__[nn - 7] + *sigma;
    z__[(*n0 << 2) - 3] = z__[nn - 3] + *sigma;
    *n0 += -2;
    goto L10;
L50:
    if (*pp == 2) {
        *pp = 0;
    }
    if (*dmin__ <= 0. || *n0 < n0in) {
        if (z__[(*i0 << 2) + *pp - 3] * 1.5 < z__[(*n0 << 2) + *pp - 3]) {
            ipn4 = *i0 + *n0 << 2;
            i__1 = *i0 + *n0 - 1 << 1;
            for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
                temp = z__[j4 - 3];
                z__[j4 - 3] = z__[ipn4 - j4 - 3];
                z__[ipn4 - j4 - 3] = temp;
                temp = z__[j4 - 2];
                z__[j4 - 2] = z__[ipn4 - j4 - 2];
                z__[ipn4 - j4 - 2] = temp;
                temp = z__[j4 - 1];
                z__[j4 - 1] = z__[ipn4 - j4 - 5];
                z__[ipn4 - j4 - 5] = temp;
                temp = z__[j4];
                z__[j4] = z__[ipn4 - j4 - 4];
                z__[ipn4 - j4 - 4] = temp;
            }
            if (*n0 - *i0 <= 4) {
                z__[(*n0 << 2) + *pp - 1] = z__[(*i0 << 2) + *pp - 1];
                z__[(*n0 << 2) - *pp] = z__[(*i0 << 2) - *pp];
            }
            d__1 = *dmin2, d__2 = z__[(*n0 << 2) + *pp - 1];
            *dmin2 = min(d__1, d__2);
            d__1 = z__[(*n0 << 2) + *pp - 1], d__2 = z__[(*i0 << 2) + *pp - 1],
            d__1 = min(d__1, d__2), d__2 = z__[(*i0 << 2) + *pp + 3];
            z__[(*n0 << 2) + *pp - 1] = min(d__1, d__2);
            d__1 = z__[(*n0 << 2) - *pp], d__2 = z__[(*i0 << 2) - *pp], d__1 = min(d__1, d__2),
            d__2 = z__[(*i0 << 2) - *pp + 4];
            z__[(*n0 << 2) - *pp] = min(d__1, d__2);
            d__1 = *qmax, d__2 = z__[(*i0 << 2) + *pp - 3], d__1 = max(d__1, d__2),
            d__2 = z__[(*i0 << 2) + *pp + 1];
            *qmax = max(d__1, d__2);
            *dmin__ = -0.;
        }
    }
    dlasq4_(i0, n0, &z__[1], pp, &n0in, dmin__, dmin1, dmin2, dn, dn1, dn2, tau, ttype, g);
L70:
    dlasq5_(i0, n0, &z__[1], pp, tau, sigma, dmin__, dmin1, dmin2, dn, dn1, dn2, ieee, &eps);
    *ndiv += *n0 - *i0 + 2;
    ++(*iter);
    if (*dmin__ >= 0. && *dmin1 >= 0.) {
        goto L90;
    } else if (*dmin__ < 0. && *dmin1 > 0. && z__[(*n0 - 1 << 2) - *pp] < tol * (*sigma + *dn1) &&
               abs(*dn) < tol * *sigma) {
        z__[(*n0 - 1 << 2) - *pp + 2] = 0.;
        *dmin__ = 0.;
        goto L90;
    } else if (*dmin__ < 0.) {
        ++(*nfail);
        if (*ttype < -22) {
            *tau = 0.;
        } else if (*dmin1 > 0.) {
            *tau = (*tau + *dmin__) * (1. - eps * 2.);
            *ttype += -11;
        } else {
            *tau *= .25;
            *ttype += -12;
        }
        goto L70;
    } else if (disnan_(dmin__)) {
        if (*tau == 0.) {
            goto L80;
        } else {
            *tau = 0.;
            goto L70;
        }
    } else {
        goto L80;
    }
L80:
    dlasq6_(i0, n0, &z__[1], pp, dmin__, dmin1, dmin2, dn, dn1, dn2);
    *ndiv += *n0 - *i0 + 2;
    ++(*iter);
    *tau = 0.;
L90:
    if (*tau < *sigma) {
        *desig += *tau;
        t = *sigma + *desig;
        *desig -= t - *sigma;
    } else {
        t = *sigma + *tau;
        *desig = *sigma - (t - *tau) + *desig;
    }
    *sigma = t;
    return 0;
}
#ifdef __cplusplus
}
#endif
