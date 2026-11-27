#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlarrk_(integer *n, integer *iw, doublereal *gl, doublereal *gu, doublereal *d__,
            doublereal *e2, doublereal *pivmin, doublereal *reltol, doublereal *w, doublereal *werr,
            integer *info)
{
    integer i__1;
    doublereal d__1, d__2;
    double log(doublereal);
    integer i__, it;
    doublereal mid, eps, tmp1, tmp2, left, atoli, right;
    integer itmax;
    doublereal rtoli, tnorm;
    extern doublereal dlamch_(char *, ftnlen);
    integer negcnt;
    --e2;
    --d__;
    if (*n <= 0) {
        *info = 0;
        return 0;
    }
    eps = dlamch_((char *)"P", (ftnlen)1);
    d__1 = abs(*gl), d__2 = abs(*gu);
    tnorm = max(d__1, d__2);
    rtoli = *reltol;
    atoli = *pivmin * 4.;
    itmax = (integer)((log(tnorm + *pivmin) - log(*pivmin)) / log(2.)) + 2;
    *info = -1;
    left = *gl - tnorm * 2. * eps * *n - *pivmin * 4.;
    right = *gu + tnorm * 2. * eps * *n + *pivmin * 4.;
    it = 0;
L10:
    tmp1 = (d__1 = right - left, abs(d__1));
    d__1 = abs(right), d__2 = abs(left);
    tmp2 = max(d__1, d__2);
    d__1 = max(atoli, *pivmin), d__2 = rtoli * tmp2;
    if (tmp1 < max(d__1, d__2)) {
        *info = 0;
        goto L30;
    }
    if (it > itmax) {
        goto L30;
    }
    ++it;
    mid = (left + right) * .5;
    negcnt = 0;
    tmp1 = d__[1] - mid;
    if (abs(tmp1) < *pivmin) {
        tmp1 = -(*pivmin);
    }
    if (tmp1 <= 0.) {
        ++negcnt;
    }
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
        tmp1 = d__[i__] - e2[i__ - 1] / tmp1 - mid;
        if (abs(tmp1) < *pivmin) {
            tmp1 = -(*pivmin);
        }
        if (tmp1 <= 0.) {
            ++negcnt;
        }
    }
    if (negcnt >= *iw) {
        right = mid;
    } else {
        left = mid;
    }
    goto L10;
L30:
    *w = (left + right) * .5;
    *werr = (d__1 = right - left, abs(d__1)) * .5;
    return 0;
}
#ifdef __cplusplus
}
#endif
