#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlarrj_(integer *n, doublereal *d__, doublereal *e2, integer *ifirst, integer *ilast,
            doublereal *rtol, integer *offset, doublereal *w, doublereal *werr, doublereal *work,
            integer *iwork, doublereal *pivmin, doublereal *spdiam, integer *info)
{
    integer i__1, i__2;
    doublereal d__1, d__2;
    double log(doublereal);
    integer i__, j, k, p;
    doublereal s;
    integer i1, i2, ii;
    doublereal fac, mid;
    integer cnt;
    doublereal tmp, left;
    integer iter, nint, prev, next, savi1;
    doublereal right, width, dplus;
    integer olnint, maxitr;
    --iwork;
    --work;
    --werr;
    --w;
    --e2;
    --d__;
    *info = 0;
    if (*n <= 0) {
        return 0;
    }
    maxitr = (integer)((log(*spdiam + *pivmin) - log(*pivmin)) / log(2.)) + 2;
    i1 = *ifirst;
    i2 = *ilast;
    nint = 0;
    prev = 0;
    i__1 = i2;
    for (i__ = i1; i__ <= i__1; ++i__) {
        k = i__ << 1;
        ii = i__ - *offset;
        left = w[ii] - werr[ii];
        mid = w[ii];
        right = w[ii] + werr[ii];
        width = right - mid;
        d__1 = abs(left), d__2 = abs(right);
        tmp = max(d__1, d__2);
        if (width < *rtol * tmp) {
            iwork[k - 1] = -1;
            if (i__ == i1 && i__ < i2) {
                i1 = i__ + 1;
            }
            if (prev >= i1 && i__ <= i2) {
                iwork[(prev << 1) - 1] = i__ + 1;
            }
        } else {
            prev = i__;
            fac = 1.;
        L20:
            cnt = 0;
            s = left;
            dplus = d__[1] - s;
            if (dplus < 0.) {
                ++cnt;
            }
            i__2 = *n;
            for (j = 2; j <= i__2; ++j) {
                dplus = d__[j] - s - e2[j - 1] / dplus;
                if (dplus < 0.) {
                    ++cnt;
                }
            }
            if (cnt > i__ - 1) {
                left -= werr[ii] * fac;
                fac *= 2.;
                goto L20;
            }
            fac = 1.;
        L50:
            cnt = 0;
            s = right;
            dplus = d__[1] - s;
            if (dplus < 0.) {
                ++cnt;
            }
            i__2 = *n;
            for (j = 2; j <= i__2; ++j) {
                dplus = d__[j] - s - e2[j - 1] / dplus;
                if (dplus < 0.) {
                    ++cnt;
                }
            }
            if (cnt < i__) {
                right += werr[ii] * fac;
                fac *= 2.;
                goto L50;
            }
            ++nint;
            iwork[k - 1] = i__ + 1;
            iwork[k] = cnt;
        }
        work[k - 1] = left;
        work[k] = right;
    }
    savi1 = i1;
    iter = 0;
L80:
    prev = i1 - 1;
    i__ = i1;
    olnint = nint;
    i__1 = olnint;
    for (p = 1; p <= i__1; ++p) {
        k = i__ << 1;
        ii = i__ - *offset;
        next = iwork[k - 1];
        left = work[k - 1];
        right = work[k];
        mid = (left + right) * .5;
        width = right - mid;
        d__1 = abs(left), d__2 = abs(right);
        tmp = max(d__1, d__2);
        if (width < *rtol * tmp || iter == maxitr) {
            --nint;
            iwork[k - 1] = 0;
            if (i1 == i__) {
                i1 = next;
            } else {
                if (prev >= i1) {
                    iwork[(prev << 1) - 1] = next;
                }
            }
            i__ = next;
            goto L100;
        }
        prev = i__;
        cnt = 0;
        s = mid;
        dplus = d__[1] - s;
        if (dplus < 0.) {
            ++cnt;
        }
        i__2 = *n;
        for (j = 2; j <= i__2; ++j) {
            dplus = d__[j] - s - e2[j - 1] / dplus;
            if (dplus < 0.) {
                ++cnt;
            }
        }
        if (cnt <= i__ - 1) {
            work[k - 1] = mid;
        } else {
            work[k] = mid;
        }
        i__ = next;
    L100:;
    }
    ++iter;
    if (nint > 0 && iter <= maxitr) {
        goto L80;
    }
    i__1 = *ilast;
    for (i__ = savi1; i__ <= i__1; ++i__) {
        k = i__ << 1;
        ii = i__ - *offset;
        if (iwork[k - 1] == 0) {
            w[ii] = (work[k - 1] + work[k]) * .5;
            werr[ii] = work[k] - w[ii];
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
