#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlarrb_(integer *n, doublereal *d__, doublereal *lld, integer *ifirst, integer *ilast,
            doublereal *rtol1, doublereal *rtol2, integer *offset, doublereal *w, doublereal *wgap,
            doublereal *werr, doublereal *work, integer *iwork, doublereal *pivmin,
            doublereal *spdiam, integer *twist, integer *info)
{
    integer i__1;
    doublereal d__1, d__2;
    double log(doublereal);
    integer i__, k, r__, i1, ii, ip;
    doublereal gap, mid, tmp, back, lgap, rgap, left;
    integer iter, nint, prev, next;
    doublereal cvrgd, right, width;
    extern integer dlaneg_(integer *, doublereal *, doublereal *, doublereal *, doublereal *,
                           integer *);
    integer negcnt;
    doublereal mnwdth;
    integer olnint, maxitr;
    --iwork;
    --work;
    --werr;
    --wgap;
    --w;
    --lld;
    --d__;
    *info = 0;
    if (*n <= 0) {
        return 0;
    }
    maxitr = (integer)((log(*spdiam + *pivmin) - log(*pivmin)) / log(2.)) + 2;
    mnwdth = *pivmin * 2.;
    r__ = *twist;
    if (r__ < 1 || r__ > *n) {
        r__ = *n;
    }
    i1 = *ifirst;
    nint = 0;
    prev = 0;
    rgap = wgap[i1 - *offset];
    i__1 = *ilast;
    for (i__ = i1; i__ <= i__1; ++i__) {
        k = i__ << 1;
        ii = i__ - *offset;
        left = w[ii] - werr[ii];
        right = w[ii] + werr[ii];
        lgap = rgap;
        rgap = wgap[ii];
        gap = min(lgap, rgap);
        back = werr[ii];
    L20:
        negcnt = dlaneg_(n, &d__[1], &lld[1], &left, pivmin, &r__);
        if (negcnt > i__ - 1) {
            left -= back;
            back *= 2.;
            goto L20;
        }
        back = werr[ii];
    L50:
        negcnt = dlaneg_(n, &d__[1], &lld[1], &right, pivmin, &r__);
        if (negcnt < i__) {
            right += back;
            back *= 2.;
            goto L50;
        }
        width = (d__1 = left - right, abs(d__1)) * .5;
        d__1 = abs(left), d__2 = abs(right);
        tmp = max(d__1, d__2);
        d__1 = *rtol1 * gap, d__2 = *rtol2 * tmp;
        cvrgd = max(d__1, d__2);
        if (width <= cvrgd || width <= mnwdth) {
            iwork[k - 1] = -1;
            if (i__ == i1 && i__ < *ilast) {
                i1 = i__ + 1;
            }
            if (prev >= i1 && i__ <= *ilast) {
                iwork[(prev << 1) - 1] = i__ + 1;
            }
        } else {
            prev = i__;
            ++nint;
            iwork[k - 1] = i__ + 1;
            iwork[k] = negcnt;
        }
        work[k - 1] = left;
        work[k] = right;
    }
    iter = 0;
L80:
    prev = i1 - 1;
    i__ = i1;
    olnint = nint;
    i__1 = olnint;
    for (ip = 1; ip <= i__1; ++ip) {
        k = i__ << 1;
        ii = i__ - *offset;
        rgap = wgap[ii];
        lgap = rgap;
        if (ii > 1) {
            lgap = wgap[ii - 1];
        }
        gap = min(lgap, rgap);
        next = iwork[k - 1];
        left = work[k - 1];
        right = work[k];
        mid = (left + right) * .5;
        width = right - mid;
        d__1 = abs(left), d__2 = abs(right);
        tmp = max(d__1, d__2);
        d__1 = *rtol1 * gap, d__2 = *rtol2 * tmp;
        cvrgd = max(d__1, d__2);
        if (width <= cvrgd || width <= mnwdth || iter == maxitr) {
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
        negcnt = dlaneg_(n, &d__[1], &lld[1], &mid, pivmin, &r__);
        if (negcnt <= i__ - 1) {
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
    for (i__ = *ifirst; i__ <= i__1; ++i__) {
        k = i__ << 1;
        ii = i__ - *offset;
        if (iwork[k - 1] == 0) {
            w[ii] = (work[k - 1] + work[k]) * .5;
            werr[ii] = work[k] - w[ii];
        }
    }
    i__1 = *ilast;
    for (i__ = *ifirst + 1; i__ <= i__1; ++i__) {
        k = i__ << 1;
        ii = i__ - *offset;
        d__1 = 0., d__2 = w[ii] - werr[ii] - w[ii - 1] - werr[ii - 1];
        wgap[ii - 1] = max(d__1, d__2);
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
