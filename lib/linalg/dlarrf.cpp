#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
int dlarrf_(integer *n, doublereal *d__, doublereal *l, doublereal *ld, integer *clstrt,
            integer *clend, doublereal *w, doublereal *wgap, doublereal *werr, doublereal *spdiam,
            doublereal *clgapl, doublereal *clgapr, doublereal *pivmin, doublereal *sigma,
            doublereal *dplus, doublereal *lplus, doublereal *work, integer *info)
{
    integer i__1;
    doublereal d__1, d__2, d__3;
    double sqrt(doublereal);
    integer i__;
    doublereal s, bestshift, smlgrowth, eps, tmp, max1, max2, rrr1, rrr2, znm2, growthbound, fail,
        fact, oldp;
    integer indx;
    doublereal prod;
    integer ktry;
    doublereal fail2, avgap, ldmax, rdmax;
    integer shift;
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    logical dorrr1;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal ldelta;
    logical nofail;
    doublereal mingap, lsigma, rdelta;
    extern logical disnan_(doublereal *);
    logical forcer;
    doublereal rsigma, clwdth;
    logical sawnan1, sawnan2, tryrrr1;
    --work;
    --lplus;
    --dplus;
    --werr;
    --wgap;
    --w;
    --ld;
    --l;
    --d__;
    *info = 0;
    if (*n <= 0) {
        return 0;
    }
    fact = 2.;
    eps = dlamch_((char *)"Precision", (ftnlen)9);
    shift = 0;
    forcer = FALSE_;
    nofail = FALSE_;
    clwdth = (d__1 = w[*clend] - w[*clstrt], abs(d__1)) + werr[*clend] + werr[*clstrt];
    avgap = clwdth / (doublereal)(*clend - *clstrt);
    mingap = min(*clgapl, *clgapr);
    d__1 = w[*clstrt], d__2 = w[*clend];
    lsigma = min(d__1, d__2) - werr[*clstrt];
    d__1 = w[*clstrt], d__2 = w[*clend];
    rsigma = max(d__1, d__2) + werr[*clend];
    lsigma -= abs(lsigma) * 4. * eps;
    rsigma += abs(rsigma) * 4. * eps;
    ldmax = mingap * .25 + *pivmin * 2.;
    rdmax = mingap * .25 + *pivmin * 2.;
    d__1 = avgap, d__2 = wgap[*clstrt];
    ldelta = max(d__1, d__2) / fact;
    d__1 = avgap, d__2 = wgap[*clend - 1];
    rdelta = max(d__1, d__2) / fact;
    s = dlamch_((char *)"S", (ftnlen)1);
    smlgrowth = 1. / s;
    fail = (doublereal)(*n - 1) * mingap / (*spdiam * eps);
    fail2 = (doublereal)(*n - 1) * mingap / (*spdiam * sqrt(eps));
    bestshift = lsigma;
    ktry = 0;
    growthbound = *spdiam * 8.;
L5:
    sawnan1 = FALSE_;
    sawnan2 = FALSE_;
    ldelta = min(ldmax, ldelta);
    rdelta = min(rdmax, rdelta);
    s = -lsigma;
    dplus[1] = d__[1] + s;
    if (abs(dplus[1]) < *pivmin) {
        dplus[1] = -(*pivmin);
        sawnan1 = TRUE_;
    }
    max1 = abs(dplus[1]);
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
        lplus[i__] = ld[i__] / dplus[i__];
        s = s * lplus[i__] * l[i__] - lsigma;
        dplus[i__ + 1] = d__[i__ + 1] + s;
        if ((d__1 = dplus[i__ + 1], abs(d__1)) < *pivmin) {
            dplus[i__ + 1] = -(*pivmin);
            sawnan1 = TRUE_;
        }
        d__2 = max1, d__3 = (d__1 = dplus[i__ + 1], abs(d__1));
        max1 = max(d__2, d__3);
    }
    sawnan1 = sawnan1 || disnan_(&max1);
    if (forcer || max1 <= growthbound && !sawnan1) {
        *sigma = lsigma;
        shift = 1;
        goto L100;
    }
    s = -rsigma;
    work[1] = d__[1] + s;
    if (abs(work[1]) < *pivmin) {
        work[1] = -(*pivmin);
        sawnan2 = TRUE_;
    }
    max2 = abs(work[1]);
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
        work[*n + i__] = ld[i__] / work[i__];
        s = s * work[*n + i__] * l[i__] - rsigma;
        work[i__ + 1] = d__[i__ + 1] + s;
        if ((d__1 = work[i__ + 1], abs(d__1)) < *pivmin) {
            work[i__ + 1] = -(*pivmin);
            sawnan2 = TRUE_;
        }
        d__2 = max2, d__3 = (d__1 = work[i__ + 1], abs(d__1));
        max2 = max(d__2, d__3);
    }
    sawnan2 = sawnan2 || disnan_(&max2);
    if (forcer || max2 <= growthbound && !sawnan2) {
        *sigma = rsigma;
        shift = 2;
        goto L100;
    }
    if (sawnan1 && sawnan2) {
        goto L50;
    } else {
        if (!sawnan1) {
            indx = 1;
            if (max1 <= smlgrowth) {
                smlgrowth = max1;
                bestshift = lsigma;
            }
        }
        if (!sawnan2) {
            if (sawnan1 || max2 <= max1) {
                indx = 2;
            }
            if (max2 <= smlgrowth) {
                smlgrowth = max2;
                bestshift = rsigma;
            }
        }
    }
    if (clwdth < mingap / 128. && min(max1, max2) < fail2 && !sawnan1 && !sawnan2) {
        dorrr1 = TRUE_;
    } else {
        dorrr1 = FALSE_;
    }
    tryrrr1 = TRUE_;
    if (tryrrr1 && dorrr1) {
        if (indx == 1) {
            tmp = (d__1 = dplus[*n], abs(d__1));
            znm2 = 1.;
            prod = 1.;
            oldp = 1.;
            for (i__ = *n - 1; i__ >= 1; --i__) {
                if (prod <= eps) {
                    prod =
                        dplus[i__ + 1] * work[*n + i__ + 1] / (dplus[i__] * work[*n + i__]) * oldp;
                } else {
                    prod *= (d__1 = work[*n + i__], abs(d__1));
                }
                oldp = prod;
                d__1 = prod;
                znm2 += d__1 * d__1;
                d__2 = tmp, d__3 = (d__1 = dplus[i__] * prod, abs(d__1));
                tmp = max(d__2, d__3);
            }
            rrr1 = tmp / (*spdiam * sqrt(znm2));
            if (rrr1 <= 8.) {
                *sigma = lsigma;
                shift = 1;
                goto L100;
            }
        } else if (indx == 2) {
            tmp = (d__1 = work[*n], abs(d__1));
            znm2 = 1.;
            prod = 1.;
            oldp = 1.;
            for (i__ = *n - 1; i__ >= 1; --i__) {
                if (prod <= eps) {
                    prod = work[i__ + 1] * lplus[i__ + 1] / (work[i__] * lplus[i__]) * oldp;
                } else {
                    prod *= (d__1 = lplus[i__], abs(d__1));
                }
                oldp = prod;
                d__1 = prod;
                znm2 += d__1 * d__1;
                d__2 = tmp, d__3 = (d__1 = work[i__] * prod, abs(d__1));
                tmp = max(d__2, d__3);
            }
            rrr2 = tmp / (*spdiam * sqrt(znm2));
            if (rrr2 <= 8.) {
                *sigma = rsigma;
                shift = 2;
                goto L100;
            }
        }
    }
L50:
    if (ktry < 1) {
        d__1 = lsigma - ldelta, d__2 = lsigma - ldmax;
        lsigma = max(d__1, d__2);
        d__1 = rsigma + rdelta, d__2 = rsigma + rdmax;
        rsigma = min(d__1, d__2);
        ldelta *= 2.;
        rdelta *= 2.;
        ++ktry;
        goto L5;
    } else {
        if (smlgrowth < fail || nofail) {
            lsigma = bestshift;
            rsigma = bestshift;
            forcer = TRUE_;
            goto L5;
        } else {
            *info = 1;
            return 0;
        }
    }
L100:
    if (shift == 1) {
    } else if (shift == 2) {
        dcopy_(n, &work[1], &c__1, &dplus[1], &c__1);
        i__1 = *n - 1;
        dcopy_(&i__1, &work[*n + 1], &c__1, &lplus[1], &c__1);
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
