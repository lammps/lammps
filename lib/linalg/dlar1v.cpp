#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlar1v_(integer *n, integer *b1, integer *bn, doublereal *lambda, doublereal *d__,
            doublereal *l, doublereal *ld, doublereal *lld, doublereal *pivmin, doublereal *gaptol,
            doublereal *z__, logical *wantnc, integer *negcnt, doublereal *ztz, doublereal *mingma,
            integer *r__, integer *isuppz, doublereal *nrminv, doublereal *resid,
            doublereal *rqcorr, doublereal *work)
{
    integer i__1;
    doublereal d__1, d__2, d__3;
    double sqrt(doublereal);
    integer i__;
    doublereal s;
    integer r1, r2;
    doublereal eps, tmp;
    integer neg1, neg2, indp, inds;
    doublereal dplus;
    extern doublereal dlamch_(char *, ftnlen);
    extern logical disnan_(doublereal *);
    integer indlpl, indumn;
    doublereal dminus;
    logical sawnan1, sawnan2;
    --work;
    --isuppz;
    --z__;
    --lld;
    --ld;
    --l;
    --d__;
    eps = dlamch_((char *)"Precision", (ftnlen)9);
    if (*r__ == 0) {
        r1 = *b1;
        r2 = *bn;
    } else {
        r1 = *r__;
        r2 = *r__;
    }
    indlpl = 0;
    indumn = *n;
    inds = (*n << 1) + 1;
    indp = *n * 3 + 1;
    if (*b1 == 1) {
        work[inds] = 0.;
    } else {
        work[inds + *b1 - 1] = lld[*b1 - 1];
    }
    sawnan1 = FALSE_;
    neg1 = 0;
    s = work[inds + *b1 - 1] - *lambda;
    i__1 = r1 - 1;
    for (i__ = *b1; i__ <= i__1; ++i__) {
        dplus = d__[i__] + s;
        work[indlpl + i__] = ld[i__] / dplus;
        if (dplus < 0.) {
            ++neg1;
        }
        work[inds + i__] = s * work[indlpl + i__] * l[i__];
        s = work[inds + i__] - *lambda;
    }
    sawnan1 = disnan_(&s);
    if (sawnan1) {
        goto L60;
    }
    i__1 = r2 - 1;
    for (i__ = r1; i__ <= i__1; ++i__) {
        dplus = d__[i__] + s;
        work[indlpl + i__] = ld[i__] / dplus;
        work[inds + i__] = s * work[indlpl + i__] * l[i__];
        s = work[inds + i__] - *lambda;
    }
    sawnan1 = disnan_(&s);
L60:
    if (sawnan1) {
        neg1 = 0;
        s = work[inds + *b1 - 1] - *lambda;
        i__1 = r1 - 1;
        for (i__ = *b1; i__ <= i__1; ++i__) {
            dplus = d__[i__] + s;
            if (abs(dplus) < *pivmin) {
                dplus = -(*pivmin);
            }
            work[indlpl + i__] = ld[i__] / dplus;
            if (dplus < 0.) {
                ++neg1;
            }
            work[inds + i__] = s * work[indlpl + i__] * l[i__];
            if (work[indlpl + i__] == 0.) {
                work[inds + i__] = lld[i__];
            }
            s = work[inds + i__] - *lambda;
        }
        i__1 = r2 - 1;
        for (i__ = r1; i__ <= i__1; ++i__) {
            dplus = d__[i__] + s;
            if (abs(dplus) < *pivmin) {
                dplus = -(*pivmin);
            }
            work[indlpl + i__] = ld[i__] / dplus;
            work[inds + i__] = s * work[indlpl + i__] * l[i__];
            if (work[indlpl + i__] == 0.) {
                work[inds + i__] = lld[i__];
            }
            s = work[inds + i__] - *lambda;
        }
    }
    sawnan2 = FALSE_;
    neg2 = 0;
    work[indp + *bn - 1] = d__[*bn] - *lambda;
    i__1 = r1;
    for (i__ = *bn - 1; i__ >= i__1; --i__) {
        dminus = lld[i__] + work[indp + i__];
        tmp = d__[i__] / dminus;
        if (dminus < 0.) {
            ++neg2;
        }
        work[indumn + i__] = l[i__] * tmp;
        work[indp + i__ - 1] = work[indp + i__] * tmp - *lambda;
    }
    tmp = work[indp + r1 - 1];
    sawnan2 = disnan_(&tmp);
    if (sawnan2) {
        neg2 = 0;
        i__1 = r1;
        for (i__ = *bn - 1; i__ >= i__1; --i__) {
            dminus = lld[i__] + work[indp + i__];
            if (abs(dminus) < *pivmin) {
                dminus = -(*pivmin);
            }
            tmp = d__[i__] / dminus;
            if (dminus < 0.) {
                ++neg2;
            }
            work[indumn + i__] = l[i__] * tmp;
            work[indp + i__ - 1] = work[indp + i__] * tmp - *lambda;
            if (tmp == 0.) {
                work[indp + i__ - 1] = d__[i__] - *lambda;
            }
        }
    }
    *mingma = work[inds + r1 - 1] + work[indp + r1 - 1];
    if (*mingma < 0.) {
        ++neg1;
    }
    if (*wantnc) {
        *negcnt = neg1 + neg2;
    } else {
        *negcnt = -1;
    }
    if (abs(*mingma) == 0.) {
        *mingma = eps * work[inds + r1 - 1];
    }
    *r__ = r1;
    i__1 = r2 - 1;
    for (i__ = r1; i__ <= i__1; ++i__) {
        tmp = work[inds + i__] + work[indp + i__];
        if (tmp == 0.) {
            tmp = eps * work[inds + i__];
        }
        if (abs(tmp) <= abs(*mingma)) {
            *mingma = tmp;
            *r__ = i__ + 1;
        }
    }
    isuppz[1] = *b1;
    isuppz[2] = *bn;
    z__[*r__] = 1.;
    *ztz = 1.;
    if (!sawnan1 && !sawnan2) {
        i__1 = *b1;
        for (i__ = *r__ - 1; i__ >= i__1; --i__) {
            z__[i__] = -(work[indlpl + i__] * z__[i__ + 1]);
            if (((d__1 = z__[i__], abs(d__1)) + (d__2 = z__[i__ + 1], abs(d__2))) *
                    (d__3 = ld[i__], abs(d__3)) <
                *gaptol) {
                z__[i__] = 0.;
                isuppz[1] = i__ + 1;
                goto L220;
            }
            *ztz += z__[i__] * z__[i__];
        }
    L220:;
    } else {
        i__1 = *b1;
        for (i__ = *r__ - 1; i__ >= i__1; --i__) {
            if (z__[i__ + 1] == 0.) {
                z__[i__] = -(ld[i__ + 1] / ld[i__]) * z__[i__ + 2];
            } else {
                z__[i__] = -(work[indlpl + i__] * z__[i__ + 1]);
            }
            if (((d__1 = z__[i__], abs(d__1)) + (d__2 = z__[i__ + 1], abs(d__2))) *
                    (d__3 = ld[i__], abs(d__3)) <
                *gaptol) {
                z__[i__] = 0.;
                isuppz[1] = i__ + 1;
                goto L240;
            }
            *ztz += z__[i__] * z__[i__];
        }
    L240:;
    }
    if (!sawnan1 && !sawnan2) {
        i__1 = *bn - 1;
        for (i__ = *r__; i__ <= i__1; ++i__) {
            z__[i__ + 1] = -(work[indumn + i__] * z__[i__]);
            if (((d__1 = z__[i__], abs(d__1)) + (d__2 = z__[i__ + 1], abs(d__2))) *
                    (d__3 = ld[i__], abs(d__3)) <
                *gaptol) {
                z__[i__ + 1] = 0.;
                isuppz[2] = i__;
                goto L260;
            }
            *ztz += z__[i__ + 1] * z__[i__ + 1];
        }
    L260:;
    } else {
        i__1 = *bn - 1;
        for (i__ = *r__; i__ <= i__1; ++i__) {
            if (z__[i__] == 0.) {
                z__[i__ + 1] = -(ld[i__ - 1] / ld[i__]) * z__[i__ - 1];
            } else {
                z__[i__ + 1] = -(work[indumn + i__] * z__[i__]);
            }
            if (((d__1 = z__[i__], abs(d__1)) + (d__2 = z__[i__ + 1], abs(d__2))) *
                    (d__3 = ld[i__], abs(d__3)) <
                *gaptol) {
                z__[i__ + 1] = 0.;
                isuppz[2] = i__;
                goto L280;
            }
            *ztz += z__[i__ + 1] * z__[i__ + 1];
        }
    L280:;
    }
    tmp = 1. / *ztz;
    *nrminv = sqrt(tmp);
    *resid = abs(*mingma) * *nrminv;
    *rqcorr = *mingma * tmp;
    return 0;
}
#ifdef __cplusplus
}
#endif
