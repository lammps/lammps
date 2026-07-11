#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b33 = 1.;
int dsterf_(integer *n, doublereal *d__, doublereal *e, integer *info)
{
    integer i__1;
    doublereal d__1, d__2, d__3;
    double sqrt(doublereal), d_lmp_sign(doublereal *, doublereal *);
    doublereal c__;
    integer i__, l, m;
    doublereal p, r__, s;
    integer l1;
    doublereal bb, rt1, rt2, eps, rte;
    integer lsv;
    doublereal eps2, oldc;
    integer lend;
    doublereal rmax;
    integer jtot;
    extern int dlae2_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
    doublereal gamma, alpha, sigma, anorm;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, ftnlen);
    integer iscale;
    extern int dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *,
                       integer *, doublereal *, integer *, integer *, ftnlen);
    doublereal oldgam, safmin;
    extern int xerbla_(char *, integer *, ftnlen);
    doublereal safmax;
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *, ftnlen);
    extern int dlasrt_(char *, integer *, doublereal *, integer *, ftnlen);
    integer lendsv;
    doublereal ssfmin;
    integer nmaxit;
    doublereal ssfmax;
    --e;
    --d__;
    *info = 0;
    if (*n < 0) {
        *info = -1;
        i__1 = -(*info);
        xerbla_((char *)"DSTERF", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n <= 1) {
        return 0;
    }
    eps = dlamch_((char *)"E", (ftnlen)1);
    d__1 = eps;
    eps2 = d__1 * d__1;
    safmin = dlamch_((char *)"S", (ftnlen)1);
    safmax = 1. / safmin;
    ssfmax = sqrt(safmax) / 3.;
    ssfmin = sqrt(safmin) / eps2;
    rmax = dlamch_((char *)"O", (ftnlen)1);
    nmaxit = *n * 30;
    sigma = 0.;
    jtot = 0;
    l1 = 1;
L10:
    if (l1 > *n) {
        goto L170;
    }
    if (l1 > 1) {
        e[l1 - 1] = 0.;
    }
    i__1 = *n - 1;
    for (m = l1; m <= i__1; ++m) {
        if ((d__3 = e[m], abs(d__3)) <=
            sqrt((d__1 = d__[m], abs(d__1))) * sqrt((d__2 = d__[m + 1], abs(d__2))) * eps) {
            e[m] = 0.;
            goto L30;
        }
    }
    m = *n;
L30:
    l = l1;
    lsv = l;
    lend = m;
    lendsv = lend;
    l1 = m + 1;
    if (lend == l) {
        goto L10;
    }
    i__1 = lend - l + 1;
    anorm = dlanst_((char *)"M", &i__1, &d__[l], &e[l], (ftnlen)1);
    iscale = 0;
    if (anorm == 0.) {
        goto L10;
    }
    if (anorm > ssfmax) {
        iscale = 1;
        i__1 = lend - l + 1;
        dlascl_((char *)"G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l], n, info, (ftnlen)1);
        i__1 = lend - l;
        dlascl_((char *)"G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n, info, (ftnlen)1);
    } else if (anorm < ssfmin) {
        iscale = 2;
        i__1 = lend - l + 1;
        dlascl_((char *)"G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l], n, info, (ftnlen)1);
        i__1 = lend - l;
        dlascl_((char *)"G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n, info, (ftnlen)1);
    }
    i__1 = lend - 1;
    for (i__ = l; i__ <= i__1; ++i__) {
        d__1 = e[i__];
        e[i__] = d__1 * d__1;
    }
    if ((d__1 = d__[lend], abs(d__1)) < (d__2 = d__[l], abs(d__2))) {
        lend = lsv;
        l = lendsv;
    }
    if (lend >= l) {
    L50:
        if (l != lend) {
            i__1 = lend - 1;
            for (m = l; m <= i__1; ++m) {
                if ((d__2 = e[m], abs(d__2)) <= eps2 * (d__1 = d__[m] * d__[m + 1], abs(d__1))) {
                    goto L70;
                }
            }
        }
        m = lend;
    L70:
        if (m < lend) {
            e[m] = 0.;
        }
        p = d__[l];
        if (m == l) {
            goto L90;
        }
        if (m == l + 1) {
            rte = sqrt(e[l]);
            dlae2_(&d__[l], &rte, &d__[l + 1], &rt1, &rt2);
            d__[l] = rt1;
            d__[l + 1] = rt2;
            e[l] = 0.;
            l += 2;
            if (l <= lend) {
                goto L50;
            }
            goto L150;
        }
        if (jtot == nmaxit) {
            goto L150;
        }
        ++jtot;
        rte = sqrt(e[l]);
        sigma = (d__[l + 1] - p) / (rte * 2.);
        r__ = dlapy2_(&sigma, &c_b33);
        sigma = p - rte / (sigma + d_lmp_sign(&r__, &sigma));
        c__ = 1.;
        s = 0.;
        gamma = d__[m] - sigma;
        p = gamma * gamma;
        i__1 = l;
        for (i__ = m - 1; i__ >= i__1; --i__) {
            bb = e[i__];
            r__ = p + bb;
            if (i__ != m - 1) {
                e[i__ + 1] = s * r__;
            }
            oldc = c__;
            c__ = p / r__;
            s = bb / r__;
            oldgam = gamma;
            alpha = d__[i__];
            gamma = c__ * (alpha - sigma) - s * oldgam;
            d__[i__ + 1] = oldgam + (alpha - gamma);
            if (c__ != 0.) {
                p = gamma * gamma / c__;
            } else {
                p = oldc * bb;
            }
        }
        e[l] = s * p;
        d__[l] = sigma + gamma;
        goto L50;
    L90:
        d__[l] = p;
        ++l;
        if (l <= lend) {
            goto L50;
        }
        goto L150;
    } else {
    L100:
        i__1 = lend + 1;
        for (m = l; m >= i__1; --m) {
            if ((d__2 = e[m - 1], abs(d__2)) <= eps2 * (d__1 = d__[m] * d__[m - 1], abs(d__1))) {
                goto L120;
            }
        }
        m = lend;
    L120:
        if (m > lend) {
            e[m - 1] = 0.;
        }
        p = d__[l];
        if (m == l) {
            goto L140;
        }
        if (m == l - 1) {
            rte = sqrt(e[l - 1]);
            dlae2_(&d__[l], &rte, &d__[l - 1], &rt1, &rt2);
            d__[l] = rt1;
            d__[l - 1] = rt2;
            e[l - 1] = 0.;
            l += -2;
            if (l >= lend) {
                goto L100;
            }
            goto L150;
        }
        if (jtot == nmaxit) {
            goto L150;
        }
        ++jtot;
        rte = sqrt(e[l - 1]);
        sigma = (d__[l - 1] - p) / (rte * 2.);
        r__ = dlapy2_(&sigma, &c_b33);
        sigma = p - rte / (sigma + d_lmp_sign(&r__, &sigma));
        c__ = 1.;
        s = 0.;
        gamma = d__[m] - sigma;
        p = gamma * gamma;
        i__1 = l - 1;
        for (i__ = m; i__ <= i__1; ++i__) {
            bb = e[i__];
            r__ = p + bb;
            if (i__ != m) {
                e[i__ - 1] = s * r__;
            }
            oldc = c__;
            c__ = p / r__;
            s = bb / r__;
            oldgam = gamma;
            alpha = d__[i__ + 1];
            gamma = c__ * (alpha - sigma) - s * oldgam;
            d__[i__] = oldgam + (alpha - gamma);
            if (c__ != 0.) {
                p = gamma * gamma / c__;
            } else {
                p = oldc * bb;
            }
        }
        e[l - 1] = s * p;
        d__[l] = sigma + gamma;
        goto L100;
    L140:
        d__[l] = p;
        --l;
        if (l >= lend) {
            goto L100;
        }
        goto L150;
    }
L150:
    if (iscale == 1) {
        i__1 = lendsv - lsv + 1;
        dlascl_((char *)"G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d__[lsv], n, info, (ftnlen)1);
    }
    if (iscale == 2) {
        i__1 = lendsv - lsv + 1;
        dlascl_((char *)"G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d__[lsv], n, info, (ftnlen)1);
    }
    if (jtot < nmaxit) {
        goto L10;
    }
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
        if (e[i__] != 0.) {
            ++(*info);
        }
    }
    goto L180;
L170:
    dlasrt_((char *)"I", n, &d__[1], info, (ftnlen)1);
L180:
    return 0;
}
#ifdef __cplusplus
}
#endif
