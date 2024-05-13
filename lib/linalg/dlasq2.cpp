#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__10 = 10;
static integer c__3 = 3;
static integer c__4 = 4;
int dlasq2_(integer *n, doublereal *z__, integer *info)
{
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;
    double sqrt(doublereal);
    doublereal d__, e, g;
    integer k;
    doublereal s, t;
    integer i0, i1, i4, n0, n1;
    doublereal dn;
    integer pp;
    doublereal dn1, dn2, dee, eps, tau, tol;
    integer ipn4;
    doublereal tol2;
    logical ieee;
    integer nbig;
    doublereal dmin__, emin, emax;
    integer kmin, ndiv, iter;
    doublereal qmin, temp, qmax, zmax;
    integer splt;
    doublereal dmin1, dmin2;
    integer nfail;
    doublereal desig, trace, sigma;
    integer iinfo;
    doublereal tempe, tempq;
    integer ttype;
    extern int dlasq3_(integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                       doublereal *, doublereal *, integer *, integer *, integer *, logical *,
                       integer *, doublereal *, doublereal *, doublereal *, doublereal *,
                       doublereal *, doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    doublereal deemin;
    integer iwhila, iwhilb;
    doublereal oldemn, safmin;
    extern int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    extern int dlasrt_(char *, integer *, doublereal *, integer *, ftnlen);
    --z__;
    *info = 0;
    eps = dlamch_((char *)"Precision", (ftnlen)9);
    safmin = dlamch_((char *)"Safe minimum", (ftnlen)12);
    tol = eps * 100.;
    d__1 = tol;
    tol2 = d__1 * d__1;
    if (*n < 0) {
        *info = -1;
        xerbla_((char *)"DLASQ2", &c__1, (ftnlen)6);
        return 0;
    } else if (*n == 0) {
        return 0;
    } else if (*n == 1) {
        if (z__[1] < 0.) {
            *info = -201;
            xerbla_((char *)"DLASQ2", &c__2, (ftnlen)6);
        }
        return 0;
    } else if (*n == 2) {
        if (z__[1] < 0.) {
            *info = -201;
            xerbla_((char *)"DLASQ2", &c__2, (ftnlen)6);
            return 0;
        } else if (z__[2] < 0.) {
            *info = -202;
            xerbla_((char *)"DLASQ2", &c__2, (ftnlen)6);
            return 0;
        } else if (z__[3] < 0.) {
            *info = -203;
            xerbla_((char *)"DLASQ2", &c__2, (ftnlen)6);
            return 0;
        } else if (z__[3] > z__[1]) {
            d__ = z__[3];
            z__[3] = z__[1];
            z__[1] = d__;
        }
        z__[5] = z__[1] + z__[2] + z__[3];
        if (z__[2] > z__[3] * tol2) {
            t = (z__[1] - z__[3] + z__[2]) * .5;
            s = z__[3] * (z__[2] / t);
            if (s <= t) {
                s = z__[3] * (z__[2] / (t * (sqrt(s / t + 1.) + 1.)));
            } else {
                s = z__[3] * (z__[2] / (t + sqrt(t) * sqrt(t + s)));
            }
            t = z__[1] + (s + z__[2]);
            z__[3] *= z__[1] / t;
            z__[1] = t;
        }
        z__[2] = z__[3];
        z__[6] = z__[2] + z__[1];
        return 0;
    }
    z__[*n * 2] = 0.;
    emin = z__[2];
    qmax = 0.;
    zmax = 0.;
    d__ = 0.;
    e = 0.;
    i__1 = *n - 1 << 1;
    for (k = 1; k <= i__1; k += 2) {
        if (z__[k] < 0.) {
            *info = -(k + 200);
            xerbla_((char *)"DLASQ2", &c__2, (ftnlen)6);
            return 0;
        } else if (z__[k + 1] < 0.) {
            *info = -(k + 201);
            xerbla_((char *)"DLASQ2", &c__2, (ftnlen)6);
            return 0;
        }
        d__ += z__[k];
        e += z__[k + 1];
        d__1 = qmax, d__2 = z__[k];
        qmax = max(d__1, d__2);
        d__1 = emin, d__2 = z__[k + 1];
        emin = min(d__1, d__2);
        d__1 = max(qmax, zmax), d__2 = z__[k + 1];
        zmax = max(d__1, d__2);
    }
    if (z__[(*n << 1) - 1] < 0.) {
        *info = -((*n << 1) + 199);
        xerbla_((char *)"DLASQ2", &c__2, (ftnlen)6);
        return 0;
    }
    d__ += z__[(*n << 1) - 1];
    d__1 = qmax, d__2 = z__[(*n << 1) - 1];
    qmax = max(d__1, d__2);
    zmax = max(qmax, zmax);
    if (e == 0.) {
        i__1 = *n;
        for (k = 2; k <= i__1; ++k) {
            z__[k] = z__[(k << 1) - 1];
        }
        dlasrt_((char *)"D", n, &z__[1], &iinfo, (ftnlen)1);
        z__[(*n << 1) - 1] = d__;
        return 0;
    }
    trace = d__ + e;
    if (trace == 0.) {
        z__[(*n << 1) - 1] = 0.;
        return 0;
    }
    ieee = ilaenv_(&c__10, (char *)"DLASQ2", (char *)"N", &c__1, &c__2, &c__3, &c__4, (ftnlen)6, (ftnlen)1) == 1;
    for (k = *n << 1; k >= 2; k += -2) {
        z__[k * 2] = 0.;
        z__[(k << 1) - 1] = z__[k];
        z__[(k << 1) - 2] = 0.;
        z__[(k << 1) - 3] = z__[k - 1];
    }
    i0 = 1;
    n0 = *n;
    if (z__[(i0 << 2) - 3] * 1.5 < z__[(n0 << 2) - 3]) {
        ipn4 = i0 + n0 << 2;
        i__1 = i0 + n0 - 1 << 1;
        for (i4 = i0 << 2; i4 <= i__1; i4 += 4) {
            temp = z__[i4 - 3];
            z__[i4 - 3] = z__[ipn4 - i4 - 3];
            z__[ipn4 - i4 - 3] = temp;
            temp = z__[i4 - 1];
            z__[i4 - 1] = z__[ipn4 - i4 - 5];
            z__[ipn4 - i4 - 5] = temp;
        }
    }
    pp = 0;
    for (k = 1; k <= 2; ++k) {
        d__ = z__[(n0 << 2) + pp - 3];
        i__1 = (i0 << 2) + pp;
        for (i4 = (n0 - 1 << 2) + pp; i4 >= i__1; i4 += -4) {
            if (z__[i4 - 1] <= tol2 * d__) {
                z__[i4 - 1] = -0.;
                d__ = z__[i4 - 3];
            } else {
                d__ = z__[i4 - 3] * (d__ / (d__ + z__[i4 - 1]));
            }
        }
        emin = z__[(i0 << 2) + pp + 1];
        d__ = z__[(i0 << 2) + pp - 3];
        i__1 = (n0 - 1 << 2) + pp;
        for (i4 = (i0 << 2) + pp; i4 <= i__1; i4 += 4) {
            z__[i4 - (pp << 1) - 2] = d__ + z__[i4 - 1];
            if (z__[i4 - 1] <= tol2 * d__) {
                z__[i4 - 1] = -0.;
                z__[i4 - (pp << 1) - 2] = d__;
                z__[i4 - (pp << 1)] = 0.;
                d__ = z__[i4 + 1];
            } else if (safmin * z__[i4 + 1] < z__[i4 - (pp << 1) - 2] &&
                       safmin * z__[i4 - (pp << 1) - 2] < z__[i4 + 1]) {
                temp = z__[i4 + 1] / z__[i4 - (pp << 1) - 2];
                z__[i4 - (pp << 1)] = z__[i4 - 1] * temp;
                d__ *= temp;
            } else {
                z__[i4 - (pp << 1)] = z__[i4 + 1] * (z__[i4 - 1] / z__[i4 - (pp << 1) - 2]);
                d__ = z__[i4 + 1] * (d__ / z__[i4 - (pp << 1) - 2]);
            }
            d__1 = emin, d__2 = z__[i4 - (pp << 1)];
            emin = min(d__1, d__2);
        }
        z__[(n0 << 2) - pp - 2] = d__;
        qmax = z__[(i0 << 2) - pp - 2];
        i__1 = (n0 << 2) - pp - 2;
        for (i4 = (i0 << 2) - pp + 2; i4 <= i__1; i4 += 4) {
            d__1 = qmax, d__2 = z__[i4];
            qmax = max(d__1, d__2);
        }
        pp = 1 - pp;
    }
    ttype = 0;
    dmin1 = 0.;
    dmin2 = 0.;
    dn = 0.;
    dn1 = 0.;
    dn2 = 0.;
    g = 0.;
    tau = 0.;
    iter = 2;
    nfail = 0;
    ndiv = n0 - i0 << 1;
    i__1 = *n + 1;
    for (iwhila = 1; iwhila <= i__1; ++iwhila) {
        if (n0 < 1) {
            goto L170;
        }
        desig = 0.;
        if (n0 == *n) {
            sigma = 0.;
        } else {
            sigma = -z__[(n0 << 2) - 1];
        }
        if (sigma < 0.) {
            *info = 1;
            return 0;
        }
        emax = 0.;
        if (n0 > i0) {
            emin = (d__1 = z__[(n0 << 2) - 5], abs(d__1));
        } else {
            emin = 0.;
        }
        qmin = z__[(n0 << 2) - 3];
        qmax = qmin;
        for (i4 = n0 << 2; i4 >= 8; i4 += -4) {
            if (z__[i4 - 5] <= 0.) {
                goto L100;
            }
            if (qmin >= emax * 4.) {
                d__1 = qmin, d__2 = z__[i4 - 3];
                qmin = min(d__1, d__2);
                d__1 = emax, d__2 = z__[i4 - 5];
                emax = max(d__1, d__2);
            }
            d__1 = qmax, d__2 = z__[i4 - 7] + z__[i4 - 5];
            qmax = max(d__1, d__2);
            d__1 = emin, d__2 = z__[i4 - 5];
            emin = min(d__1, d__2);
        }
        i4 = 4;
    L100:
        i0 = i4 / 4;
        pp = 0;
        if (n0 - i0 > 1) {
            dee = z__[(i0 << 2) - 3];
            deemin = dee;
            kmin = i0;
            i__2 = (n0 << 2) - 3;
            for (i4 = (i0 << 2) + 1; i4 <= i__2; i4 += 4) {
                dee = z__[i4] * (dee / (dee + z__[i4 - 2]));
                if (dee <= deemin) {
                    deemin = dee;
                    kmin = (i4 + 3) / 4;
                }
            }
            if (kmin - i0 << 1 < n0 - kmin && deemin <= z__[(n0 << 2) - 3] * .5) {
                ipn4 = i0 + n0 << 2;
                pp = 2;
                i__2 = i0 + n0 - 1 << 1;
                for (i4 = i0 << 2; i4 <= i__2; i4 += 4) {
                    temp = z__[i4 - 3];
                    z__[i4 - 3] = z__[ipn4 - i4 - 3];
                    z__[ipn4 - i4 - 3] = temp;
                    temp = z__[i4 - 2];
                    z__[i4 - 2] = z__[ipn4 - i4 - 2];
                    z__[ipn4 - i4 - 2] = temp;
                    temp = z__[i4 - 1];
                    z__[i4 - 1] = z__[ipn4 - i4 - 5];
                    z__[ipn4 - i4 - 5] = temp;
                    temp = z__[i4];
                    z__[i4] = z__[ipn4 - i4 - 4];
                    z__[ipn4 - i4 - 4] = temp;
                }
            }
        }
        d__1 = 0., d__2 = qmin - sqrt(qmin) * 2. * sqrt(emax);
        dmin__ = -max(d__1, d__2);
        nbig = (n0 - i0 + 1) * 100;
        i__2 = nbig;
        for (iwhilb = 1; iwhilb <= i__2; ++iwhilb) {
            if (i0 > n0) {
                goto L150;
            }
            dlasq3_(&i0, &n0, &z__[1], &pp, &dmin__, &sigma, &desig, &qmax, &nfail, &iter, &ndiv,
                    &ieee, &ttype, &dmin1, &dmin2, &dn, &dn1, &dn2, &g, &tau);
            pp = 1 - pp;
            if (pp == 0 && n0 - i0 >= 3) {
                if (z__[n0 * 4] <= tol2 * qmax || z__[(n0 << 2) - 1] <= tol2 * sigma) {
                    splt = i0 - 1;
                    qmax = z__[(i0 << 2) - 3];
                    emin = z__[(i0 << 2) - 1];
                    oldemn = z__[i0 * 4];
                    i__3 = n0 - 3 << 2;
                    for (i4 = i0 << 2; i4 <= i__3; i4 += 4) {
                        if (z__[i4] <= tol2 * z__[i4 - 3] || z__[i4 - 1] <= tol2 * sigma) {
                            z__[i4 - 1] = -sigma;
                            splt = i4 / 4;
                            qmax = 0.;
                            emin = z__[i4 + 3];
                            oldemn = z__[i4 + 4];
                        } else {
                            d__1 = qmax, d__2 = z__[i4 + 1];
                            qmax = max(d__1, d__2);
                            d__1 = emin, d__2 = z__[i4 - 1];
                            emin = min(d__1, d__2);
                            d__1 = oldemn, d__2 = z__[i4];
                            oldemn = min(d__1, d__2);
                        }
                    }
                    z__[(n0 << 2) - 1] = emin;
                    z__[n0 * 4] = oldemn;
                    i0 = splt + 1;
                }
            }
        }
        *info = 2;
        i1 = i0;
        n1 = n0;
    L145:
        tempq = z__[(i0 << 2) - 3];
        z__[(i0 << 2) - 3] += sigma;
        i__2 = n0;
        for (k = i0 + 1; k <= i__2; ++k) {
            tempe = z__[(k << 2) - 5];
            z__[(k << 2) - 5] *= tempq / z__[(k << 2) - 7];
            tempq = z__[(k << 2) - 3];
            z__[(k << 2) - 3] = z__[(k << 2) - 3] + sigma + tempe - z__[(k << 2) - 5];
        }
        if (i1 > 1) {
            n1 = i1 - 1;
            while (i1 >= 2 && z__[(i1 << 2) - 5] >= 0.) {
                --i1;
            }
            sigma = -z__[(n1 << 2) - 1];
            goto L145;
        }
        i__2 = *n;
        for (k = 1; k <= i__2; ++k) {
            z__[(k << 1) - 1] = z__[(k << 2) - 3];
            if (k < n0) {
                z__[k * 2] = z__[(k << 2) - 1];
            } else {
                z__[k * 2] = 0.;
            }
        }
        return 0;
    L150:;
    }
    *info = 3;
    return 0;
L170:
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
        z__[k] = z__[(k << 2) - 3];
    }
    dlasrt_((char *)"D", n, &z__[1], &iinfo, (ftnlen)1);
    e = 0.;
    for (k = *n; k >= 1; --k) {
        e += z__[k];
    }
    z__[(*n << 1) + 1] = trace;
    z__[(*n << 1) + 2] = e;
    z__[(*n << 1) + 3] = (doublereal)iter;
    i__1 = *n;
    z__[(*n << 1) + 4] = (doublereal)ndiv / (doublereal)(i__1 * i__1);
    z__[(*n << 1) + 5] = nfail * 100. / (doublereal)iter;
    return 0;
}
#ifdef __cplusplus
}
#endif
