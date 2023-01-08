#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b15 = -.125;
static integer c__1 = 1;
static doublereal c_b49 = 1.;
static doublereal c_b72 = -1.;
int dbdsqr_(char *uplo, integer *n, integer *ncvt, integer *nru, integer *ncc, doublereal *d__,
            doublereal *e, doublereal *vt, integer *ldvt, doublereal *u, integer *ldu,
            doublereal *c__, integer *ldc, doublereal *work, integer *info, ftnlen uplo_len)
{
    integer c_dim1, c_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;
    double pow_lmp_dd(doublereal *, doublereal *), sqrt(doublereal), d_lmp_sign(doublereal *, doublereal *);
    integer iterdivn;
    doublereal f, g, h__;
    integer i__, j, m;
    doublereal r__;
    integer maxitdivn;
    doublereal cs;
    integer ll;
    doublereal sn, mu;
    integer nm1, nm12, nm13, lll;
    doublereal eps, sll, tol, abse;
    integer idir;
    doublereal abss;
    integer oldm;
    doublereal cosl;
    integer isub, iter;
    doublereal unfl, sinl, cosr, smin, smax, sinr;
    extern int drot_(integer *, doublereal *, integer *, doublereal *, integer *, doublereal *,
                     doublereal *),
        dlas2_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *),
        dscal_(integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    doublereal oldcs;
    extern int dlasr_(char *, char *, char *, integer *, integer *, doublereal *, doublereal *,
                      doublereal *, integer *, ftnlen, ftnlen, ftnlen);
    integer oldll;
    doublereal shift, sigmn, oldsn;
    extern int dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
    doublereal sminl, sigmx;
    logical lower;
    extern int dlasq1_(integer *, doublereal *, doublereal *, doublereal *, integer *),
        dlasv2_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                doublereal *, doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern int dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *),
        xerbla_(char *, integer *, ftnlen);
    doublereal sminoa, thresh;
    logical rotate;
    doublereal tolmul;
    --d__;
    --e;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    *info = 0;
    lower = lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1);
    if (!lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1) && !lower) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*ncvt < 0) {
        *info = -3;
    } else if (*nru < 0) {
        *info = -4;
    } else if (*ncc < 0) {
        *info = -5;
    } else if (*ncvt == 0 && *ldvt < 1 || *ncvt > 0 && *ldvt < max(1, *n)) {
        *info = -9;
    } else if (*ldu < max(1, *nru)) {
        *info = -11;
    } else if (*ncc == 0 && *ldc < 1 || *ncc > 0 && *ldc < max(1, *n)) {
        *info = -13;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DBDSQR", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    if (*n == 1) {
        goto L160;
    }
    rotate = *ncvt > 0 || *nru > 0 || *ncc > 0;
    if (!rotate) {
        dlasq1_(n, &d__[1], &e[1], &work[1], info);
        if (*info != 2) {
            return 0;
        }
        *info = 0;
    }
    nm1 = *n - 1;
    nm12 = nm1 + nm1;
    nm13 = nm12 + nm1;
    idir = 0;
    eps = dlamch_((char *)"Epsilon", (ftnlen)7);
    unfl = dlamch_((char *)"Safe minimum", (ftnlen)12);
    if (lower) {
        i__1 = *n - 1;
        for (i__ = 1; i__ <= i__1; ++i__) {
            dlartg_(&d__[i__], &e[i__], &cs, &sn, &r__);
            d__[i__] = r__;
            e[i__] = sn * d__[i__ + 1];
            d__[i__ + 1] = cs * d__[i__ + 1];
            work[i__] = cs;
            work[nm1 + i__] = sn;
        }
        if (*nru > 0) {
            dlasr_((char *)"R", (char *)"V", (char *)"F", nru, n, &work[1], &work[*n], &u[u_offset], ldu, (ftnlen)1,
                   (ftnlen)1, (ftnlen)1);
        }
        if (*ncc > 0) {
            dlasr_((char *)"L", (char *)"V", (char *)"F", n, ncc, &work[1], &work[*n], &c__[c_offset], ldc, (ftnlen)1,
                   (ftnlen)1, (ftnlen)1);
        }
    }
    d__3 = 100., d__4 = pow_lmp_dd(&eps, &c_b15);
    d__1 = 10., d__2 = min(d__3, d__4);
    tolmul = max(d__1, d__2);
    tol = tolmul * eps;
    smax = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        d__2 = smax, d__3 = (d__1 = d__[i__], abs(d__1));
        smax = max(d__2, d__3);
    }
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
        d__2 = smax, d__3 = (d__1 = e[i__], abs(d__1));
        smax = max(d__2, d__3);
    }
    sminl = 0.;
    if (tol >= 0.) {
        sminoa = abs(d__[1]);
        if (sminoa == 0.) {
            goto L50;
        }
        mu = sminoa;
        i__1 = *n;
        for (i__ = 2; i__ <= i__1; ++i__) {
            mu = (d__2 = d__[i__], abs(d__2)) * (mu / (mu + (d__1 = e[i__ - 1], abs(d__1))));
            sminoa = min(sminoa, mu);
            if (sminoa == 0.) {
                goto L50;
            }
        }
    L50:
        sminoa /= sqrt((doublereal)(*n));
        d__1 = tol * sminoa, d__2 = *n * (*n * unfl) * 6;
        thresh = max(d__1, d__2);
    } else {
        d__1 = abs(tol) * smax, d__2 = *n * (*n * unfl) * 6;
        thresh = max(d__1, d__2);
    }
    maxitdivn = *n * 6;
    iterdivn = 0;
    iter = -1;
    oldll = -1;
    oldm = -1;
    m = *n;
L60:
    if (m <= 1) {
        goto L160;
    }
    if (iter >= *n) {
        iter -= *n;
        ++iterdivn;
        if (iterdivn >= maxitdivn) {
            goto L200;
        }
    }
    if (tol < 0. && (d__1 = d__[m], abs(d__1)) <= thresh) {
        d__[m] = 0.;
    }
    smax = (d__1 = d__[m], abs(d__1));
    smin = smax;
    i__1 = m - 1;
    for (lll = 1; lll <= i__1; ++lll) {
        ll = m - lll;
        abss = (d__1 = d__[ll], abs(d__1));
        abse = (d__1 = e[ll], abs(d__1));
        if (tol < 0. && abss <= thresh) {
            d__[ll] = 0.;
        }
        if (abse <= thresh) {
            goto L80;
        }
        smin = min(smin, abss);
        d__1 = max(smax, abss);
        smax = max(d__1, abse);
    }
    ll = 0;
    goto L90;
L80:
    e[ll] = 0.;
    if (ll == m - 1) {
        --m;
        goto L60;
    }
L90:
    ++ll;
    if (ll == m - 1) {
        dlasv2_(&d__[m - 1], &e[m - 1], &d__[m], &sigmn, &sigmx, &sinr, &cosr, &sinl, &cosl);
        d__[m - 1] = sigmx;
        e[m - 1] = 0.;
        d__[m] = sigmn;
        if (*ncvt > 0) {
            drot_(ncvt, &vt[m - 1 + vt_dim1], ldvt, &vt[m + vt_dim1], ldvt, &cosr, &sinr);
        }
        if (*nru > 0) {
            drot_(nru, &u[(m - 1) * u_dim1 + 1], &c__1, &u[m * u_dim1 + 1], &c__1, &cosl, &sinl);
        }
        if (*ncc > 0) {
            drot_(ncc, &c__[m - 1 + c_dim1], ldc, &c__[m + c_dim1], ldc, &cosl, &sinl);
        }
        m += -2;
        goto L60;
    }
    if (ll > oldm || m < oldll) {
        if ((d__1 = d__[ll], abs(d__1)) >= (d__2 = d__[m], abs(d__2))) {
            idir = 1;
        } else {
            idir = 2;
        }
    }
    if (idir == 1) {
        if ((d__2 = e[m - 1], abs(d__2)) <= abs(tol) * (d__1 = d__[m], abs(d__1)) ||
            tol < 0. && (d__3 = e[m - 1], abs(d__3)) <= thresh) {
            e[m - 1] = 0.;
            goto L60;
        }
        if (tol >= 0.) {
            mu = (d__1 = d__[ll], abs(d__1));
            sminl = mu;
            i__1 = m - 1;
            for (lll = ll; lll <= i__1; ++lll) {
                if ((d__1 = e[lll], abs(d__1)) <= tol * mu) {
                    e[lll] = 0.;
                    goto L60;
                }
                mu = (d__2 = d__[lll + 1], abs(d__2)) * (mu / (mu + (d__1 = e[lll], abs(d__1))));
                sminl = min(sminl, mu);
            }
        }
    } else {
        if ((d__2 = e[ll], abs(d__2)) <= abs(tol) * (d__1 = d__[ll], abs(d__1)) ||
            tol < 0. && (d__3 = e[ll], abs(d__3)) <= thresh) {
            e[ll] = 0.;
            goto L60;
        }
        if (tol >= 0.) {
            mu = (d__1 = d__[m], abs(d__1));
            sminl = mu;
            i__1 = ll;
            for (lll = m - 1; lll >= i__1; --lll) {
                if ((d__1 = e[lll], abs(d__1)) <= tol * mu) {
                    e[lll] = 0.;
                    goto L60;
                }
                mu = (d__2 = d__[lll], abs(d__2)) * (mu / (mu + (d__1 = e[lll], abs(d__1))));
                sminl = min(sminl, mu);
            }
        }
    }
    oldll = ll;
    oldm = m;
    d__1 = eps, d__2 = tol * .01;
    if (tol >= 0. && *n * tol * (sminl / smax) <= max(d__1, d__2)) {
        shift = 0.;
    } else {
        if (idir == 1) {
            sll = (d__1 = d__[ll], abs(d__1));
            dlas2_(&d__[m - 1], &e[m - 1], &d__[m], &shift, &r__);
        } else {
            sll = (d__1 = d__[m], abs(d__1));
            dlas2_(&d__[ll], &e[ll], &d__[ll + 1], &shift, &r__);
        }
        if (sll > 0.) {
            d__1 = shift / sll;
            if (d__1 * d__1 < eps) {
                shift = 0.;
            }
        }
    }
    iter = iter + m - ll;
    if (shift == 0.) {
        if (idir == 1) {
            cs = 1.;
            oldcs = 1.;
            i__1 = m - 1;
            for (i__ = ll; i__ <= i__1; ++i__) {
                d__1 = d__[i__] * cs;
                dlartg_(&d__1, &e[i__], &cs, &sn, &r__);
                if (i__ > ll) {
                    e[i__ - 1] = oldsn * r__;
                }
                d__1 = oldcs * r__;
                d__2 = d__[i__ + 1] * sn;
                dlartg_(&d__1, &d__2, &oldcs, &oldsn, &d__[i__]);
                work[i__ - ll + 1] = cs;
                work[i__ - ll + 1 + nm1] = sn;
                work[i__ - ll + 1 + nm12] = oldcs;
                work[i__ - ll + 1 + nm13] = oldsn;
            }
            h__ = d__[m] * cs;
            d__[m] = h__ * oldcs;
            e[m - 1] = h__ * oldsn;
            if (*ncvt > 0) {
                i__1 = m - ll + 1;
                dlasr_((char *)"L", (char *)"V", (char *)"F", &i__1, ncvt, &work[1], &work[*n], &vt[ll + vt_dim1], ldvt,
                       (ftnlen)1, (ftnlen)1, (ftnlen)1);
            }
            if (*nru > 0) {
                i__1 = m - ll + 1;
                dlasr_((char *)"R", (char *)"V", (char *)"F", nru, &i__1, &work[nm12 + 1], &work[nm13 + 1],
                       &u[ll * u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            }
            if (*ncc > 0) {
                i__1 = m - ll + 1;
                dlasr_((char *)"L", (char *)"V", (char *)"F", &i__1, ncc, &work[nm12 + 1], &work[nm13 + 1],
                       &c__[ll + c_dim1], ldc, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            }
            if ((d__1 = e[m - 1], abs(d__1)) <= thresh) {
                e[m - 1] = 0.;
            }
        } else {
            cs = 1.;
            oldcs = 1.;
            i__1 = ll + 1;
            for (i__ = m; i__ >= i__1; --i__) {
                d__1 = d__[i__] * cs;
                dlartg_(&d__1, &e[i__ - 1], &cs, &sn, &r__);
                if (i__ < m) {
                    e[i__] = oldsn * r__;
                }
                d__1 = oldcs * r__;
                d__2 = d__[i__ - 1] * sn;
                dlartg_(&d__1, &d__2, &oldcs, &oldsn, &d__[i__]);
                work[i__ - ll] = cs;
                work[i__ - ll + nm1] = -sn;
                work[i__ - ll + nm12] = oldcs;
                work[i__ - ll + nm13] = -oldsn;
            }
            h__ = d__[ll] * cs;
            d__[ll] = h__ * oldcs;
            e[ll] = h__ * oldsn;
            if (*ncvt > 0) {
                i__1 = m - ll + 1;
                dlasr_((char *)"L", (char *)"V", (char *)"B", &i__1, ncvt, &work[nm12 + 1], &work[nm13 + 1],
                       &vt[ll + vt_dim1], ldvt, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            }
            if (*nru > 0) {
                i__1 = m - ll + 1;
                dlasr_((char *)"R", (char *)"V", (char *)"B", nru, &i__1, &work[1], &work[*n], &u[ll * u_dim1 + 1], ldu,
                       (ftnlen)1, (ftnlen)1, (ftnlen)1);
            }
            if (*ncc > 0) {
                i__1 = m - ll + 1;
                dlasr_((char *)"L", (char *)"V", (char *)"B", &i__1, ncc, &work[1], &work[*n], &c__[ll + c_dim1], ldc,
                       (ftnlen)1, (ftnlen)1, (ftnlen)1);
            }
            if ((d__1 = e[ll], abs(d__1)) <= thresh) {
                e[ll] = 0.;
            }
        }
    } else {
        if (idir == 1) {
            f = ((d__1 = d__[ll], abs(d__1)) - shift) *
                (d_lmp_sign(&c_b49, &d__[ll]) + shift / d__[ll]);
            g = e[ll];
            i__1 = m - 1;
            for (i__ = ll; i__ <= i__1; ++i__) {
                dlartg_(&f, &g, &cosr, &sinr, &r__);
                if (i__ > ll) {
                    e[i__ - 1] = r__;
                }
                f = cosr * d__[i__] + sinr * e[i__];
                e[i__] = cosr * e[i__] - sinr * d__[i__];
                g = sinr * d__[i__ + 1];
                d__[i__ + 1] = cosr * d__[i__ + 1];
                dlartg_(&f, &g, &cosl, &sinl, &r__);
                d__[i__] = r__;
                f = cosl * e[i__] + sinl * d__[i__ + 1];
                d__[i__ + 1] = cosl * d__[i__ + 1] - sinl * e[i__];
                if (i__ < m - 1) {
                    g = sinl * e[i__ + 1];
                    e[i__ + 1] = cosl * e[i__ + 1];
                }
                work[i__ - ll + 1] = cosr;
                work[i__ - ll + 1 + nm1] = sinr;
                work[i__ - ll + 1 + nm12] = cosl;
                work[i__ - ll + 1 + nm13] = sinl;
            }
            e[m - 1] = f;
            if (*ncvt > 0) {
                i__1 = m - ll + 1;
                dlasr_((char *)"L", (char *)"V", (char *)"F", &i__1, ncvt, &work[1], &work[*n], &vt[ll + vt_dim1], ldvt,
                       (ftnlen)1, (ftnlen)1, (ftnlen)1);
            }
            if (*nru > 0) {
                i__1 = m - ll + 1;
                dlasr_((char *)"R", (char *)"V", (char *)"F", nru, &i__1, &work[nm12 + 1], &work[nm13 + 1],
                       &u[ll * u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            }
            if (*ncc > 0) {
                i__1 = m - ll + 1;
                dlasr_((char *)"L", (char *)"V", (char *)"F", &i__1, ncc, &work[nm12 + 1], &work[nm13 + 1],
                       &c__[ll + c_dim1], ldc, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            }
            if ((d__1 = e[m - 1], abs(d__1)) <= thresh) {
                e[m - 1] = 0.;
            }
        } else {
            f = ((d__1 = d__[m], abs(d__1)) - shift) * (d_lmp_sign(&c_b49, &d__[m]) + shift / d__[m]);
            g = e[m - 1];
            i__1 = ll + 1;
            for (i__ = m; i__ >= i__1; --i__) {
                dlartg_(&f, &g, &cosr, &sinr, &r__);
                if (i__ < m) {
                    e[i__] = r__;
                }
                f = cosr * d__[i__] + sinr * e[i__ - 1];
                e[i__ - 1] = cosr * e[i__ - 1] - sinr * d__[i__];
                g = sinr * d__[i__ - 1];
                d__[i__ - 1] = cosr * d__[i__ - 1];
                dlartg_(&f, &g, &cosl, &sinl, &r__);
                d__[i__] = r__;
                f = cosl * e[i__ - 1] + sinl * d__[i__ - 1];
                d__[i__ - 1] = cosl * d__[i__ - 1] - sinl * e[i__ - 1];
                if (i__ > ll + 1) {
                    g = sinl * e[i__ - 2];
                    e[i__ - 2] = cosl * e[i__ - 2];
                }
                work[i__ - ll] = cosr;
                work[i__ - ll + nm1] = -sinr;
                work[i__ - ll + nm12] = cosl;
                work[i__ - ll + nm13] = -sinl;
            }
            e[ll] = f;
            if ((d__1 = e[ll], abs(d__1)) <= thresh) {
                e[ll] = 0.;
            }
            if (*ncvt > 0) {
                i__1 = m - ll + 1;
                dlasr_((char *)"L", (char *)"V", (char *)"B", &i__1, ncvt, &work[nm12 + 1], &work[nm13 + 1],
                       &vt[ll + vt_dim1], ldvt, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            }
            if (*nru > 0) {
                i__1 = m - ll + 1;
                dlasr_((char *)"R", (char *)"V", (char *)"B", nru, &i__1, &work[1], &work[*n], &u[ll * u_dim1 + 1], ldu,
                       (ftnlen)1, (ftnlen)1, (ftnlen)1);
            }
            if (*ncc > 0) {
                i__1 = m - ll + 1;
                dlasr_((char *)"L", (char *)"V", (char *)"B", &i__1, ncc, &work[1], &work[*n], &c__[ll + c_dim1], ldc,
                       (ftnlen)1, (ftnlen)1, (ftnlen)1);
            }
        }
    }
    goto L60;
L160:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        if (d__[i__] < 0.) {
            d__[i__] = -d__[i__];
            if (*ncvt > 0) {
                dscal_(ncvt, &c_b72, &vt[i__ + vt_dim1], ldvt);
            }
        }
    }
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
        isub = 1;
        smin = d__[1];
        i__2 = *n + 1 - i__;
        for (j = 2; j <= i__2; ++j) {
            if (d__[j] <= smin) {
                isub = j;
                smin = d__[j];
            }
        }
        if (isub != *n + 1 - i__) {
            d__[isub] = d__[*n + 1 - i__];
            d__[*n + 1 - i__] = smin;
            if (*ncvt > 0) {
                dswap_(ncvt, &vt[isub + vt_dim1], ldvt, &vt[*n + 1 - i__ + vt_dim1], ldvt);
            }
            if (*nru > 0) {
                dswap_(nru, &u[isub * u_dim1 + 1], &c__1, &u[(*n + 1 - i__) * u_dim1 + 1], &c__1);
            }
            if (*ncc > 0) {
                dswap_(ncc, &c__[isub + c_dim1], ldc, &c__[*n + 1 - i__ + c_dim1], ldc);
            }
        }
    }
    goto L220;
L200:
    *info = 0;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
        if (e[i__] != 0.) {
            ++(*info);
        }
    }
L220:
    return 0;
}
#ifdef __cplusplus
}
#endif
