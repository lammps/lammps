#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b9 = 0.;
static doublereal c_b10 = 1.;
static integer c__0 = 0;
static integer c__1 = 1;
static integer c__2 = 2;
int dsteqr_(char *compz, integer *n, doublereal *d__, doublereal *e, doublereal *z__, integer *ldz,
            doublereal *work, integer *info, ftnlen compz_len)
{
    integer z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;
    double sqrt(doublereal), d_lmp_sign(doublereal *, doublereal *);
    doublereal b, c__, f, g;
    integer i__, j, k, l, m;
    doublereal p, r__, s;
    integer l1, ii, mm, lm1, mm1, nm1;
    doublereal rt1, rt2, eps;
    integer lsv;
    doublereal tst, eps2;
    integer lend, jtot;
    extern int dlae2_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int dlasr_(char *, char *, char *, integer *, integer *, doublereal *, doublereal *,
                      doublereal *, integer *, ftnlen, ftnlen, ftnlen);
    doublereal anorm;
    extern int dswap_(integer *, doublereal *, integer *, doublereal *, integer *),
        dlaev2_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                doublereal *);
    integer lendm1, lendp1;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, ftnlen);
    integer iscale;
    extern int dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *,
                       integer *, doublereal *, integer *, integer *, ftnlen),
        dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                ftnlen);
    doublereal safmin;
    extern int dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
    doublereal safmax;
    extern int xerbla_(char *, integer *, ftnlen);
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *, ftnlen);
    extern int dlasrt_(char *, integer *, doublereal *, integer *, ftnlen);
    integer lendsv;
    doublereal ssfmin;
    integer nmaxit, icompz;
    doublereal ssfmax;
    --d__;
    --e;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    *info = 0;
    if (lsame_(compz, (char *)"N", (ftnlen)1, (ftnlen)1)) {
        icompz = 0;
    } else if (lsame_(compz, (char *)"V", (ftnlen)1, (ftnlen)1)) {
        icompz = 1;
    } else if (lsame_(compz, (char *)"I", (ftnlen)1, (ftnlen)1)) {
        icompz = 2;
    } else {
        icompz = -1;
    }
    if (icompz < 0) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*ldz < 1 || icompz > 0 && *ldz < max(1, *n)) {
        *info = -6;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DSTEQR", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    if (*n == 1) {
        if (icompz == 2) {
            z__[z_dim1 + 1] = 1.;
        }
        return 0;
    }
    eps = dlamch_((char *)"E", (ftnlen)1);
    d__1 = eps;
    eps2 = d__1 * d__1;
    safmin = dlamch_((char *)"S", (ftnlen)1);
    safmax = 1. / safmin;
    ssfmax = sqrt(safmax) / 3.;
    ssfmin = sqrt(safmin) / eps2;
    if (icompz == 2) {
        dlaset_((char *)"Full", n, n, &c_b9, &c_b10, &z__[z_offset], ldz, (ftnlen)4);
    }
    nmaxit = *n * 30;
    jtot = 0;
    l1 = 1;
    nm1 = *n - 1;
L10:
    if (l1 > *n) {
        goto L160;
    }
    if (l1 > 1) {
        e[l1 - 1] = 0.;
    }
    if (l1 <= nm1) {
        i__1 = nm1;
        for (m = l1; m <= i__1; ++m) {
            tst = (d__1 = e[m], abs(d__1));
            if (tst == 0.) {
                goto L30;
            }
            if (tst <=
                sqrt((d__1 = d__[m], abs(d__1))) * sqrt((d__2 = d__[m + 1], abs(d__2))) * eps) {
                e[m] = 0.;
                goto L30;
            }
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
    if ((d__1 = d__[lend], abs(d__1)) < (d__2 = d__[l], abs(d__2))) {
        lend = lsv;
        l = lendsv;
    }
    if (lend > l) {
    L40:
        if (l != lend) {
            lendm1 = lend - 1;
            i__1 = lendm1;
            for (m = l; m <= i__1; ++m) {
                d__2 = (d__1 = e[m], abs(d__1));
                tst = d__2 * d__2;
                if (tst <=
                    eps2 * (d__1 = d__[m], abs(d__1)) * (d__2 = d__[m + 1], abs(d__2)) + safmin) {
                    goto L60;
                }
            }
        }
        m = lend;
    L60:
        if (m < lend) {
            e[m] = 0.;
        }
        p = d__[l];
        if (m == l) {
            goto L80;
        }
        if (m == l + 1) {
            if (icompz > 0) {
                dlaev2_(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2, &c__, &s);
                work[l] = c__;
                work[*n - 1 + l] = s;
                dlasr_((char *)"R", (char *)"V", (char *)"B", n, &c__2, &work[l], &work[*n - 1 + l], &z__[l * z_dim1 + 1],
                       ldz, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            } else {
                dlae2_(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2);
            }
            d__[l] = rt1;
            d__[l + 1] = rt2;
            e[l] = 0.;
            l += 2;
            if (l <= lend) {
                goto L40;
            }
            goto L140;
        }
        if (jtot == nmaxit) {
            goto L140;
        }
        ++jtot;
        g = (d__[l + 1] - p) / (e[l] * 2.);
        r__ = dlapy2_(&g, &c_b10);
        g = d__[m] - p + e[l] / (g + d_lmp_sign(&r__, &g));
        s = 1.;
        c__ = 1.;
        p = 0.;
        mm1 = m - 1;
        i__1 = l;
        for (i__ = mm1; i__ >= i__1; --i__) {
            f = s * e[i__];
            b = c__ * e[i__];
            dlartg_(&g, &f, &c__, &s, &r__);
            if (i__ != m - 1) {
                e[i__ + 1] = r__;
            }
            g = d__[i__ + 1] - p;
            r__ = (d__[i__] - g) * s + c__ * 2. * b;
            p = s * r__;
            d__[i__ + 1] = g + p;
            g = c__ * r__ - b;
            if (icompz > 0) {
                work[i__] = c__;
                work[*n - 1 + i__] = -s;
            }
        }
        if (icompz > 0) {
            mm = m - l + 1;
            dlasr_((char *)"R", (char *)"V", (char *)"B", n, &mm, &work[l], &work[*n - 1 + l], &z__[l * z_dim1 + 1], ldz,
                   (ftnlen)1, (ftnlen)1, (ftnlen)1);
        }
        d__[l] -= p;
        e[l] = g;
        goto L40;
    L80:
        d__[l] = p;
        ++l;
        if (l <= lend) {
            goto L40;
        }
        goto L140;
    } else {
    L90:
        if (l != lend) {
            lendp1 = lend + 1;
            i__1 = lendp1;
            for (m = l; m >= i__1; --m) {
                d__2 = (d__1 = e[m - 1], abs(d__1));
                tst = d__2 * d__2;
                if (tst <=
                    eps2 * (d__1 = d__[m], abs(d__1)) * (d__2 = d__[m - 1], abs(d__2)) + safmin) {
                    goto L110;
                }
            }
        }
        m = lend;
    L110:
        if (m > lend) {
            e[m - 1] = 0.;
        }
        p = d__[l];
        if (m == l) {
            goto L130;
        }
        if (m == l - 1) {
            if (icompz > 0) {
                dlaev2_(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2, &c__, &s);
                work[m] = c__;
                work[*n - 1 + m] = s;
                dlasr_((char *)"R", (char *)"V", (char *)"F", n, &c__2, &work[m], &work[*n - 1 + m],
                       &z__[(l - 1) * z_dim1 + 1], ldz, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            } else {
                dlae2_(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2);
            }
            d__[l - 1] = rt1;
            d__[l] = rt2;
            e[l - 1] = 0.;
            l += -2;
            if (l >= lend) {
                goto L90;
            }
            goto L140;
        }
        if (jtot == nmaxit) {
            goto L140;
        }
        ++jtot;
        g = (d__[l - 1] - p) / (e[l - 1] * 2.);
        r__ = dlapy2_(&g, &c_b10);
        g = d__[m] - p + e[l - 1] / (g + d_lmp_sign(&r__, &g));
        s = 1.;
        c__ = 1.;
        p = 0.;
        lm1 = l - 1;
        i__1 = lm1;
        for (i__ = m; i__ <= i__1; ++i__) {
            f = s * e[i__];
            b = c__ * e[i__];
            dlartg_(&g, &f, &c__, &s, &r__);
            if (i__ != m) {
                e[i__ - 1] = r__;
            }
            g = d__[i__] - p;
            r__ = (d__[i__ + 1] - g) * s + c__ * 2. * b;
            p = s * r__;
            d__[i__] = g + p;
            g = c__ * r__ - b;
            if (icompz > 0) {
                work[i__] = c__;
                work[*n - 1 + i__] = s;
            }
        }
        if (icompz > 0) {
            mm = l - m + 1;
            dlasr_((char *)"R", (char *)"V", (char *)"F", n, &mm, &work[m], &work[*n - 1 + m], &z__[m * z_dim1 + 1], ldz,
                   (ftnlen)1, (ftnlen)1, (ftnlen)1);
        }
        d__[l] -= p;
        e[lm1] = g;
        goto L90;
    L130:
        d__[l] = p;
        --l;
        if (l >= lend) {
            goto L90;
        }
        goto L140;
    }
L140:
    if (iscale == 1) {
        i__1 = lendsv - lsv + 1;
        dlascl_((char *)"G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d__[lsv], n, info, (ftnlen)1);
        i__1 = lendsv - lsv;
        dlascl_((char *)"G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &e[lsv], n, info, (ftnlen)1);
    } else if (iscale == 2) {
        i__1 = lendsv - lsv + 1;
        dlascl_((char *)"G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d__[lsv], n, info, (ftnlen)1);
        i__1 = lendsv - lsv;
        dlascl_((char *)"G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &e[lsv], n, info, (ftnlen)1);
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
    goto L190;
L160:
    if (icompz == 0) {
        dlasrt_((char *)"I", n, &d__[1], info, (ftnlen)1);
    } else {
        i__1 = *n;
        for (ii = 2; ii <= i__1; ++ii) {
            i__ = ii - 1;
            k = i__;
            p = d__[i__];
            i__2 = *n;
            for (j = ii; j <= i__2; ++j) {
                if (d__[j] < p) {
                    k = j;
                    p = d__[j];
                }
            }
            if (k != i__) {
                d__[k] = d__[i__];
                d__[i__] = p;
                dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[k * z_dim1 + 1], &c__1);
            }
        }
    }
L190:
    return 0;
}
#ifdef __cplusplus
}
#endif
