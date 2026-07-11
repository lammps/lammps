#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b3 = -1.;
static integer c__1 = 1;
int dlaed2_(integer *k, integer *n, integer *n1, doublereal *d__, doublereal *q, integer *ldq,
            integer *indxq, doublereal *rho, doublereal *z__, doublereal *dlambda, doublereal *w,
            doublereal *q2, integer *indx, integer *indxc, integer *indxp, integer *coltyp,
            integer *info)
{
    integer q_dim1, q_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;
    double sqrt(doublereal);
    doublereal c__;
    integer i__, j;
    doublereal s, t;
    integer k2, n2, ct, nj, pj, js, iq1, iq2, n1p1;
    doublereal eps, tau, tol;
    integer psm[4], imax, jmax;
    extern int drot_(integer *, doublereal *, integer *, doublereal *, integer *, doublereal *,
                     doublereal *);
    integer ctot[4];
    extern int dscal_(integer *, doublereal *, doublereal *, integer *),
        dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern int dlamrg_(integer *, integer *, doublereal *, integer *, integer *, integer *),
        dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *,
                ftnlen),
        xerbla_(char *, integer *, ftnlen);
    --d__;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --indxq;
    --z__;
    --dlambda;
    --w;
    --q2;
    --indx;
    --indxc;
    --indxp;
    --coltyp;
    *info = 0;
    if (*n < 0) {
        *info = -2;
    } else if (*ldq < max(1, *n)) {
        *info = -6;
    } else {
        i__1 = 1, i__2 = *n / 2;
        if (min(i__1, i__2) > *n1 || *n / 2 < *n1) {
            *info = -3;
        }
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLAED2", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    n2 = *n - *n1;
    n1p1 = *n1 + 1;
    if (*rho < 0.) {
        dscal_(&n2, &c_b3, &z__[n1p1], &c__1);
    }
    t = 1. / sqrt(2.);
    dscal_(n, &t, &z__[1], &c__1);
    *rho = (d__1 = *rho * 2., abs(d__1));
    i__1 = *n;
    for (i__ = n1p1; i__ <= i__1; ++i__) {
        indxq[i__] += *n1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        dlambda[i__] = d__[indxq[i__]];
    }
    dlamrg_(n1, &n2, &dlambda[1], &c__1, &c__1, &indxc[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        indx[i__] = indxq[indxc[i__]];
    }
    imax = idamax_(n, &z__[1], &c__1);
    jmax = idamax_(n, &d__[1], &c__1);
    eps = dlamch_((char *)"Epsilon", (ftnlen)7);
    d__3 = (d__1 = d__[jmax], abs(d__1)), d__4 = (d__2 = z__[imax], abs(d__2));
    tol = eps * 8. * max(d__3, d__4);
    if (*rho * (d__1 = z__[imax], abs(d__1)) <= tol) {
        *k = 0;
        iq2 = 1;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__ = indx[j];
            dcopy_(n, &q[i__ * q_dim1 + 1], &c__1, &q2[iq2], &c__1);
            dlambda[j] = d__[i__];
            iq2 += *n;
        }
        dlacpy_((char *)"A", n, n, &q2[1], n, &q[q_offset], ldq, (ftnlen)1);
        dcopy_(n, &dlambda[1], &c__1, &d__[1], &c__1);
        goto L190;
    }
    i__1 = *n1;
    for (i__ = 1; i__ <= i__1; ++i__) {
        coltyp[i__] = 1;
    }
    i__1 = *n;
    for (i__ = n1p1; i__ <= i__1; ++i__) {
        coltyp[i__] = 3;
    }
    *k = 0;
    k2 = *n + 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
        nj = indx[j];
        if (*rho * (d__1 = z__[nj], abs(d__1)) <= tol) {
            --k2;
            coltyp[nj] = 4;
            indxp[k2] = nj;
            if (j == *n) {
                goto L100;
            }
        } else {
            pj = nj;
            goto L80;
        }
    }
L80:
    ++j;
    nj = indx[j];
    if (j > *n) {
        goto L100;
    }
    if (*rho * (d__1 = z__[nj], abs(d__1)) <= tol) {
        --k2;
        coltyp[nj] = 4;
        indxp[k2] = nj;
    } else {
        s = z__[pj];
        c__ = z__[nj];
        tau = dlapy2_(&c__, &s);
        t = d__[nj] - d__[pj];
        c__ /= tau;
        s = -s / tau;
        if ((d__1 = t * c__ * s, abs(d__1)) <= tol) {
            z__[nj] = tau;
            z__[pj] = 0.;
            if (coltyp[nj] != coltyp[pj]) {
                coltyp[nj] = 2;
            }
            coltyp[pj] = 4;
            drot_(n, &q[pj * q_dim1 + 1], &c__1, &q[nj * q_dim1 + 1], &c__1, &c__, &s);
            d__1 = c__;
            d__2 = s;
            t = d__[pj] * (d__1 * d__1) + d__[nj] * (d__2 * d__2);
            d__1 = s;
            d__2 = c__;
            d__[nj] = d__[pj] * (d__1 * d__1) + d__[nj] * (d__2 * d__2);
            d__[pj] = t;
            --k2;
            i__ = 1;
        L90:
            if (k2 + i__ <= *n) {
                if (d__[pj] < d__[indxp[k2 + i__]]) {
                    indxp[k2 + i__ - 1] = indxp[k2 + i__];
                    indxp[k2 + i__] = pj;
                    ++i__;
                    goto L90;
                } else {
                    indxp[k2 + i__ - 1] = pj;
                }
            } else {
                indxp[k2 + i__ - 1] = pj;
            }
            pj = nj;
        } else {
            ++(*k);
            dlambda[*k] = d__[pj];
            w[*k] = z__[pj];
            indxp[*k] = pj;
            pj = nj;
        }
    }
    goto L80;
L100:
    ++(*k);
    dlambda[*k] = d__[pj];
    w[*k] = z__[pj];
    indxp[*k] = pj;
    for (j = 1; j <= 4; ++j) {
        ctot[j - 1] = 0;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
        ct = coltyp[j];
        ++ctot[ct - 1];
    }
    psm[0] = 1;
    psm[1] = ctot[0] + 1;
    psm[2] = psm[1] + ctot[1];
    psm[3] = psm[2] + ctot[2];
    *k = *n - ctot[3];
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
        js = indxp[j];
        ct = coltyp[js];
        indx[psm[ct - 1]] = js;
        indxc[psm[ct - 1]] = j;
        ++psm[ct - 1];
    }
    i__ = 1;
    iq1 = 1;
    iq2 = (ctot[0] + ctot[1]) * *n1 + 1;
    i__1 = ctot[0];
    for (j = 1; j <= i__1; ++j) {
        js = indx[i__];
        dcopy_(n1, &q[js * q_dim1 + 1], &c__1, &q2[iq1], &c__1);
        z__[i__] = d__[js];
        ++i__;
        iq1 += *n1;
    }
    i__1 = ctot[1];
    for (j = 1; j <= i__1; ++j) {
        js = indx[i__];
        dcopy_(n1, &q[js * q_dim1 + 1], &c__1, &q2[iq1], &c__1);
        dcopy_(&n2, &q[*n1 + 1 + js * q_dim1], &c__1, &q2[iq2], &c__1);
        z__[i__] = d__[js];
        ++i__;
        iq1 += *n1;
        iq2 += n2;
    }
    i__1 = ctot[2];
    for (j = 1; j <= i__1; ++j) {
        js = indx[i__];
        dcopy_(&n2, &q[*n1 + 1 + js * q_dim1], &c__1, &q2[iq2], &c__1);
        z__[i__] = d__[js];
        ++i__;
        iq2 += n2;
    }
    iq1 = iq2;
    i__1 = ctot[3];
    for (j = 1; j <= i__1; ++j) {
        js = indx[i__];
        dcopy_(n, &q[js * q_dim1 + 1], &c__1, &q2[iq2], &c__1);
        iq2 += *n;
        z__[i__] = d__[js];
        ++i__;
    }
    if (*k < *n) {
        dlacpy_((char *)"A", n, &ctot[3], &q2[iq1], n, &q[(*k + 1) * q_dim1 + 1], ldq, (ftnlen)1);
        i__1 = *n - *k;
        dcopy_(&i__1, &z__[*k + 1], &c__1, &d__[*k + 1], &c__1);
    }
    for (j = 1; j <= 4; ++j) {
        coltyp[j] = ctot[j - 1];
    }
L190:
    return 0;
}
#ifdef __cplusplus
}
#endif
