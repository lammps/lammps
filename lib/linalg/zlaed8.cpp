#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b3 = -1.;
static integer c__1 = 1;
int zlaed8_(integer *k, integer *n, integer *qsiz, doublecomplex *q, integer *ldq, doublereal *d__,
            doublereal *rho, integer *cutpnt, doublereal *z__, doublereal *dlamda,
            doublecomplex *q2, integer *ldq2, doublereal *w, integer *indxp, integer *indx,
            integer *indxq, integer *perm, integer *givptr, integer *givcol, doublereal *givnum,
            integer *info)
{
    integer q_dim1, q_offset, q2_dim1, q2_offset, i__1;
    doublereal d__1;
    double sqrt(doublereal);
    doublereal c__;
    integer i__, j;
    doublereal s, t;
    integer k2, n1, n2, jp, n1p1;
    doublereal eps, tau, tol;
    integer jlam, imax, jmax;
    extern int dscal_(integer *, doublereal *, doublereal *, integer *),
        dcopy_(integer *, doublereal *, integer *, doublereal *, integer *),
        zdrot_(integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublereal *,
               doublereal *),
        zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, integer *);
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern int dlamrg_(integer *, integer *, doublereal *, integer *, integer *, integer *),
        xerbla_(char *, integer *, ftnlen),
        zlacpy_(char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *,
                integer *, ftnlen);
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --d__;
    --z__;
    --dlamda;
    q2_dim1 = *ldq2;
    q2_offset = 1 + q2_dim1;
    q2 -= q2_offset;
    --w;
    --indxp;
    --indx;
    --indxq;
    --perm;
    givcol -= 3;
    givnum -= 3;
    *info = 0;
    if (*n < 0) {
        *info = -2;
    } else if (*qsiz < *n) {
        *info = -3;
    } else if (*ldq < max(1, *n)) {
        *info = -5;
    } else if (*cutpnt < min(1, *n) || *cutpnt > *n) {
        *info = -8;
    } else if (*ldq2 < max(1, *n)) {
        *info = -12;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"ZLAED8", &i__1, (ftnlen)6);
        return 0;
    }
    *givptr = 0;
    if (*n == 0) {
        return 0;
    }
    n1 = *cutpnt;
    n2 = *n - n1;
    n1p1 = n1 + 1;
    if (*rho < 0.) {
        dscal_(&n2, &c_b3, &z__[n1p1], &c__1);
    }
    t = 1. / sqrt(2.);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
        indx[j] = j;
    }
    dscal_(n, &t, &z__[1], &c__1);
    *rho = (d__1 = *rho * 2., abs(d__1));
    i__1 = *n;
    for (i__ = *cutpnt + 1; i__ <= i__1; ++i__) {
        indxq[i__] += *cutpnt;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        dlamda[i__] = d__[indxq[i__]];
        w[i__] = z__[indxq[i__]];
    }
    i__ = 1;
    j = *cutpnt + 1;
    dlamrg_(&n1, &n2, &dlamda[1], &c__1, &c__1, &indx[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        d__[i__] = dlamda[indx[i__]];
        z__[i__] = w[indx[i__]];
    }
    imax = idamax_(n, &z__[1], &c__1);
    jmax = idamax_(n, &d__[1], &c__1);
    eps = dlamch_((char *)"Epsilon", (ftnlen)7);
    tol = eps * 8. * (d__1 = d__[jmax], abs(d__1));
    if (*rho * (d__1 = z__[imax], abs(d__1)) <= tol) {
        *k = 0;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            perm[j] = indxq[indx[j]];
            zcopy_(qsiz, &q[perm[j] * q_dim1 + 1], &c__1, &q2[j * q2_dim1 + 1], &c__1);
        }
        zlacpy_((char *)"A", qsiz, n, &q2[q2_dim1 + 1], ldq2, &q[q_dim1 + 1], ldq, (ftnlen)1);
        return 0;
    }
    *k = 0;
    k2 = *n + 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
        if (*rho * (d__1 = z__[j], abs(d__1)) <= tol) {
            --k2;
            indxp[k2] = j;
            if (j == *n) {
                goto L100;
            }
        } else {
            jlam = j;
            goto L70;
        }
    }
L70:
    ++j;
    if (j > *n) {
        goto L90;
    }
    if (*rho * (d__1 = z__[j], abs(d__1)) <= tol) {
        --k2;
        indxp[k2] = j;
    } else {
        s = z__[jlam];
        c__ = z__[j];
        tau = dlapy2_(&c__, &s);
        t = d__[j] - d__[jlam];
        c__ /= tau;
        s = -s / tau;
        if ((d__1 = t * c__ * s, abs(d__1)) <= tol) {
            z__[j] = tau;
            z__[jlam] = 0.;
            ++(*givptr);
            givcol[(*givptr << 1) + 1] = indxq[indx[jlam]];
            givcol[(*givptr << 1) + 2] = indxq[indx[j]];
            givnum[(*givptr << 1) + 1] = c__;
            givnum[(*givptr << 1) + 2] = s;
            zdrot_(qsiz, &q[indxq[indx[jlam]] * q_dim1 + 1], &c__1, &q[indxq[indx[j]] * q_dim1 + 1],
                   &c__1, &c__, &s);
            t = d__[jlam] * c__ * c__ + d__[j] * s * s;
            d__[j] = d__[jlam] * s * s + d__[j] * c__ * c__;
            d__[jlam] = t;
            --k2;
            i__ = 1;
        L80:
            if (k2 + i__ <= *n) {
                if (d__[jlam] < d__[indxp[k2 + i__]]) {
                    indxp[k2 + i__ - 1] = indxp[k2 + i__];
                    indxp[k2 + i__] = jlam;
                    ++i__;
                    goto L80;
                } else {
                    indxp[k2 + i__ - 1] = jlam;
                }
            } else {
                indxp[k2 + i__ - 1] = jlam;
            }
            jlam = j;
        } else {
            ++(*k);
            w[*k] = z__[jlam];
            dlamda[*k] = d__[jlam];
            indxp[*k] = jlam;
            jlam = j;
        }
    }
    goto L70;
L90:
    ++(*k);
    w[*k] = z__[jlam];
    dlamda[*k] = d__[jlam];
    indxp[*k] = jlam;
L100:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
        jp = indxp[j];
        dlamda[j] = d__[jp];
        perm[j] = indxq[indx[jp]];
        zcopy_(qsiz, &q[perm[j] * q_dim1 + 1], &c__1, &q2[j * q2_dim1 + 1], &c__1);
    }
    if (*k < *n) {
        i__1 = *n - *k;
        dcopy_(&i__1, &dlamda[*k + 1], &c__1, &d__[*k + 1], &c__1);
        i__1 = *n - *k;
        zlacpy_((char *)"A", qsiz, &i__1, &q2[(*k + 1) * q2_dim1 + 1], ldq2, &q[(*k + 1) * q_dim1 + 1], ldq,
                (ftnlen)1);
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
