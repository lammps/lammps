#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static doublereal c_b30 = 0.;
int dlasd2_(integer *nl, integer *nr, integer *sqre, integer *k, doublereal *d__, doublereal *z__,
            doublereal *alpha, doublereal *beta, doublereal *u, integer *ldu, doublereal *vt,
            integer *ldvt, doublereal *dsigma, doublereal *u2, integer *ldu2, doublereal *vt2,
            integer *ldvt2, integer *idxp, integer *idx, integer *idxc, integer *idxq,
            integer *coltyp, integer *info)
{
    integer u_dim1, u_offset, u2_dim1, u2_offset, vt_dim1, vt_offset, vt2_dim1, vt2_offset, i__1;
    doublereal d__1, d__2;
    doublereal c__;
    integer i__, j, m, n;
    doublereal s;
    integer k2;
    doublereal z1;
    integer ct, jp;
    doublereal eps, tau, tol;
    integer psm[4], nlp1, nlp2, idxi, idxj;
    extern int drot_(integer *, doublereal *, integer *, doublereal *, integer *, doublereal *,
                     doublereal *);
    integer ctot[4], idxjp;
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer jprev;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, ftnlen);
    extern int dlamrg_(integer *, integer *, doublereal *, integer *, integer *, integer *),
        dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *,
                ftnlen),
        dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                ftnlen),
        xerbla_(char *, integer *, ftnlen);
    doublereal hlftol;
    --d__;
    --z__;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --dsigma;
    u2_dim1 = *ldu2;
    u2_offset = 1 + u2_dim1;
    u2 -= u2_offset;
    vt2_dim1 = *ldvt2;
    vt2_offset = 1 + vt2_dim1;
    vt2 -= vt2_offset;
    --idxp;
    --idx;
    --idxc;
    --idxq;
    --coltyp;
    *info = 0;
    if (*nl < 1) {
        *info = -1;
    } else if (*nr < 1) {
        *info = -2;
    } else if (*sqre != 1 && *sqre != 0) {
        *info = -3;
    }
    n = *nl + *nr + 1;
    m = n + *sqre;
    if (*ldu < n) {
        *info = -10;
    } else if (*ldvt < m) {
        *info = -12;
    } else if (*ldu2 < n) {
        *info = -15;
    } else if (*ldvt2 < m) {
        *info = -17;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLASD2", &i__1, (ftnlen)6);
        return 0;
    }
    nlp1 = *nl + 1;
    nlp2 = *nl + 2;
    z1 = *alpha * vt[nlp1 + nlp1 * vt_dim1];
    z__[1] = z1;
    for (i__ = *nl; i__ >= 1; --i__) {
        z__[i__ + 1] = *alpha * vt[i__ + nlp1 * vt_dim1];
        d__[i__ + 1] = d__[i__];
        idxq[i__ + 1] = idxq[i__] + 1;
    }
    i__1 = m;
    for (i__ = nlp2; i__ <= i__1; ++i__) {
        z__[i__] = *beta * vt[i__ + nlp2 * vt_dim1];
    }
    i__1 = nlp1;
    for (i__ = 2; i__ <= i__1; ++i__) {
        coltyp[i__] = 1;
    }
    i__1 = n;
    for (i__ = nlp2; i__ <= i__1; ++i__) {
        coltyp[i__] = 2;
    }
    i__1 = n;
    for (i__ = nlp2; i__ <= i__1; ++i__) {
        idxq[i__] += nlp1;
    }
    i__1 = n;
    for (i__ = 2; i__ <= i__1; ++i__) {
        dsigma[i__] = d__[idxq[i__]];
        u2[i__ + u2_dim1] = z__[idxq[i__]];
        idxc[i__] = coltyp[idxq[i__]];
    }
    dlamrg_(nl, nr, &dsigma[2], &c__1, &c__1, &idx[2]);
    i__1 = n;
    for (i__ = 2; i__ <= i__1; ++i__) {
        idxi = idx[i__] + 1;
        d__[i__] = dsigma[idxi];
        z__[i__] = u2[idxi + u2_dim1];
        coltyp[i__] = idxc[idxi];
    }
    eps = dlamch_((char *)"Epsilon", (ftnlen)7);
    d__1 = abs(*alpha), d__2 = abs(*beta);
    tol = max(d__1, d__2);
    d__2 = (d__1 = d__[n], abs(d__1));
    tol = eps * 8. * max(d__2, tol);
    *k = 1;
    k2 = n + 1;
    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
        if ((d__1 = z__[j], abs(d__1)) <= tol) {
            --k2;
            idxp[k2] = j;
            coltyp[j] = 4;
            if (j == n) {
                goto L120;
            }
        } else {
            jprev = j;
            goto L90;
        }
    }
L90:
    j = jprev;
L100:
    ++j;
    if (j > n) {
        goto L110;
    }
    if ((d__1 = z__[j], abs(d__1)) <= tol) {
        --k2;
        idxp[k2] = j;
        coltyp[j] = 4;
    } else {
        if ((d__1 = d__[j] - d__[jprev], abs(d__1)) <= tol) {
            s = z__[jprev];
            c__ = z__[j];
            tau = dlapy2_(&c__, &s);
            c__ /= tau;
            s = -s / tau;
            z__[j] = tau;
            z__[jprev] = 0.;
            idxjp = idxq[idx[jprev] + 1];
            idxj = idxq[idx[j] + 1];
            if (idxjp <= nlp1) {
                --idxjp;
            }
            if (idxj <= nlp1) {
                --idxj;
            }
            drot_(&n, &u[idxjp * u_dim1 + 1], &c__1, &u[idxj * u_dim1 + 1], &c__1, &c__, &s);
            drot_(&m, &vt[idxjp + vt_dim1], ldvt, &vt[idxj + vt_dim1], ldvt, &c__, &s);
            if (coltyp[j] != coltyp[jprev]) {
                coltyp[j] = 3;
            }
            coltyp[jprev] = 4;
            --k2;
            idxp[k2] = jprev;
            jprev = j;
        } else {
            ++(*k);
            u2[*k + u2_dim1] = z__[jprev];
            dsigma[*k] = d__[jprev];
            idxp[*k] = jprev;
            jprev = j;
        }
    }
    goto L100;
L110:
    ++(*k);
    u2[*k + u2_dim1] = z__[jprev];
    dsigma[*k] = d__[jprev];
    idxp[*k] = jprev;
L120:
    for (j = 1; j <= 4; ++j) {
        ctot[j - 1] = 0;
    }
    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
        ct = coltyp[j];
        ++ctot[ct - 1];
    }
    psm[0] = 2;
    psm[1] = ctot[0] + 2;
    psm[2] = psm[1] + ctot[1];
    psm[3] = psm[2] + ctot[2];
    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
        jp = idxp[j];
        ct = coltyp[jp];
        idxc[psm[ct - 1]] = j;
        ++psm[ct - 1];
    }
    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
        jp = idxp[j];
        dsigma[j] = d__[jp];
        idxj = idxq[idx[idxp[idxc[j]]] + 1];
        if (idxj <= nlp1) {
            --idxj;
        }
        dcopy_(&n, &u[idxj * u_dim1 + 1], &c__1, &u2[j * u2_dim1 + 1], &c__1);
        dcopy_(&m, &vt[idxj + vt_dim1], ldvt, &vt2[j + vt2_dim1], ldvt2);
    }
    dsigma[1] = 0.;
    hlftol = tol / 2.;
    if (abs(dsigma[2]) <= hlftol) {
        dsigma[2] = hlftol;
    }
    if (m > n) {
        z__[1] = dlapy2_(&z1, &z__[m]);
        if (z__[1] <= tol) {
            c__ = 1.;
            s = 0.;
            z__[1] = tol;
        } else {
            c__ = z1 / z__[1];
            s = z__[m] / z__[1];
        }
    } else {
        if (abs(z1) <= tol) {
            z__[1] = tol;
        } else {
            z__[1] = z1;
        }
    }
    i__1 = *k - 1;
    dcopy_(&i__1, &u2[u2_dim1 + 2], &c__1, &z__[2], &c__1);
    dlaset_((char *)"A", &n, &c__1, &c_b30, &c_b30, &u2[u2_offset], ldu2, (ftnlen)1);
    u2[nlp1 + u2_dim1] = 1.;
    if (m > n) {
        i__1 = nlp1;
        for (i__ = 1; i__ <= i__1; ++i__) {
            vt[m + i__ * vt_dim1] = -s * vt[nlp1 + i__ * vt_dim1];
            vt2[i__ * vt2_dim1 + 1] = c__ * vt[nlp1 + i__ * vt_dim1];
        }
        i__1 = m;
        for (i__ = nlp2; i__ <= i__1; ++i__) {
            vt2[i__ * vt2_dim1 + 1] = s * vt[m + i__ * vt_dim1];
            vt[m + i__ * vt_dim1] = c__ * vt[m + i__ * vt_dim1];
        }
    } else {
        dcopy_(&m, &vt[nlp1 + vt_dim1], ldvt, &vt2[vt2_dim1 + 1], ldvt2);
    }
    if (m > n) {
        dcopy_(&m, &vt[m + vt_dim1], ldvt, &vt2[m + vt2_dim1], ldvt2);
    }
    if (n > *k) {
        i__1 = n - *k;
        dcopy_(&i__1, &dsigma[*k + 1], &c__1, &d__[*k + 1], &c__1);
        i__1 = n - *k;
        dlacpy_((char *)"A", &n, &i__1, &u2[(*k + 1) * u2_dim1 + 1], ldu2, &u[(*k + 1) * u_dim1 + 1], ldu,
                (ftnlen)1);
        i__1 = n - *k;
        dlacpy_((char *)"A", &i__1, &m, &vt2[*k + 1 + vt2_dim1], ldvt2, &vt[*k + 1 + vt_dim1], ldvt,
                (ftnlen)1);
    }
    for (j = 1; j <= 4; ++j) {
        coltyp[j] = ctot[j - 1];
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
