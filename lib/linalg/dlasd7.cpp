#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
int dlasd7_(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *k, doublereal *d__,
            doublereal *z__, doublereal *zw, doublereal *vf, doublereal *vfw, doublereal *vl,
            doublereal *vlw, doublereal *alpha, doublereal *beta, doublereal *dsigma, integer *idx,
            integer *idxp, integer *idxq, integer *perm, integer *givptr, integer *givcol,
            integer *ldgcol, doublereal *givnum, integer *ldgnum, doublereal *c__, doublereal *s,
            integer *info)
{
    integer givcol_dim1, givcol_offset, givnum_dim1, givnum_offset, i__1;
    doublereal d__1, d__2;
    integer i__, j, m, n, k2;
    doublereal z1;
    integer jp;
    doublereal eps, tau, tol;
    integer nlp1, nlp2, idxi, idxj;
    extern int drot_(integer *, doublereal *, integer *, doublereal *, integer *, doublereal *,
                     doublereal *);
    integer idxjp;
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer jprev;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, ftnlen);
    extern int dlamrg_(integer *, integer *, doublereal *, integer *, integer *, integer *),
        xerbla_(char *, integer *, ftnlen);
    doublereal hlftol;
    --d__;
    --z__;
    --zw;
    --vf;
    --vfw;
    --vl;
    --vlw;
    --dsigma;
    --idx;
    --idxp;
    --idxq;
    --perm;
    givcol_dim1 = *ldgcol;
    givcol_offset = 1 + givcol_dim1;
    givcol -= givcol_offset;
    givnum_dim1 = *ldgnum;
    givnum_offset = 1 + givnum_dim1;
    givnum -= givnum_offset;
    *info = 0;
    n = *nl + *nr + 1;
    m = n + *sqre;
    if (*icompq < 0 || *icompq > 1) {
        *info = -1;
    } else if (*nl < 1) {
        *info = -2;
    } else if (*nr < 1) {
        *info = -3;
    } else if (*sqre < 0 || *sqre > 1) {
        *info = -4;
    } else if (*ldgcol < n) {
        *info = -22;
    } else if (*ldgnum < n) {
        *info = -24;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLASD7", &i__1, (ftnlen)6);
        return 0;
    }
    nlp1 = *nl + 1;
    nlp2 = *nl + 2;
    if (*icompq == 1) {
        *givptr = 0;
    }
    z1 = *alpha * vl[nlp1];
    vl[nlp1] = 0.;
    tau = vf[nlp1];
    for (i__ = *nl; i__ >= 1; --i__) {
        z__[i__ + 1] = *alpha * vl[i__];
        vl[i__] = 0.;
        vf[i__ + 1] = vf[i__];
        d__[i__ + 1] = d__[i__];
        idxq[i__ + 1] = idxq[i__] + 1;
    }
    vf[1] = tau;
    i__1 = m;
    for (i__ = nlp2; i__ <= i__1; ++i__) {
        z__[i__] = *beta * vf[i__];
        vf[i__] = 0.;
    }
    i__1 = n;
    for (i__ = nlp2; i__ <= i__1; ++i__) {
        idxq[i__] += nlp1;
    }
    i__1 = n;
    for (i__ = 2; i__ <= i__1; ++i__) {
        dsigma[i__] = d__[idxq[i__]];
        zw[i__] = z__[idxq[i__]];
        vfw[i__] = vf[idxq[i__]];
        vlw[i__] = vl[idxq[i__]];
    }
    dlamrg_(nl, nr, &dsigma[2], &c__1, &c__1, &idx[2]);
    i__1 = n;
    for (i__ = 2; i__ <= i__1; ++i__) {
        idxi = idx[i__] + 1;
        d__[i__] = dsigma[idxi];
        z__[i__] = zw[idxi];
        vf[i__] = vfw[idxi];
        vl[i__] = vlw[idxi];
    }
    eps = dlamch_((char *)"Epsilon", (ftnlen)7);
    d__1 = abs(*alpha), d__2 = abs(*beta);
    tol = max(d__1, d__2);
    d__2 = (d__1 = d__[n], abs(d__1));
    tol = eps * 64. * max(d__2, tol);
    *k = 1;
    k2 = n + 1;
    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
        if ((d__1 = z__[j], abs(d__1)) <= tol) {
            --k2;
            idxp[k2] = j;
            if (j == n) {
                goto L100;
            }
        } else {
            jprev = j;
            goto L70;
        }
    }
L70:
    j = jprev;
L80:
    ++j;
    if (j > n) {
        goto L90;
    }
    if ((d__1 = z__[j], abs(d__1)) <= tol) {
        --k2;
        idxp[k2] = j;
    } else {
        if ((d__1 = d__[j] - d__[jprev], abs(d__1)) <= tol) {
            *s = z__[jprev];
            *c__ = z__[j];
            tau = dlapy2_(c__, s);
            z__[j] = tau;
            z__[jprev] = 0.;
            *c__ /= tau;
            *s = -(*s) / tau;
            if (*icompq == 1) {
                ++(*givptr);
                idxjp = idxq[idx[jprev] + 1];
                idxj = idxq[idx[j] + 1];
                if (idxjp <= nlp1) {
                    --idxjp;
                }
                if (idxj <= nlp1) {
                    --idxj;
                }
                givcol[*givptr + (givcol_dim1 << 1)] = idxjp;
                givcol[*givptr + givcol_dim1] = idxj;
                givnum[*givptr + (givnum_dim1 << 1)] = *c__;
                givnum[*givptr + givnum_dim1] = *s;
            }
            drot_(&c__1, &vf[jprev], &c__1, &vf[j], &c__1, c__, s);
            drot_(&c__1, &vl[jprev], &c__1, &vl[j], &c__1, c__, s);
            --k2;
            idxp[k2] = jprev;
            jprev = j;
        } else {
            ++(*k);
            zw[*k] = z__[jprev];
            dsigma[*k] = d__[jprev];
            idxp[*k] = jprev;
            jprev = j;
        }
    }
    goto L80;
L90:
    ++(*k);
    zw[*k] = z__[jprev];
    dsigma[*k] = d__[jprev];
    idxp[*k] = jprev;
L100:
    i__1 = n;
    for (j = 2; j <= i__1; ++j) {
        jp = idxp[j];
        dsigma[j] = d__[jp];
        vfw[j] = vf[jp];
        vlw[j] = vl[jp];
    }
    if (*icompq == 1) {
        i__1 = n;
        for (j = 2; j <= i__1; ++j) {
            jp = idxp[j];
            perm[j] = idxq[idx[jp] + 1];
            if (perm[j] <= nlp1) {
                --perm[j];
            }
        }
    }
    i__1 = n - *k;
    dcopy_(&i__1, &dsigma[*k + 1], &c__1, &d__[*k + 1], &c__1);
    dsigma[1] = 0.;
    hlftol = tol / 2.;
    if (abs(dsigma[2]) <= hlftol) {
        dsigma[2] = hlftol;
    }
    if (m > n) {
        z__[1] = dlapy2_(&z1, &z__[m]);
        if (z__[1] <= tol) {
            *c__ = 1.;
            *s = 0.;
            z__[1] = tol;
        } else {
            *c__ = z1 / z__[1];
            *s = -z__[m] / z__[1];
        }
        drot_(&c__1, &vf[m], &c__1, &vf[1], &c__1, c__, s);
        drot_(&c__1, &vl[m], &c__1, &vl[1], &c__1, c__, s);
    } else {
        if (abs(z1) <= tol) {
            z__[1] = tol;
        } else {
            z__[1] = z1;
        }
    }
    i__1 = *k - 1;
    dcopy_(&i__1, &zw[2], &c__1, &z__[2], &c__1);
    i__1 = n - 1;
    dcopy_(&i__1, &vfw[2], &c__1, &vf[2], &c__1);
    i__1 = n - 1;
    dcopy_(&i__1, &vlw[2], &c__1, &vl[2], &c__1);
    return 0;
}
#ifdef __cplusplus
}
#endif
