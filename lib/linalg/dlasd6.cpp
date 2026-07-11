#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__0 = 0;
static doublereal c_b7 = 1.;
static integer c__1 = 1;
static integer c_n1 = -1;
int dlasd6_(integer *icompq, integer *nl, integer *nr, integer *sqre, doublereal *d__,
            doublereal *vf, doublereal *vl, doublereal *alpha, doublereal *beta, integer *idxq,
            integer *perm, integer *givptr, integer *givcol, integer *ldgcol, doublereal *givnum,
            integer *ldgnum, doublereal *poles, doublereal *difl, doublereal *difr, doublereal *z__,
            integer *k, doublereal *c__, doublereal *s, doublereal *work, integer *iwork,
            integer *info)
{
    integer givcol_dim1, givcol_offset, givnum_dim1, givnum_offset, poles_dim1, poles_offset, i__1;
    doublereal d__1, d__2;
    integer i__, m, n, n1, n2, iw, idx, idxc, idxp, ivfw, ivlw;
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *),
        dlasd7_(integer *, integer *, integer *, integer *, integer *, doublereal *, doublereal *,
                doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                doublereal *, doublereal *, integer *, integer *, integer *, integer *, integer *,
                integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                integer *),
        dlasd8_(integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *,
                doublereal *, doublereal *, integer *, doublereal *, doublereal *, integer *),
        dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *,
                doublereal *, integer *, integer *, ftnlen),
        dlamrg_(integer *, integer *, doublereal *, integer *, integer *, integer *);
    integer isigma;
    extern int xerbla_(char *, integer *, ftnlen);
    doublereal orgnrm;
    --d__;
    --vf;
    --vl;
    --idxq;
    --perm;
    givcol_dim1 = *ldgcol;
    givcol_offset = 1 + givcol_dim1;
    givcol -= givcol_offset;
    poles_dim1 = *ldgnum;
    poles_offset = 1 + poles_dim1;
    poles -= poles_offset;
    givnum_dim1 = *ldgnum;
    givnum_offset = 1 + givnum_dim1;
    givnum -= givnum_offset;
    --difl;
    --difr;
    --z__;
    --work;
    --iwork;
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
        *info = -14;
    } else if (*ldgnum < n) {
        *info = -16;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLASD6", &i__1, (ftnlen)6);
        return 0;
    }
    isigma = 1;
    iw = isigma + n;
    ivfw = iw + m;
    ivlw = ivfw + m;
    idx = 1;
    idxc = idx + n;
    idxp = idxc + n;
    d__1 = abs(*alpha), d__2 = abs(*beta);
    orgnrm = max(d__1, d__2);
    d__[*nl + 1] = 0.;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        if ((d__1 = d__[i__], abs(d__1)) > orgnrm) {
            orgnrm = (d__1 = d__[i__], abs(d__1));
        }
    }
    dlascl_((char *)"G", &c__0, &c__0, &orgnrm, &c_b7, &n, &c__1, &d__[1], &n, info, (ftnlen)1);
    *alpha /= orgnrm;
    *beta /= orgnrm;
    dlasd7_(icompq, nl, nr, sqre, k, &d__[1], &z__[1], &work[iw], &vf[1], &work[ivfw], &vl[1],
            &work[ivlw], alpha, beta, &work[isigma], &iwork[idx], &iwork[idxp], &idxq[1], &perm[1],
            givptr, &givcol[givcol_offset], ldgcol, &givnum[givnum_offset], ldgnum, c__, s, info);
    dlasd8_(icompq, k, &d__[1], &z__[1], &vf[1], &vl[1], &difl[1], &difr[1], ldgnum, &work[isigma],
            &work[iw], info);
    if (*info != 0) {
        return 0;
    }
    if (*icompq == 1) {
        dcopy_(k, &d__[1], &c__1, &poles[poles_dim1 + 1], &c__1);
        dcopy_(k, &work[isigma], &c__1, &poles[(poles_dim1 << 1) + 1], &c__1);
    }
    dlascl_((char *)"G", &c__0, &c__0, &c_b7, &orgnrm, &n, &c__1, &d__[1], &n, info, (ftnlen)1);
    n1 = *k;
    n2 = n - *k;
    dlamrg_(&n1, &n2, &d__[1], &c__1, &c_n1, &idxq[1]);
    return 0;
}
#ifdef __cplusplus
}
#endif
