#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__0 = 0;
static doublereal c_b7 = 1.;
static integer c__1 = 1;
static integer c_n1 = -1;
int dlasd1_(integer *nl, integer *nr, integer *sqre, doublereal *d__, doublereal *alpha,
            doublereal *beta, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt,
            integer *idxq, integer *iwork, doublereal *work, integer *info)
{
    integer u_dim1, u_offset, vt_dim1, vt_offset, i__1;
    doublereal d__1, d__2;
    integer i__, k, m, n, n1, n2, iq, iz, iu2, ldq, idx, ldu2, ivt2, idxc, idxp, ldvt2;
    extern int dlasd2_(integer *, integer *, integer *, integer *, doublereal *, doublereal *,
                       doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *,
                       doublereal *, doublereal *, integer *, doublereal *, integer *, integer *,
                       integer *, integer *, integer *, integer *, integer *),
        dlasd3_(integer *, integer *, integer *, integer *, doublereal *, doublereal *, integer *,
                doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *,
                integer *, doublereal *, integer *, integer *, integer *, doublereal *, integer *),
        dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *,
                doublereal *, integer *, integer *, ftnlen),
        dlamrg_(integer *, integer *, doublereal *, integer *, integer *, integer *);
    integer isigma;
    extern int xerbla_(char *, integer *, ftnlen);
    doublereal orgnrm;
    integer coltyp;
    --d__;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --idxq;
    --iwork;
    --work;
    *info = 0;
    if (*nl < 1) {
        *info = -1;
    } else if (*nr < 1) {
        *info = -2;
    } else if (*sqre < 0 || *sqre > 1) {
        *info = -3;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLASD1", &i__1, (ftnlen)6);
        return 0;
    }
    n = *nl + *nr + 1;
    m = n + *sqre;
    ldu2 = n;
    ldvt2 = m;
    iz = 1;
    isigma = iz + m;
    iu2 = isigma + n;
    ivt2 = iu2 + ldu2 * n;
    iq = ivt2 + ldvt2 * m;
    idx = 1;
    idxc = idx + n;
    coltyp = idxc + n;
    idxp = coltyp + n;
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
    dlasd2_(nl, nr, sqre, &k, &d__[1], &work[iz], alpha, beta, &u[u_offset], ldu, &vt[vt_offset],
            ldvt, &work[isigma], &work[iu2], &ldu2, &work[ivt2], &ldvt2, &iwork[idxp], &iwork[idx],
            &iwork[idxc], &idxq[1], &iwork[coltyp], info);
    ldq = k;
    dlasd3_(nl, nr, sqre, &k, &d__[1], &work[iq], &ldq, &work[isigma], &u[u_offset], ldu,
            &work[iu2], &ldu2, &vt[vt_offset], ldvt, &work[ivt2], &ldvt2, &iwork[idxc],
            &iwork[coltyp], &work[iz], info);
    if (*info != 0) {
        return 0;
    }
    dlascl_((char *)"G", &c__0, &c__0, &c_b7, &orgnrm, &n, &c__1, &d__[1], &n, info, (ftnlen)1);
    n1 = k;
    n2 = n - k;
    dlamrg_(&n1, &n2, &d__[1], &c__1, &c_n1, &idxq[1]);
    return 0;
}
#ifdef __cplusplus
}
#endif
