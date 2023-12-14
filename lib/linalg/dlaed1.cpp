#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
int dlaed1_(integer *n, doublereal *d__, doublereal *q, integer *ldq, integer *indxq,
            doublereal *rho, integer *cutpnt, doublereal *work, integer *iwork, integer *info)
{
    integer q_dim1, q_offset, i__1, i__2;
    integer i__, k, n1, n2, is, iw, iz, iq2, zpp1, indx, indxc;
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer indxp;
    extern int dlaed2_(integer *, integer *, integer *, doublereal *, doublereal *, integer *,
                       integer *, doublereal *, doublereal *, doublereal *, doublereal *,
                       doublereal *, integer *, integer *, integer *, integer *, integer *),
        dlaed3_(integer *, integer *, integer *, doublereal *, doublereal *, integer *,
                doublereal *, doublereal *, doublereal *, integer *, integer *, doublereal *,
                doublereal *, integer *);
    integer idlmda;
    extern int dlamrg_(integer *, integer *, doublereal *, integer *, integer *, integer *),
        xerbla_(char *, integer *, ftnlen);
    integer coltyp;
    --d__;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --indxq;
    --work;
    --iwork;
    *info = 0;
    if (*n < 0) {
        *info = -1;
    } else if (*ldq < max(1, *n)) {
        *info = -4;
    } else {
        i__1 = 1, i__2 = *n / 2;
        if (min(i__1, i__2) > *cutpnt || *n / 2 < *cutpnt) {
            *info = -7;
        }
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLAED1", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    iz = 1;
    idlmda = iz + *n;
    iw = idlmda + *n;
    iq2 = iw + *n;
    indx = 1;
    indxc = indx + *n;
    coltyp = indxc + *n;
    indxp = coltyp + *n;
    dcopy_(cutpnt, &q[*cutpnt + q_dim1], ldq, &work[iz], &c__1);
    zpp1 = *cutpnt + 1;
    i__1 = *n - *cutpnt;
    dcopy_(&i__1, &q[zpp1 + zpp1 * q_dim1], ldq, &work[iz + *cutpnt], &c__1);
    dlaed2_(&k, n, cutpnt, &d__[1], &q[q_offset], ldq, &indxq[1], rho, &work[iz], &work[idlmda],
            &work[iw], &work[iq2], &iwork[indx], &iwork[indxc], &iwork[indxp], &iwork[coltyp],
            info);
    if (*info != 0) {
        goto L20;
    }
    if (k != 0) {
        is = (iwork[coltyp] + iwork[coltyp + 1]) * *cutpnt +
             (iwork[coltyp + 1] + iwork[coltyp + 2]) * (*n - *cutpnt) + iq2;
        dlaed3_(&k, n, cutpnt, &d__[1], &q[q_offset], ldq, rho, &work[idlmda], &work[iq2],
                &iwork[indxc], &iwork[coltyp], &work[iw], &work[is], info);
        if (*info != 0) {
            goto L20;
        }
        n1 = k;
        n2 = *n - k;
        dlamrg_(&n1, &n2, &d__[1], &c__1, &c_n1, &indxq[1]);
    } else {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            indxq[i__] = i__;
        }
    }
L20:
    return 0;
}
#ifdef __cplusplus
}
#endif
