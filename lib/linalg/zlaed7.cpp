#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__2 = 2;
static integer c__1 = 1;
static integer c_n1 = -1;
int zlaed7_(integer *n, integer *cutpnt, integer *qsiz, integer *tlvls, integer *curlvl,
            integer *curpbm, doublereal *d__, doublecomplex *q, integer *ldq, doublereal *rho,
            integer *indxq, doublereal *qstore, integer *qptr, integer *prmptr, integer *perm,
            integer *givptr, integer *givcol, doublereal *givnum, doublecomplex *work,
            doublereal *rwork, integer *iwork, integer *info)
{
    integer q_dim1, q_offset, i__1, i__2;
    integer pow_lmp_ii(integer *, integer *);
    integer i__, k, n1, n2, iq, iw, iz, ptr, indx, curr, indxc, indxp;
    extern int dlaed9_(integer *, integer *, integer *, integer *, doublereal *, doublereal *,
                       integer *, doublereal *, doublereal *, doublereal *, doublereal *, integer *,
                       integer *),
        zlaed8_(integer *, integer *, integer *, doublecomplex *, integer *, doublereal *,
                doublereal *, integer *, doublereal *, doublereal *, doublecomplex *, integer *,
                doublereal *, integer *, integer *, integer *, integer *, integer *, integer *,
                doublereal *, integer *),
        dlaeda_(integer *, integer *, integer *, integer *, integer *, integer *, integer *,
                integer *, doublereal *, doublereal *, integer *, doublereal *, doublereal *,
                integer *);
    integer idlmda;
    extern int dlamrg_(integer *, integer *, doublereal *, integer *, integer *, integer *),
        xerbla_(char *, integer *, ftnlen),
        zlacrm_(integer *, integer *, doublecomplex *, integer *, doublereal *, integer *,
                doublecomplex *, integer *, doublereal *);
    integer coltyp;
    --d__;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --indxq;
    --qstore;
    --qptr;
    --prmptr;
    --perm;
    --givptr;
    givcol -= 3;
    givnum -= 3;
    --work;
    --rwork;
    --iwork;
    *info = 0;
    if (*n < 0) {
        *info = -1;
    } else if (min(1, *n) > *cutpnt || *n < *cutpnt) {
        *info = -2;
    } else if (*qsiz < *n) {
        *info = -3;
    } else if (*ldq < max(1, *n)) {
        *info = -9;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"ZLAED7", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    iz = 1;
    idlmda = iz + *n;
    iw = idlmda + *n;
    iq = iw + *n;
    indx = 1;
    indxc = indx + *n;
    coltyp = indxc + *n;
    indxp = coltyp + *n;
    ptr = pow_lmp_ii(&c__2, tlvls) + 1;
    i__1 = *curlvl - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
        i__2 = *tlvls - i__;
        ptr += pow_lmp_ii(&c__2, &i__2);
    }
    curr = ptr + *curpbm;
    dlaeda_(n, tlvls, curlvl, curpbm, &prmptr[1], &perm[1], &givptr[1], &givcol[3], &givnum[3],
            &qstore[1], &qptr[1], &rwork[iz], &rwork[iz + *n], info);
    if (*curlvl == *tlvls) {
        qptr[curr] = 1;
        prmptr[curr] = 1;
        givptr[curr] = 1;
    }
    zlaed8_(&k, n, qsiz, &q[q_offset], ldq, &d__[1], rho, cutpnt, &rwork[iz], &rwork[idlmda],
            &work[1], qsiz, &rwork[iw], &iwork[indxp], &iwork[indx], &indxq[1], &perm[prmptr[curr]],
            &givptr[curr + 1], &givcol[(givptr[curr] << 1) + 1], &givnum[(givptr[curr] << 1) + 1],
            info);
    prmptr[curr + 1] = prmptr[curr] + *n;
    givptr[curr + 1] += givptr[curr];
    if (k != 0) {
        dlaed9_(&k, &c__1, &k, n, &d__[1], &rwork[iq], &k, rho, &rwork[idlmda], &rwork[iw],
                &qstore[qptr[curr]], &k, info);
        zlacrm_(qsiz, &k, &work[1], qsiz, &qstore[qptr[curr]], &k, &q[q_offset], ldq, &rwork[iq]);
        i__1 = k;
        qptr[curr + 1] = qptr[curr] + i__1 * i__1;
        if (*info != 0) {
            return 0;
        }
        n1 = k;
        n2 = *n - k;
        dlamrg_(&n1, &n2, &d__[1], &c__1, &c_n1, &indxq[1]);
    } else {
        qptr[curr + 1] = qptr[curr];
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            indxq[i__] = i__;
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
