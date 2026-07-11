#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__2 = 2;
static integer c__1 = 1;
static doublereal c_b10 = 1.;
static doublereal c_b11 = 0.;
static integer c_n1 = -1;
int dlaed7_(integer *icompq, integer *n, integer *qsiz, integer *tlvls, integer *curlvl,
            integer *curpbm, doublereal *d__, doublereal *q, integer *ldq, integer *indxq,
            doublereal *rho, integer *cutpnt, doublereal *qstore, integer *qptr, integer *prmptr,
            integer *perm, integer *givptr, integer *givcol, doublereal *givnum, doublereal *work,
            integer *iwork, integer *info)
{
    integer q_dim1, q_offset, i__1, i__2;
    integer pow_lmp_ii(integer *, integer *);
    integer i__, k, n1, n2, is, iw, iz, iq2, ptr, ldq2, indx, curr;
    extern int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
                      integer *, doublereal *, integer *, doublereal *, doublereal *, integer *,
                      ftnlen, ftnlen);
    integer indxc, indxp;
    extern int dlaed8_(integer *, integer *, integer *, integer *, doublereal *, doublereal *,
                       integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                       doublereal *, integer *, doublereal *, integer *, integer *, integer *,
                       doublereal *, integer *, integer *, integer *),
        dlaed9_(integer *, integer *, integer *, integer *, doublereal *, doublereal *, integer *,
                doublereal *, doublereal *, doublereal *, doublereal *, integer *, integer *),
        dlaeda_(integer *, integer *, integer *, integer *, integer *, integer *, integer *,
                integer *, doublereal *, doublereal *, integer *, doublereal *, doublereal *,
                integer *);
    integer idlmda;
    extern int dlamrg_(integer *, integer *, doublereal *, integer *, integer *, integer *),
        xerbla_(char *, integer *, ftnlen);
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
    --iwork;
    *info = 0;
    if (*icompq < 0 || *icompq > 1) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*icompq == 1 && *qsiz < *n) {
        *info = -3;
    } else if (*ldq < max(1, *n)) {
        *info = -9;
    } else if (min(1, *n) > *cutpnt || *n < *cutpnt) {
        *info = -12;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLAED7", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    if (*icompq == 1) {
        ldq2 = *qsiz;
    } else {
        ldq2 = *n;
    }
    iz = 1;
    idlmda = iz + *n;
    iw = idlmda + *n;
    iq2 = iw + *n;
    is = iq2 + *n * ldq2;
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
            &qstore[1], &qptr[1], &work[iz], &work[iz + *n], info);
    if (*curlvl == *tlvls) {
        qptr[curr] = 1;
        prmptr[curr] = 1;
        givptr[curr] = 1;
    }
    dlaed8_(icompq, &k, n, qsiz, &d__[1], &q[q_offset], ldq, &indxq[1], rho, cutpnt, &work[iz],
            &work[idlmda], &work[iq2], &ldq2, &work[iw], &perm[prmptr[curr]], &givptr[curr + 1],
            &givcol[(givptr[curr] << 1) + 1], &givnum[(givptr[curr] << 1) + 1], &iwork[indxp],
            &iwork[indx], info);
    prmptr[curr + 1] = prmptr[curr] + *n;
    givptr[curr + 1] += givptr[curr];
    if (k != 0) {
        dlaed9_(&k, &c__1, &k, n, &d__[1], &work[is], &k, rho, &work[idlmda], &work[iw],
                &qstore[qptr[curr]], &k, info);
        if (*info != 0) {
            goto L30;
        }
        if (*icompq == 1) {
            dgemm_((char *)"N", (char *)"N", qsiz, &k, &k, &c_b10, &work[iq2], &ldq2, &qstore[qptr[curr]], &k,
                   &c_b11, &q[q_offset], ldq, (ftnlen)1, (ftnlen)1);
        }
        i__1 = k;
        qptr[curr + 1] = qptr[curr] + i__1 * i__1;
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
L30:
    return 0;
}
#ifdef __cplusplus
}
#endif
