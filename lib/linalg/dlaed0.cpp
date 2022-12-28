#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__9 = 9;
static integer c__0 = 0;
static integer c__2 = 2;
static doublereal c_b23 = 1.;
static doublereal c_b24 = 0.;
static integer c__1 = 1;
int dlaed0_(integer *icompq, integer *qsiz, integer *n, doublereal *d__, doublereal *e,
            doublereal *q, integer *ldq, doublereal *qstore, integer *ldqs, doublereal *work,
            integer *iwork, integer *info)
{
    integer q_dim1, q_offset, qstore_dim1, qstore_offset, i__1, i__2;
    doublereal d__1;
    double log(doublereal);
    integer pow_lmp_ii(integer *, integer *);
    integer i__, j, k, iq, lgn, msd2, smm1, spm1, spm2;
    doublereal temp;
    integer curr;
    extern int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
                      integer *, doublereal *, integer *, doublereal *, doublereal *, integer *,
                      ftnlen, ftnlen);
    integer iperm;
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer indxq, iwrem;
    extern int dlaed1_(integer *, doublereal *, doublereal *, integer *, integer *, doublereal *,
                       integer *, doublereal *, integer *, integer *);
    integer iqptr;
    extern int dlaed7_(integer *, integer *, integer *, integer *, integer *, integer *,
                       doublereal *, doublereal *, integer *, integer *, doublereal *, integer *,
                       doublereal *, integer *, integer *, integer *, integer *, integer *,
                       doublereal *, doublereal *, integer *, integer *);
    integer tlvls;
    extern int dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *,
                       integer *, ftnlen);
    integer igivcl;
    extern int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    integer igivnm, submat, curprb, subpbs, igivpt;
    extern int dsteqr_(char *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                       doublereal *, integer *, ftnlen);
    integer curlvl, matsiz, iprmpt, smlsiz;
    --d__;
    --e;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    qstore_dim1 = *ldqs;
    qstore_offset = 1 + qstore_dim1;
    qstore -= qstore_offset;
    --work;
    --iwork;
    *info = 0;
    if (*icompq < 0 || *icompq > 2) {
        *info = -1;
    } else if (*icompq == 1 && *qsiz < max(0, *n)) {
        *info = -2;
    } else if (*n < 0) {
        *info = -3;
    } else if (*ldq < max(1, *n)) {
        *info = -7;
    } else if (*ldqs < max(1, *n)) {
        *info = -9;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLAED0", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    smlsiz = ilaenv_(&c__9, (char *)"DLAED0", (char *)" ", &c__0, &c__0, &c__0, &c__0, (ftnlen)6, (ftnlen)1);
    iwork[1] = *n;
    subpbs = 1;
    tlvls = 0;
L10:
    if (iwork[subpbs] > smlsiz) {
        for (j = subpbs; j >= 1; --j) {
            iwork[j * 2] = (iwork[j] + 1) / 2;
            iwork[(j << 1) - 1] = iwork[j] / 2;
        }
        ++tlvls;
        subpbs <<= 1;
        goto L10;
    }
    i__1 = subpbs;
    for (j = 2; j <= i__1; ++j) {
        iwork[j] += iwork[j - 1];
    }
    spm1 = subpbs - 1;
    i__1 = spm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
        submat = iwork[i__] + 1;
        smm1 = submat - 1;
        d__[smm1] -= (d__1 = e[smm1], abs(d__1));
        d__[submat] -= (d__1 = e[smm1], abs(d__1));
    }
    indxq = (*n << 2) + 3;
    if (*icompq != 2) {
        temp = log((doublereal)(*n)) / log(2.);
        lgn = (integer)temp;
        if (pow_lmp_ii(&c__2, &lgn) < *n) {
            ++lgn;
        }
        if (pow_lmp_ii(&c__2, &lgn) < *n) {
            ++lgn;
        }
        iprmpt = indxq + *n + 1;
        iperm = iprmpt + *n * lgn;
        iqptr = iperm + *n * lgn;
        igivpt = iqptr + *n + 2;
        igivcl = igivpt + *n * lgn;
        igivnm = 1;
        iq = igivnm + (*n << 1) * lgn;
        i__1 = *n;
        iwrem = iq + i__1 * i__1 + 1;
        i__1 = subpbs;
        for (i__ = 0; i__ <= i__1; ++i__) {
            iwork[iprmpt + i__] = 1;
            iwork[igivpt + i__] = 1;
        }
        iwork[iqptr] = 1;
    }
    curr = 0;
    i__1 = spm1;
    for (i__ = 0; i__ <= i__1; ++i__) {
        if (i__ == 0) {
            submat = 1;
            matsiz = iwork[1];
        } else {
            submat = iwork[i__] + 1;
            matsiz = iwork[i__ + 1] - iwork[i__];
        }
        if (*icompq == 2) {
            dsteqr_((char *)"I", &matsiz, &d__[submat], &e[submat], &q[submat + submat * q_dim1], ldq,
                    &work[1], info, (ftnlen)1);
            if (*info != 0) {
                goto L130;
            }
        } else {
            dsteqr_((char *)"I", &matsiz, &d__[submat], &e[submat], &work[iq - 1 + iwork[iqptr + curr]],
                    &matsiz, &work[1], info, (ftnlen)1);
            if (*info != 0) {
                goto L130;
            }
            if (*icompq == 1) {
                dgemm_((char *)"N", (char *)"N", qsiz, &matsiz, &matsiz, &c_b23, &q[submat * q_dim1 + 1], ldq,
                       &work[iq - 1 + iwork[iqptr + curr]], &matsiz, &c_b24,
                       &qstore[submat * qstore_dim1 + 1], ldqs, (ftnlen)1, (ftnlen)1);
            }
            i__2 = matsiz;
            iwork[iqptr + curr + 1] = iwork[iqptr + curr] + i__2 * i__2;
            ++curr;
        }
        k = 1;
        i__2 = iwork[i__ + 1];
        for (j = submat; j <= i__2; ++j) {
            iwork[indxq + j] = k;
            ++k;
        }
    }
    curlvl = 1;
L80:
    if (subpbs > 1) {
        spm2 = subpbs - 2;
        i__1 = spm2;
        for (i__ = 0; i__ <= i__1; i__ += 2) {
            if (i__ == 0) {
                submat = 1;
                matsiz = iwork[2];
                msd2 = iwork[1];
                curprb = 0;
            } else {
                submat = iwork[i__] + 1;
                matsiz = iwork[i__ + 2] - iwork[i__];
                msd2 = matsiz / 2;
                ++curprb;
            }
            if (*icompq == 2) {
                dlaed1_(&matsiz, &d__[submat], &q[submat + submat * q_dim1], ldq,
                        &iwork[indxq + submat], &e[submat + msd2 - 1], &msd2, &work[1],
                        &iwork[subpbs + 1], info);
            } else {
                dlaed7_(icompq, &matsiz, qsiz, &tlvls, &curlvl, &curprb, &d__[submat],
                        &qstore[submat * qstore_dim1 + 1], ldqs, &iwork[indxq + submat],
                        &e[submat + msd2 - 1], &msd2, &work[iq], &iwork[iqptr], &iwork[iprmpt],
                        &iwork[iperm], &iwork[igivpt], &iwork[igivcl], &work[igivnm], &work[iwrem],
                        &iwork[subpbs + 1], info);
            }
            if (*info != 0) {
                goto L130;
            }
            iwork[i__ / 2 + 1] = iwork[i__ + 2];
        }
        subpbs /= 2;
        ++curlvl;
        goto L80;
    }
    if (*icompq == 1) {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            j = iwork[indxq + i__];
            work[i__] = d__[j];
            dcopy_(qsiz, &qstore[j * qstore_dim1 + 1], &c__1, &q[i__ * q_dim1 + 1], &c__1);
        }
        dcopy_(n, &work[1], &c__1, &d__[1], &c__1);
    } else if (*icompq == 2) {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            j = iwork[indxq + i__];
            work[i__] = d__[j];
            dcopy_(n, &q[j * q_dim1 + 1], &c__1, &work[*n * i__ + 1], &c__1);
        }
        dcopy_(n, &work[1], &c__1, &d__[1], &c__1);
        dlacpy_((char *)"A", n, n, &work[*n + 1], n, &q[q_offset], ldq, (ftnlen)1);
    } else {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            j = iwork[indxq + i__];
            work[i__] = d__[j];
        }
        dcopy_(n, &work[1], &c__1, &d__[1], &c__1);
    }
    goto L140;
L130:
    *info = submat * (*n + 1) + submat + matsiz - 1;
L140:
    return 0;
}
#ifdef __cplusplus
}
#endif
