#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__9 = 9;
static integer c__0 = 0;
static integer c__2 = 2;
static integer c__1 = 1;
int zlaed0_(integer *qsiz, integer *n, doublereal *d__, doublereal *e, doublecomplex *q,
            integer *ldq, doublecomplex *qstore, integer *ldqs, doublereal *rwork, integer *iwork,
            integer *info)
{
    integer q_dim1, q_offset, qstore_dim1, qstore_offset, i__1, i__2;
    doublereal d__1;
    double log(doublereal);
    integer pow_lmp_ii(integer *, integer *);
    integer i__, j, k, ll, iq, lgn, msd2, smm1, spm1, spm2;
    doublereal temp;
    integer curr, iperm;
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer indxq, iwrem, iqptr, tlvls;
    extern int zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, integer *),
        zlaed7_(integer *, integer *, integer *, integer *, integer *, integer *, doublereal *,
                doublecomplex *, integer *, doublereal *, integer *, doublereal *, integer *,
                integer *, integer *, integer *, integer *, doublereal *, doublecomplex *,
                doublereal *, integer *, integer *);
    integer igivcl;
    extern int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    extern int zlacrm_(integer *, integer *, doublecomplex *, integer *, doublereal *, integer *,
                       doublecomplex *, integer *, doublereal *);
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
    --rwork;
    --iwork;
    *info = 0;
    if (*qsiz < max(0, *n)) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*ldq < max(1, *n)) {
        *info = -6;
    } else if (*ldqs < max(1, *n)) {
        *info = -8;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"ZLAED0", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    smlsiz = ilaenv_(&c__9, (char *)"ZLAED0", (char *)" ", &c__0, &c__0, &c__0, &c__0, (ftnlen)6, (ftnlen)1);
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
        ll = iq - 1 + iwork[iqptr + curr];
        dsteqr_((char *)"I", &matsiz, &d__[submat], &e[submat], &rwork[ll], &matsiz, &rwork[1], info,
                (ftnlen)1);
        zlacrm_(qsiz, &matsiz, &q[submat * q_dim1 + 1], ldq, &rwork[ll], &matsiz,
                &qstore[submat * qstore_dim1 + 1], ldqs, &rwork[iwrem]);
        i__2 = matsiz;
        iwork[iqptr + curr + 1] = iwork[iqptr + curr] + i__2 * i__2;
        ++curr;
        if (*info > 0) {
            *info = submat * (*n + 1) + submat + matsiz - 1;
            return 0;
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
            zlaed7_(&matsiz, &msd2, qsiz, &tlvls, &curlvl, &curprb, &d__[submat],
                    &qstore[submat * qstore_dim1 + 1], ldqs, &e[submat + msd2 - 1],
                    &iwork[indxq + submat], &rwork[iq], &iwork[iqptr], &iwork[iprmpt],
                    &iwork[iperm], &iwork[igivpt], &iwork[igivcl], &rwork[igivnm],
                    &q[submat * q_dim1 + 1], &rwork[iwrem], &iwork[subpbs + 1], info);
            if (*info > 0) {
                *info = submat * (*n + 1) + submat + matsiz - 1;
                return 0;
            }
            iwork[i__ / 2 + 1] = iwork[i__ + 2];
        }
        subpbs /= 2;
        ++curlvl;
        goto L80;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        j = iwork[indxq + i__];
        rwork[i__] = d__[j];
        zcopy_(qsiz, &qstore[j * qstore_dim1 + 1], &c__1, &q[i__ * q_dim1 + 1], &c__1);
    }
    dcopy_(n, &rwork[1], &c__1, &d__[1], &c__1);
    return 0;
}
#ifdef __cplusplus
}
#endif
