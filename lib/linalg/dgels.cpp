#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b33 = 0.;
static integer c__0 = 0;
int dgels_(char *trans, integer *m, integer *n, integer *nrhs, doublereal *a, integer *lda,
           doublereal *b, integer *ldb, doublereal *work, integer *lwork, integer *info,
           ftnlen trans_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    integer i__, j, nb, mn;
    doublereal anrm, bnrm;
    integer brow;
    logical tpsd;
    integer iascl, ibscl;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer wsize;
    doublereal rwork[1];
    extern doublereal dlamch_(char *, ftnlen),
        dlange_(char *, integer *, integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern int dgelqf_(integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                       integer *, integer *),
        dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *,
                doublereal *, integer *, integer *, ftnlen),
        dgeqrf_(integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                integer *, integer *),
        dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                ftnlen),
        xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    integer scllen;
    doublereal bignum;
    extern int dormlq_(char *, char *, integer *, integer *, integer *, doublereal *, integer *,
                       doublereal *, doublereal *, integer *, doublereal *, integer *, integer *,
                       ftnlen, ftnlen),
        dormqr_(char *, char *, integer *, integer *, integer *, doublereal *, integer *,
                doublereal *, doublereal *, integer *, doublereal *, integer *, integer *, ftnlen,
                ftnlen);
    doublereal smlnum;
    logical lquery;
    extern int dtrtrs_(char *, char *, char *, integer *, integer *, doublereal *, integer *,
                       doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --work;
    *info = 0;
    mn = min(*m, *n);
    lquery = *lwork == -1;
    if (!(lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1) || lsame_(trans, (char *)"T", (ftnlen)1, (ftnlen)1))) {
        *info = -1;
    } else if (*m < 0) {
        *info = -2;
    } else if (*n < 0) {
        *info = -3;
    } else if (*nrhs < 0) {
        *info = -4;
    } else if (*lda < max(1, *m)) {
        *info = -6;
    } else {
        i__1 = max(1, *m);
        if (*ldb < max(i__1, *n)) {
            *info = -8;
        } else {
            i__1 = 1, i__2 = mn + max(mn, *nrhs);
            if (*lwork < max(i__1, i__2) && !lquery) {
                *info = -10;
            }
        }
    }
    if (*info == 0 || *info == -10) {
        tpsd = TRUE_;
        if (lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1)) {
            tpsd = FALSE_;
        }
        if (*m >= *n) {
            nb = ilaenv_(&c__1, (char *)"DGEQRF", (char *)" ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
            if (tpsd) {
                i__1 = nb,
                i__2 = ilaenv_(&c__1, (char *)"DORMQR", (char *)"LN", m, nrhs, n, &c_n1, (ftnlen)6, (ftnlen)2);
                nb = max(i__1, i__2);
            } else {
                i__1 = nb,
                i__2 = ilaenv_(&c__1, (char *)"DORMQR", (char *)"LT", m, nrhs, n, &c_n1, (ftnlen)6, (ftnlen)2);
                nb = max(i__1, i__2);
            }
        } else {
            nb = ilaenv_(&c__1, (char *)"DGELQF", (char *)" ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
            if (tpsd) {
                i__1 = nb,
                i__2 = ilaenv_(&c__1, (char *)"DORMLQ", (char *)"LT", n, nrhs, m, &c_n1, (ftnlen)6, (ftnlen)2);
                nb = max(i__1, i__2);
            } else {
                i__1 = nb,
                i__2 = ilaenv_(&c__1, (char *)"DORMLQ", (char *)"LN", n, nrhs, m, &c_n1, (ftnlen)6, (ftnlen)2);
                nb = max(i__1, i__2);
            }
        }
        i__1 = 1, i__2 = mn + max(mn, *nrhs) * nb;
        wsize = max(i__1, i__2);
        work[1] = (doublereal)wsize;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DGELS ", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    i__1 = min(*m, *n);
    if (min(i__1, *nrhs) == 0) {
        i__1 = max(*m, *n);
        dlaset_((char *)"Full", &i__1, nrhs, &c_b33, &c_b33, &b[b_offset], ldb, (ftnlen)4);
        return 0;
    }
    smlnum = dlamch_((char *)"S", (ftnlen)1) / dlamch_((char *)"P", (ftnlen)1);
    bignum = 1. / smlnum;
    anrm = dlange_((char *)"M", m, n, &a[a_offset], lda, rwork, (ftnlen)1);
    iascl = 0;
    if (anrm > 0. && anrm < smlnum) {
        dlascl_((char *)"G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, info, (ftnlen)1);
        iascl = 1;
    } else if (anrm > bignum) {
        dlascl_((char *)"G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, info, (ftnlen)1);
        iascl = 2;
    } else if (anrm == 0.) {
        i__1 = max(*m, *n);
        dlaset_((char *)"F", &i__1, nrhs, &c_b33, &c_b33, &b[b_offset], ldb, (ftnlen)1);
        goto L50;
    }
    brow = *m;
    if (tpsd) {
        brow = *n;
    }
    bnrm = dlange_((char *)"M", &brow, nrhs, &b[b_offset], ldb, rwork, (ftnlen)1);
    ibscl = 0;
    if (bnrm > 0. && bnrm < smlnum) {
        dlascl_((char *)"G", &c__0, &c__0, &bnrm, &smlnum, &brow, nrhs, &b[b_offset], ldb, info, (ftnlen)1);
        ibscl = 1;
    } else if (bnrm > bignum) {
        dlascl_((char *)"G", &c__0, &c__0, &bnrm, &bignum, &brow, nrhs, &b[b_offset], ldb, info, (ftnlen)1);
        ibscl = 2;
    }
    if (*m >= *n) {
        i__1 = *lwork - mn;
        dgeqrf_(m, n, &a[a_offset], lda, &work[1], &work[mn + 1], &i__1, info);
        if (!tpsd) {
            i__1 = *lwork - mn;
            dormqr_((char *)"Left", (char *)"Transpose", m, nrhs, n, &a[a_offset], lda, &work[1], &b[b_offset], ldb,
                    &work[mn + 1], &i__1, info, (ftnlen)4, (ftnlen)9);
            dtrtrs_((char *)"Upper", (char *)"No transpose", (char *)"Non-unit", n, nrhs, &a[a_offset], lda, &b[b_offset],
                    ldb, info, (ftnlen)5, (ftnlen)12, (ftnlen)8);
            if (*info > 0) {
                return 0;
            }
            scllen = *n;
        } else {
            dtrtrs_((char *)"Upper", (char *)"Transpose", (char *)"Non-unit", n, nrhs, &a[a_offset], lda, &b[b_offset], ldb,
                    info, (ftnlen)5, (ftnlen)9, (ftnlen)8);
            if (*info > 0) {
                return 0;
            }
            i__1 = *nrhs;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *m;
                for (i__ = *n + 1; i__ <= i__2; ++i__) {
                    b[i__ + j * b_dim1] = 0.;
                }
            }
            i__1 = *lwork - mn;
            dormqr_((char *)"Left", (char *)"No transpose", m, nrhs, n, &a[a_offset], lda, &work[1], &b[b_offset],
                    ldb, &work[mn + 1], &i__1, info, (ftnlen)4, (ftnlen)12);
            scllen = *m;
        }
    } else {
        i__1 = *lwork - mn;
        dgelqf_(m, n, &a[a_offset], lda, &work[1], &work[mn + 1], &i__1, info);
        if (!tpsd) {
            dtrtrs_((char *)"Lower", (char *)"No transpose", (char *)"Non-unit", m, nrhs, &a[a_offset], lda, &b[b_offset],
                    ldb, info, (ftnlen)5, (ftnlen)12, (ftnlen)8);
            if (*info > 0) {
                return 0;
            }
            i__1 = *nrhs;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *n;
                for (i__ = *m + 1; i__ <= i__2; ++i__) {
                    b[i__ + j * b_dim1] = 0.;
                }
            }
            i__1 = *lwork - mn;
            dormlq_((char *)"Left", (char *)"Transpose", n, nrhs, m, &a[a_offset], lda, &work[1], &b[b_offset], ldb,
                    &work[mn + 1], &i__1, info, (ftnlen)4, (ftnlen)9);
            scllen = *n;
        } else {
            i__1 = *lwork - mn;
            dormlq_((char *)"Left", (char *)"No transpose", n, nrhs, m, &a[a_offset], lda, &work[1], &b[b_offset],
                    ldb, &work[mn + 1], &i__1, info, (ftnlen)4, (ftnlen)12);
            dtrtrs_((char *)"Lower", (char *)"Transpose", (char *)"Non-unit", m, nrhs, &a[a_offset], lda, &b[b_offset], ldb,
                    info, (ftnlen)5, (ftnlen)9, (ftnlen)8);
            if (*info > 0) {
                return 0;
            }
            scllen = *m;
        }
    }
    if (iascl == 1) {
        dlascl_((char *)"G", &c__0, &c__0, &anrm, &smlnum, &scllen, nrhs, &b[b_offset], ldb, info,
                (ftnlen)1);
    } else if (iascl == 2) {
        dlascl_((char *)"G", &c__0, &c__0, &anrm, &bignum, &scllen, nrhs, &b[b_offset], ldb, info,
                (ftnlen)1);
    }
    if (ibscl == 1) {
        dlascl_((char *)"G", &c__0, &c__0, &smlnum, &bnrm, &scllen, nrhs, &b[b_offset], ldb, info,
                (ftnlen)1);
    } else if (ibscl == 2) {
        dlascl_((char *)"G", &c__0, &c__0, &bignum, &bnrm, &scllen, nrhs, &b[b_offset], ldb, info,
                (ftnlen)1);
    }
L50:
    work[1] = (doublereal)wsize;
    return 0;
}
#ifdef __cplusplus
}
#endif
