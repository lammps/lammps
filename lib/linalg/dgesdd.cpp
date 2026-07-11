#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b63 = 0.;
static integer c__1 = 1;
static doublereal c_b84 = 1.;
int dgesdd_(char *jobz, integer *m, integer *n, doublereal *a, integer *lda, doublereal *s,
            doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, doublereal *work,
            integer *lwork, integer *iwork, integer *info, ftnlen jobz_len)
{
    integer a_dim1, a_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2, i__3;
    double sqrt(doublereal);
    integer lwork_dorglq_mn__, lwork_dorglq_nn__, lwork_dorgqr_mm__, lwork_dorgqr_mn__, i__, ie,
        lwork_dorgbr_p_mm__, il, lwork_dorgbr_q_nn__, ir, iu, blk;
    doublereal dum[1], eps;
    integer ivt, iscl;
    doublereal anrm;
    integer idum[1], ierr, itau, lwork_dormbr_qln_mm__, lwork_dormbr_qln_mn__,
        lwork_dormbr_qln_nn__, lwork_dormbr_prt_mm__, lwork_dormbr_prt_mn__, lwork_dormbr_prt_nn__;
    extern int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
                      integer *, doublereal *, integer *, doublereal *, doublereal *, integer *,
                      ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer chunk, minmn, wrkbl, itaup, itauq, mnthr;
    logical wntqa;
    integer nwork;
    logical wntqn, wntqo, wntqs;
    extern int dbdsdc_(char *, char *, integer *, doublereal *, doublereal *, doublereal *,
                       integer *, doublereal *, integer *, doublereal *, integer *, doublereal *,
                       integer *, integer *, ftnlen, ftnlen),
        dgebrd_(integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                doublereal *, doublereal *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen),
        dlange_(char *, integer *, integer *, doublereal *, integer *, doublereal *, ftnlen);
    integer bdspac;
    extern int dgelqf_(integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                       integer *, integer *),
        dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *,
                doublereal *, integer *, integer *, ftnlen),
        dgeqrf_(integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                integer *, integer *),
        dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *,
                ftnlen),
        dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                ftnlen),
        xerbla_(char *, integer *, ftnlen),
        dorgbr_(char *, integer *, integer *, integer *, doublereal *, integer *, doublereal *,
                doublereal *, integer *, integer *, ftnlen);
    extern logical disnan_(doublereal *);
    doublereal bignum;
    extern int dormbr_(char *, char *, char *, integer *, integer *, integer *, doublereal *,
                       integer *, doublereal *, doublereal *, integer *, doublereal *, integer *,
                       integer *, ftnlen, ftnlen, ftnlen),
        dorglq_(integer *, integer *, integer *, doublereal *, integer *, doublereal *,
                doublereal *, integer *, integer *),
        dorgqr_(integer *, integer *, integer *, doublereal *, integer *, doublereal *,
                doublereal *, integer *, integer *);
    integer ldwrkl, ldwrkr, minwrk, ldwrku, maxwrk, ldwkvt;
    doublereal smlnum;
    logical wntqas, lquery;
    extern doublereal droundup_lwork__(integer *);
    integer lwork_dgebrd_mm__, lwork_dgebrd_mn__, lwork_dgebrd_nn__, lwork_dgelqf_mn__,
        lwork_dgeqrf_mn__;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --s;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --work;
    --iwork;
    *info = 0;
    minmn = min(*m, *n);
    wntqa = lsame_(jobz, (char *)"A", (ftnlen)1, (ftnlen)1);
    wntqs = lsame_(jobz, (char *)"S", (ftnlen)1, (ftnlen)1);
    wntqas = wntqa || wntqs;
    wntqo = lsame_(jobz, (char *)"O", (ftnlen)1, (ftnlen)1);
    wntqn = lsame_(jobz, (char *)"N", (ftnlen)1, (ftnlen)1);
    lquery = *lwork == -1;
    if (!(wntqa || wntqs || wntqo || wntqn)) {
        *info = -1;
    } else if (*m < 0) {
        *info = -2;
    } else if (*n < 0) {
        *info = -3;
    } else if (*lda < max(1, *m)) {
        *info = -5;
    } else if (*ldu < 1 || wntqas && *ldu < *m || wntqo && *m < *n && *ldu < *m) {
        *info = -8;
    } else if (*ldvt < 1 || wntqa && *ldvt < *n || wntqs && *ldvt < minmn ||
               wntqo && *m >= *n && *ldvt < *n) {
        *info = -10;
    }
    if (*info == 0) {
        minwrk = 1;
        maxwrk = 1;
        bdspac = 0;
        mnthr = (integer)(minmn * 11. / 6.);
        if (*m >= *n && minmn > 0) {
            if (wntqn) {
                bdspac = *n * 7;
            } else {
                bdspac = *n * 3 * *n + (*n << 2);
            }
            dgebrd_(m, n, dum, m, dum, dum, dum, dum, dum, &c_n1, &ierr);
            lwork_dgebrd_mn__ = (integer)dum[0];
            dgebrd_(n, n, dum, n, dum, dum, dum, dum, dum, &c_n1, &ierr);
            lwork_dgebrd_nn__ = (integer)dum[0];
            dgeqrf_(m, n, dum, m, dum, dum, &c_n1, &ierr);
            lwork_dgeqrf_mn__ = (integer)dum[0];
            dorgbr_((char *)"Q", n, n, n, dum, n, dum, dum, &c_n1, &ierr, (ftnlen)1);
            lwork_dorgbr_q_nn__ = (integer)dum[0];
            dorgqr_(m, m, n, dum, m, dum, dum, &c_n1, &ierr);
            lwork_dorgqr_mm__ = (integer)dum[0];
            dorgqr_(m, n, n, dum, m, dum, dum, &c_n1, &ierr);
            lwork_dorgqr_mn__ = (integer)dum[0];
            dormbr_((char *)"P", (char *)"R", (char *)"T", n, n, n, dum, n, dum, dum, n, dum, &c_n1, &ierr, (ftnlen)1,
                    (ftnlen)1, (ftnlen)1);
            lwork_dormbr_prt_nn__ = (integer)dum[0];
            dormbr_((char *)"Q", (char *)"L", (char *)"N", n, n, n, dum, n, dum, dum, n, dum, &c_n1, &ierr, (ftnlen)1,
                    (ftnlen)1, (ftnlen)1);
            lwork_dormbr_qln_nn__ = (integer)dum[0];
            dormbr_((char *)"Q", (char *)"L", (char *)"N", m, n, n, dum, m, dum, dum, m, dum, &c_n1, &ierr, (ftnlen)1,
                    (ftnlen)1, (ftnlen)1);
            lwork_dormbr_qln_mn__ = (integer)dum[0];
            dormbr_((char *)"Q", (char *)"L", (char *)"N", m, m, n, dum, m, dum, dum, m, dum, &c_n1, &ierr, (ftnlen)1,
                    (ftnlen)1, (ftnlen)1);
            lwork_dormbr_qln_mm__ = (integer)dum[0];
            if (*m >= mnthr) {
                if (wntqn) {
                    wrkbl = *n + lwork_dgeqrf_mn__;
                    i__1 = wrkbl, i__2 = *n * 3 + lwork_dgebrd_nn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = bdspac + *n;
                    maxwrk = max(i__1, i__2);
                    minwrk = bdspac + *n;
                } else if (wntqo) {
                    wrkbl = *n + lwork_dgeqrf_mn__;
                    i__1 = wrkbl, i__2 = *n + lwork_dorgqr_mn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *n * 3 + lwork_dgebrd_nn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_qln_nn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_prt_nn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *n * 3 + bdspac;
                    wrkbl = max(i__1, i__2);
                    maxwrk = wrkbl + (*n << 1) * *n;
                    minwrk = bdspac + (*n << 1) * *n + *n * 3;
                } else if (wntqs) {
                    wrkbl = *n + lwork_dgeqrf_mn__;
                    i__1 = wrkbl, i__2 = *n + lwork_dorgqr_mn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *n * 3 + lwork_dgebrd_nn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_qln_nn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_prt_nn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *n * 3 + bdspac;
                    wrkbl = max(i__1, i__2);
                    maxwrk = wrkbl + *n * *n;
                    minwrk = bdspac + *n * *n + *n * 3;
                } else if (wntqa) {
                    wrkbl = *n + lwork_dgeqrf_mn__;
                    i__1 = wrkbl, i__2 = *n + lwork_dorgqr_mm__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *n * 3 + lwork_dgebrd_nn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_qln_nn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_prt_nn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *n * 3 + bdspac;
                    wrkbl = max(i__1, i__2);
                    maxwrk = wrkbl + *n * *n;
                    i__1 = *n * 3 + bdspac, i__2 = *n + *m;
                    minwrk = *n * *n + max(i__1, i__2);
                }
            } else {
                wrkbl = *n * 3 + lwork_dgebrd_mn__;
                if (wntqn) {
                    i__1 = wrkbl, i__2 = *n * 3 + bdspac;
                    maxwrk = max(i__1, i__2);
                    minwrk = *n * 3 + max(*m, bdspac);
                } else if (wntqo) {
                    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_prt_nn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_qln_mn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *n * 3 + bdspac;
                    wrkbl = max(i__1, i__2);
                    maxwrk = wrkbl + *m * *n;
                    i__1 = *m, i__2 = *n * *n + bdspac;
                    minwrk = *n * 3 + max(i__1, i__2);
                } else if (wntqs) {
                    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_qln_mn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_prt_nn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *n * 3 + bdspac;
                    maxwrk = max(i__1, i__2);
                    minwrk = *n * 3 + max(*m, bdspac);
                } else if (wntqa) {
                    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_qln_mm__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *n * 3 + lwork_dormbr_prt_nn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *n * 3 + bdspac;
                    maxwrk = max(i__1, i__2);
                    minwrk = *n * 3 + max(*m, bdspac);
                }
            }
        } else if (minmn > 0) {
            if (wntqn) {
                bdspac = *m * 7;
            } else {
                bdspac = *m * 3 * *m + (*m << 2);
            }
            dgebrd_(m, n, dum, m, dum, dum, dum, dum, dum, &c_n1, &ierr);
            lwork_dgebrd_mn__ = (integer)dum[0];
            dgebrd_(m, m, &a[a_offset], m, &s[1], dum, dum, dum, dum, &c_n1, &ierr);
            lwork_dgebrd_mm__ = (integer)dum[0];
            dgelqf_(m, n, &a[a_offset], m, dum, dum, &c_n1, &ierr);
            lwork_dgelqf_mn__ = (integer)dum[0];
            dorglq_(n, n, m, dum, n, dum, dum, &c_n1, &ierr);
            lwork_dorglq_nn__ = (integer)dum[0];
            dorglq_(m, n, m, &a[a_offset], m, dum, dum, &c_n1, &ierr);
            lwork_dorglq_mn__ = (integer)dum[0];
            dorgbr_((char *)"P", m, m, m, &a[a_offset], n, dum, dum, &c_n1, &ierr, (ftnlen)1);
            lwork_dorgbr_p_mm__ = (integer)dum[0];
            dormbr_((char *)"P", (char *)"R", (char *)"T", m, m, m, dum, m, dum, dum, m, dum, &c_n1, &ierr, (ftnlen)1,
                    (ftnlen)1, (ftnlen)1);
            lwork_dormbr_prt_mm__ = (integer)dum[0];
            dormbr_((char *)"P", (char *)"R", (char *)"T", m, n, m, dum, m, dum, dum, m, dum, &c_n1, &ierr, (ftnlen)1,
                    (ftnlen)1, (ftnlen)1);
            lwork_dormbr_prt_mn__ = (integer)dum[0];
            dormbr_((char *)"P", (char *)"R", (char *)"T", n, n, m, dum, n, dum, dum, n, dum, &c_n1, &ierr, (ftnlen)1,
                    (ftnlen)1, (ftnlen)1);
            lwork_dormbr_prt_nn__ = (integer)dum[0];
            dormbr_((char *)"Q", (char *)"L", (char *)"N", m, m, m, dum, m, dum, dum, m, dum, &c_n1, &ierr, (ftnlen)1,
                    (ftnlen)1, (ftnlen)1);
            lwork_dormbr_qln_mm__ = (integer)dum[0];
            if (*n >= mnthr) {
                if (wntqn) {
                    wrkbl = *m + lwork_dgelqf_mn__;
                    i__1 = wrkbl, i__2 = *m * 3 + lwork_dgebrd_mm__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = bdspac + *m;
                    maxwrk = max(i__1, i__2);
                    minwrk = bdspac + *m;
                } else if (wntqo) {
                    wrkbl = *m + lwork_dgelqf_mn__;
                    i__1 = wrkbl, i__2 = *m + lwork_dorglq_mn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *m * 3 + lwork_dgebrd_mm__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_qln_mm__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_prt_mm__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *m * 3 + bdspac;
                    wrkbl = max(i__1, i__2);
                    maxwrk = wrkbl + (*m << 1) * *m;
                    minwrk = bdspac + (*m << 1) * *m + *m * 3;
                } else if (wntqs) {
                    wrkbl = *m + lwork_dgelqf_mn__;
                    i__1 = wrkbl, i__2 = *m + lwork_dorglq_mn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *m * 3 + lwork_dgebrd_mm__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_qln_mm__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_prt_mm__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *m * 3 + bdspac;
                    wrkbl = max(i__1, i__2);
                    maxwrk = wrkbl + *m * *m;
                    minwrk = bdspac + *m * *m + *m * 3;
                } else if (wntqa) {
                    wrkbl = *m + lwork_dgelqf_mn__;
                    i__1 = wrkbl, i__2 = *m + lwork_dorglq_nn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *m * 3 + lwork_dgebrd_mm__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_qln_mm__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_prt_mm__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *m * 3 + bdspac;
                    wrkbl = max(i__1, i__2);
                    maxwrk = wrkbl + *m * *m;
                    i__1 = *m * 3 + bdspac, i__2 = *m + *n;
                    minwrk = *m * *m + max(i__1, i__2);
                }
            } else {
                wrkbl = *m * 3 + lwork_dgebrd_mn__;
                if (wntqn) {
                    i__1 = wrkbl, i__2 = *m * 3 + bdspac;
                    maxwrk = max(i__1, i__2);
                    minwrk = *m * 3 + max(*n, bdspac);
                } else if (wntqo) {
                    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_qln_mm__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_prt_mn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *m * 3 + bdspac;
                    wrkbl = max(i__1, i__2);
                    maxwrk = wrkbl + *m * *n;
                    i__1 = *n, i__2 = *m * *m + bdspac;
                    minwrk = *m * 3 + max(i__1, i__2);
                } else if (wntqs) {
                    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_qln_mm__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_prt_mn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *m * 3 + bdspac;
                    maxwrk = max(i__1, i__2);
                    minwrk = *m * 3 + max(*n, bdspac);
                } else if (wntqa) {
                    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_qln_mm__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *m * 3 + lwork_dormbr_prt_nn__;
                    wrkbl = max(i__1, i__2);
                    i__1 = wrkbl, i__2 = *m * 3 + bdspac;
                    maxwrk = max(i__1, i__2);
                    minwrk = *m * 3 + max(*n, bdspac);
                }
            }
        }
        maxwrk = max(maxwrk, minwrk);
        work[1] = droundup_lwork__(&maxwrk);
        if (*lwork < minwrk && !lquery) {
            *info = -12;
        }
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DGESDD", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    if (*m == 0 || *n == 0) {
        return 0;
    }
    eps = dlamch_((char *)"P", (ftnlen)1);
    smlnum = sqrt(dlamch_((char *)"S", (ftnlen)1)) / eps;
    bignum = 1. / smlnum;
    anrm = dlange_((char *)"M", m, n, &a[a_offset], lda, dum, (ftnlen)1);
    if (disnan_(&anrm)) {
        *info = -4;
        return 0;
    }
    iscl = 0;
    if (anrm > 0. && anrm < smlnum) {
        iscl = 1;
        dlascl_((char *)"G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, &ierr, (ftnlen)1);
    } else if (anrm > bignum) {
        iscl = 1;
        dlascl_((char *)"G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, &ierr, (ftnlen)1);
    }
    if (*m >= *n) {
        if (*m >= mnthr) {
            if (wntqn) {
                itau = 1;
                nwork = itau + *n;
                i__1 = *lwork - nwork + 1;
                dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &i__1, &ierr);
                i__1 = *n - 1;
                i__2 = *n - 1;
                dlaset_((char *)"L", &i__1, &i__2, &c_b63, &c_b63, &a[a_dim1 + 2], lda, (ftnlen)1);
                ie = 1;
                itauq = ie + *n;
                itaup = itauq + *n;
                nwork = itaup + *n;
                i__1 = *lwork - nwork + 1;
                dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &work[itaup],
                        &work[nwork], &i__1, &ierr);
                nwork = ie + *n;
                dbdsdc_((char *)"U", (char *)"N", n, &s[1], &work[ie], dum, &c__1, dum, &c__1, dum, idum,
                        &work[nwork], &iwork[1], info, (ftnlen)1, (ftnlen)1);
            } else if (wntqo) {
                ir = 1;
                if (*lwork >= *lda * *n + *n * *n + *n * 3 + bdspac) {
                    ldwrkr = *lda;
                } else {
                    ldwrkr = (*lwork - *n * *n - *n * 3 - bdspac) / *n;
                }
                itau = ir + ldwrkr * *n;
                nwork = itau + *n;
                i__1 = *lwork - nwork + 1;
                dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &i__1, &ierr);
                dlacpy_((char *)"U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, (ftnlen)1);
                i__1 = *n - 1;
                i__2 = *n - 1;
                dlaset_((char *)"L", &i__1, &i__2, &c_b63, &c_b63, &work[ir + 1], &ldwrkr, (ftnlen)1);
                i__1 = *lwork - nwork + 1;
                dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[nwork], &i__1, &ierr);
                ie = itau;
                itauq = ie + *n;
                itaup = itauq + *n;
                nwork = itaup + *n;
                i__1 = *lwork - nwork + 1;
                dgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &work[itauq], &work[itaup],
                        &work[nwork], &i__1, &ierr);
                iu = nwork;
                nwork = iu + *n * *n;
                dbdsdc_((char *)"U", (char *)"I", n, &s[1], &work[ie], &work[iu], n, &vt[vt_offset], ldvt, dum,
                        idum, &work[nwork], &iwork[1], info, (ftnlen)1, (ftnlen)1);
                i__1 = *lwork - nwork + 1;
                dormbr_((char *)"Q", (char *)"L", (char *)"N", n, n, n, &work[ir], &ldwrkr, &work[itauq], &work[iu], n,
                        &work[nwork], &i__1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                i__1 = *lwork - nwork + 1;
                dormbr_((char *)"P", (char *)"R", (char *)"T", n, n, n, &work[ir], &ldwrkr, &work[itaup], &vt[vt_offset],
                        ldvt, &work[nwork], &i__1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                i__1 = *m;
                i__2 = ldwrkr;
                for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
                    i__3 = *m - i__ + 1;
                    chunk = min(i__3, ldwrkr);
                    dgemm_((char *)"N", (char *)"N", &chunk, n, n, &c_b84, &a[i__ + a_dim1], lda, &work[iu], n,
                           &c_b63, &work[ir], &ldwrkr, (ftnlen)1, (ftnlen)1);
                    dlacpy_((char *)"F", &chunk, n, &work[ir], &ldwrkr, &a[i__ + a_dim1], lda, (ftnlen)1);
                }
            } else if (wntqs) {
                ir = 1;
                ldwrkr = *n;
                itau = ir + ldwrkr * *n;
                nwork = itau + *n;
                i__2 = *lwork - nwork + 1;
                dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &i__2, &ierr);
                dlacpy_((char *)"U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, (ftnlen)1);
                i__2 = *n - 1;
                i__1 = *n - 1;
                dlaset_((char *)"L", &i__2, &i__1, &c_b63, &c_b63, &work[ir + 1], &ldwrkr, (ftnlen)1);
                i__2 = *lwork - nwork + 1;
                dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[nwork], &i__2, &ierr);
                ie = itau;
                itauq = ie + *n;
                itaup = itauq + *n;
                nwork = itaup + *n;
                i__2 = *lwork - nwork + 1;
                dgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &work[itauq], &work[itaup],
                        &work[nwork], &i__2, &ierr);
                dbdsdc_((char *)"U", (char *)"I", n, &s[1], &work[ie], &u[u_offset], ldu, &vt[vt_offset], ldvt, dum,
                        idum, &work[nwork], &iwork[1], info, (ftnlen)1, (ftnlen)1);
                i__2 = *lwork - nwork + 1;
                dormbr_((char *)"Q", (char *)"L", (char *)"N", n, n, n, &work[ir], &ldwrkr, &work[itauq], &u[u_offset], ldu,
                        &work[nwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                i__2 = *lwork - nwork + 1;
                dormbr_((char *)"P", (char *)"R", (char *)"T", n, n, n, &work[ir], &ldwrkr, &work[itaup], &vt[vt_offset],
                        ldvt, &work[nwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                dlacpy_((char *)"F", n, n, &u[u_offset], ldu, &work[ir], &ldwrkr, (ftnlen)1);
                dgemm_((char *)"N", (char *)"N", m, n, n, &c_b84, &a[a_offset], lda, &work[ir], &ldwrkr, &c_b63,
                       &u[u_offset], ldu, (ftnlen)1, (ftnlen)1);
            } else if (wntqa) {
                iu = 1;
                ldwrku = *n;
                itau = iu + ldwrku * *n;
                nwork = itau + *n;
                i__2 = *lwork - nwork + 1;
                dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &i__2, &ierr);
                dlacpy_((char *)"L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (ftnlen)1);
                i__2 = *lwork - nwork + 1;
                dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &work[nwork], &i__2, &ierr);
                i__2 = *n - 1;
                i__1 = *n - 1;
                dlaset_((char *)"L", &i__2, &i__1, &c_b63, &c_b63, &a[a_dim1 + 2], lda, (ftnlen)1);
                ie = itau;
                itauq = ie + *n;
                itaup = itauq + *n;
                nwork = itaup + *n;
                i__2 = *lwork - nwork + 1;
                dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &work[itaup],
                        &work[nwork], &i__2, &ierr);
                dbdsdc_((char *)"U", (char *)"I", n, &s[1], &work[ie], &work[iu], n, &vt[vt_offset], ldvt, dum,
                        idum, &work[nwork], &iwork[1], info, (ftnlen)1, (ftnlen)1);
                i__2 = *lwork - nwork + 1;
                dormbr_((char *)"Q", (char *)"L", (char *)"N", n, n, n, &a[a_offset], lda, &work[itauq], &work[iu], &ldwrku,
                        &work[nwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                i__2 = *lwork - nwork + 1;
                dormbr_((char *)"P", (char *)"R", (char *)"T", n, n, n, &a[a_offset], lda, &work[itaup], &vt[vt_offset],
                        ldvt, &work[nwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                dgemm_((char *)"N", (char *)"N", m, n, n, &c_b84, &u[u_offset], ldu, &work[iu], &ldwrku, &c_b63,
                       &a[a_offset], lda, (ftnlen)1, (ftnlen)1);
                dlacpy_((char *)"F", m, n, &a[a_offset], lda, &u[u_offset], ldu, (ftnlen)1);
            }
        } else {
            ie = 1;
            itauq = ie + *n;
            itaup = itauq + *n;
            nwork = itaup + *n;
            i__2 = *lwork - nwork + 1;
            dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &work[itaup],
                    &work[nwork], &i__2, &ierr);
            if (wntqn) {
                dbdsdc_((char *)"U", (char *)"N", n, &s[1], &work[ie], dum, &c__1, dum, &c__1, dum, idum,
                        &work[nwork], &iwork[1], info, (ftnlen)1, (ftnlen)1);
            } else if (wntqo) {
                iu = nwork;
                if (*lwork >= *m * *n + *n * 3 + bdspac) {
                    ldwrku = *m;
                    nwork = iu + ldwrku * *n;
                    dlaset_((char *)"F", m, n, &c_b63, &c_b63, &work[iu], &ldwrku, (ftnlen)1);
                    ir = -1;
                } else {
                    ldwrku = *n;
                    nwork = iu + ldwrku * *n;
                    ir = nwork;
                    ldwrkr = (*lwork - *n * *n - *n * 3) / *n;
                }
                nwork = iu + ldwrku * *n;
                dbdsdc_((char *)"U", (char *)"I", n, &s[1], &work[ie], &work[iu], &ldwrku, &vt[vt_offset], ldvt,
                        dum, idum, &work[nwork], &iwork[1], info, (ftnlen)1, (ftnlen)1);
                i__2 = *lwork - nwork + 1;
                dormbr_((char *)"P", (char *)"R", (char *)"T", n, n, n, &a[a_offset], lda, &work[itaup], &vt[vt_offset],
                        ldvt, &work[nwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                if (*lwork >= *m * *n + *n * 3 + bdspac) {
                    i__2 = *lwork - nwork + 1;
                    dormbr_((char *)"Q", (char *)"L", (char *)"N", m, n, n, &a[a_offset], lda, &work[itauq], &work[iu],
                            &ldwrku, &work[nwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                    dlacpy_((char *)"F", m, n, &work[iu], &ldwrku, &a[a_offset], lda, (ftnlen)1);
                } else {
                    i__2 = *lwork - nwork + 1;
                    dorgbr_((char *)"Q", m, n, n, &a[a_offset], lda, &work[itauq], &work[nwork], &i__2,
                            &ierr, (ftnlen)1);
                    i__2 = *m;
                    i__1 = ldwrkr;
                    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
                        i__3 = *m - i__ + 1;
                        chunk = min(i__3, ldwrkr);
                        dgemm_((char *)"N", (char *)"N", &chunk, n, n, &c_b84, &a[i__ + a_dim1], lda, &work[iu],
                               &ldwrku, &c_b63, &work[ir], &ldwrkr, (ftnlen)1, (ftnlen)1);
                        dlacpy_((char *)"F", &chunk, n, &work[ir], &ldwrkr, &a[i__ + a_dim1], lda,
                                (ftnlen)1);
                    }
                }
            } else if (wntqs) {
                dlaset_((char *)"F", m, n, &c_b63, &c_b63, &u[u_offset], ldu, (ftnlen)1);
                dbdsdc_((char *)"U", (char *)"I", n, &s[1], &work[ie], &u[u_offset], ldu, &vt[vt_offset], ldvt, dum,
                        idum, &work[nwork], &iwork[1], info, (ftnlen)1, (ftnlen)1);
                i__1 = *lwork - nwork + 1;
                dormbr_((char *)"Q", (char *)"L", (char *)"N", m, n, n, &a[a_offset], lda, &work[itauq], &u[u_offset], ldu,
                        &work[nwork], &i__1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                i__1 = *lwork - nwork + 1;
                dormbr_((char *)"P", (char *)"R", (char *)"T", n, n, n, &a[a_offset], lda, &work[itaup], &vt[vt_offset],
                        ldvt, &work[nwork], &i__1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            } else if (wntqa) {
                dlaset_((char *)"F", m, m, &c_b63, &c_b63, &u[u_offset], ldu, (ftnlen)1);
                dbdsdc_((char *)"U", (char *)"I", n, &s[1], &work[ie], &u[u_offset], ldu, &vt[vt_offset], ldvt, dum,
                        idum, &work[nwork], &iwork[1], info, (ftnlen)1, (ftnlen)1);
                if (*m > *n) {
                    i__1 = *m - *n;
                    i__2 = *m - *n;
                    dlaset_((char *)"F", &i__1, &i__2, &c_b63, &c_b84, &u[*n + 1 + (*n + 1) * u_dim1], ldu,
                            (ftnlen)1);
                }
                i__1 = *lwork - nwork + 1;
                dormbr_((char *)"Q", (char *)"L", (char *)"N", m, m, n, &a[a_offset], lda, &work[itauq], &u[u_offset], ldu,
                        &work[nwork], &i__1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                i__1 = *lwork - nwork + 1;
                dormbr_((char *)"P", (char *)"R", (char *)"T", n, n, m, &a[a_offset], lda, &work[itaup], &vt[vt_offset],
                        ldvt, &work[nwork], &i__1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            }
        }
    } else {
        if (*n >= mnthr) {
            if (wntqn) {
                itau = 1;
                nwork = itau + *m;
                i__1 = *lwork - nwork + 1;
                dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &i__1, &ierr);
                i__1 = *m - 1;
                i__2 = *m - 1;
                dlaset_((char *)"U", &i__1, &i__2, &c_b63, &c_b63, &a[(a_dim1 << 1) + 1], lda, (ftnlen)1);
                ie = 1;
                itauq = ie + *m;
                itaup = itauq + *m;
                nwork = itaup + *m;
                i__1 = *lwork - nwork + 1;
                dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &work[itaup],
                        &work[nwork], &i__1, &ierr);
                nwork = ie + *m;
                dbdsdc_((char *)"U", (char *)"N", m, &s[1], &work[ie], dum, &c__1, dum, &c__1, dum, idum,
                        &work[nwork], &iwork[1], info, (ftnlen)1, (ftnlen)1);
            } else if (wntqo) {
                ivt = 1;
                il = ivt + *m * *m;
                if (*lwork >= *m * *n + *m * *m + *m * 3 + bdspac) {
                    ldwrkl = *m;
                    chunk = *n;
                } else {
                    ldwrkl = *m;
                    chunk = (*lwork - *m * *m) / *m;
                }
                itau = il + ldwrkl * *m;
                nwork = itau + *m;
                i__1 = *lwork - nwork + 1;
                dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &i__1, &ierr);
                dlacpy_((char *)"L", m, m, &a[a_offset], lda, &work[il], &ldwrkl, (ftnlen)1);
                i__1 = *m - 1;
                i__2 = *m - 1;
                dlaset_((char *)"U", &i__1, &i__2, &c_b63, &c_b63, &work[il + ldwrkl], &ldwrkl, (ftnlen)1);
                i__1 = *lwork - nwork + 1;
                dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[nwork], &i__1, &ierr);
                ie = itau;
                itauq = ie + *m;
                itaup = itauq + *m;
                nwork = itaup + *m;
                i__1 = *lwork - nwork + 1;
                dgebrd_(m, m, &work[il], &ldwrkl, &s[1], &work[ie], &work[itauq], &work[itaup],
                        &work[nwork], &i__1, &ierr);
                dbdsdc_((char *)"U", (char *)"I", m, &s[1], &work[ie], &u[u_offset], ldu, &work[ivt], m, dum, idum,
                        &work[nwork], &iwork[1], info, (ftnlen)1, (ftnlen)1);
                i__1 = *lwork - nwork + 1;
                dormbr_((char *)"Q", (char *)"L", (char *)"N", m, m, m, &work[il], &ldwrkl, &work[itauq], &u[u_offset], ldu,
                        &work[nwork], &i__1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                i__1 = *lwork - nwork + 1;
                dormbr_((char *)"P", (char *)"R", (char *)"T", m, m, m, &work[il], &ldwrkl, &work[itaup], &work[ivt], m,
                        &work[nwork], &i__1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                i__1 = *n;
                i__2 = chunk;
                for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
                    i__3 = *n - i__ + 1;
                    blk = min(i__3, chunk);
                    dgemm_((char *)"N", (char *)"N", m, &blk, m, &c_b84, &work[ivt], m, &a[i__ * a_dim1 + 1], lda,
                           &c_b63, &work[il], &ldwrkl, (ftnlen)1, (ftnlen)1);
                    dlacpy_((char *)"F", m, &blk, &work[il], &ldwrkl, &a[i__ * a_dim1 + 1], lda, (ftnlen)1);
                }
            } else if (wntqs) {
                il = 1;
                ldwrkl = *m;
                itau = il + ldwrkl * *m;
                nwork = itau + *m;
                i__2 = *lwork - nwork + 1;
                dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &i__2, &ierr);
                dlacpy_((char *)"L", m, m, &a[a_offset], lda, &work[il], &ldwrkl, (ftnlen)1);
                i__2 = *m - 1;
                i__1 = *m - 1;
                dlaset_((char *)"U", &i__2, &i__1, &c_b63, &c_b63, &work[il + ldwrkl], &ldwrkl, (ftnlen)1);
                i__2 = *lwork - nwork + 1;
                dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[nwork], &i__2, &ierr);
                ie = itau;
                itauq = ie + *m;
                itaup = itauq + *m;
                nwork = itaup + *m;
                i__2 = *lwork - nwork + 1;
                dgebrd_(m, m, &work[il], &ldwrkl, &s[1], &work[ie], &work[itauq], &work[itaup],
                        &work[nwork], &i__2, &ierr);
                dbdsdc_((char *)"U", (char *)"I", m, &s[1], &work[ie], &u[u_offset], ldu, &vt[vt_offset], ldvt, dum,
                        idum, &work[nwork], &iwork[1], info, (ftnlen)1, (ftnlen)1);
                i__2 = *lwork - nwork + 1;
                dormbr_((char *)"Q", (char *)"L", (char *)"N", m, m, m, &work[il], &ldwrkl, &work[itauq], &u[u_offset], ldu,
                        &work[nwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                i__2 = *lwork - nwork + 1;
                dormbr_((char *)"P", (char *)"R", (char *)"T", m, m, m, &work[il], &ldwrkl, &work[itaup], &vt[vt_offset],
                        ldvt, &work[nwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                dlacpy_((char *)"F", m, m, &vt[vt_offset], ldvt, &work[il], &ldwrkl, (ftnlen)1);
                dgemm_((char *)"N", (char *)"N", m, n, m, &c_b84, &work[il], &ldwrkl, &a[a_offset], lda, &c_b63,
                       &vt[vt_offset], ldvt, (ftnlen)1, (ftnlen)1);
            } else if (wntqa) {
                ivt = 1;
                ldwkvt = *m;
                itau = ivt + ldwkvt * *m;
                nwork = itau + *m;
                i__2 = *lwork - nwork + 1;
                dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &i__2, &ierr);
                dlacpy_((char *)"U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (ftnlen)1);
                i__2 = *lwork - nwork + 1;
                dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &work[nwork], &i__2, &ierr);
                i__2 = *m - 1;
                i__1 = *m - 1;
                dlaset_((char *)"U", &i__2, &i__1, &c_b63, &c_b63, &a[(a_dim1 << 1) + 1], lda, (ftnlen)1);
                ie = itau;
                itauq = ie + *m;
                itaup = itauq + *m;
                nwork = itaup + *m;
                i__2 = *lwork - nwork + 1;
                dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &work[itaup],
                        &work[nwork], &i__2, &ierr);
                dbdsdc_((char *)"U", (char *)"I", m, &s[1], &work[ie], &u[u_offset], ldu, &work[ivt], &ldwkvt, dum,
                        idum, &work[nwork], &iwork[1], info, (ftnlen)1, (ftnlen)1);
                i__2 = *lwork - nwork + 1;
                dormbr_((char *)"Q", (char *)"L", (char *)"N", m, m, m, &a[a_offset], lda, &work[itauq], &u[u_offset], ldu,
                        &work[nwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                i__2 = *lwork - nwork + 1;
                dormbr_((char *)"P", (char *)"R", (char *)"T", m, m, m, &a[a_offset], lda, &work[itaup], &work[ivt],
                        &ldwkvt, &work[nwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                dgemm_((char *)"N", (char *)"N", m, n, m, &c_b84, &work[ivt], &ldwkvt, &vt[vt_offset], ldvt, &c_b63,
                       &a[a_offset], lda, (ftnlen)1, (ftnlen)1);
                dlacpy_((char *)"F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (ftnlen)1);
            }
        } else {
            ie = 1;
            itauq = ie + *m;
            itaup = itauq + *m;
            nwork = itaup + *m;
            i__2 = *lwork - nwork + 1;
            dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &work[itaup],
                    &work[nwork], &i__2, &ierr);
            if (wntqn) {
                dbdsdc_((char *)"L", (char *)"N", m, &s[1], &work[ie], dum, &c__1, dum, &c__1, dum, idum,
                        &work[nwork], &iwork[1], info, (ftnlen)1, (ftnlen)1);
            } else if (wntqo) {
                ldwkvt = *m;
                ivt = nwork;
                if (*lwork >= *m * *n + *m * 3 + bdspac) {
                    dlaset_((char *)"F", m, n, &c_b63, &c_b63, &work[ivt], &ldwkvt, (ftnlen)1);
                    nwork = ivt + ldwkvt * *n;
                    il = -1;
                } else {
                    nwork = ivt + ldwkvt * *m;
                    il = nwork;
                    chunk = (*lwork - *m * *m - *m * 3) / *m;
                }
                dbdsdc_((char *)"L", (char *)"I", m, &s[1], &work[ie], &u[u_offset], ldu, &work[ivt], &ldwkvt, dum,
                        idum, &work[nwork], &iwork[1], info, (ftnlen)1, (ftnlen)1);
                i__2 = *lwork - nwork + 1;
                dormbr_((char *)"Q", (char *)"L", (char *)"N", m, m, n, &a[a_offset], lda, &work[itauq], &u[u_offset], ldu,
                        &work[nwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                if (*lwork >= *m * *n + *m * 3 + bdspac) {
                    i__2 = *lwork - nwork + 1;
                    dormbr_((char *)"P", (char *)"R", (char *)"T", m, n, m, &a[a_offset], lda, &work[itaup], &work[ivt],
                            &ldwkvt, &work[nwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                    dlacpy_((char *)"F", m, n, &work[ivt], &ldwkvt, &a[a_offset], lda, (ftnlen)1);
                } else {
                    i__2 = *lwork - nwork + 1;
                    dorgbr_((char *)"P", m, n, m, &a[a_offset], lda, &work[itaup], &work[nwork], &i__2,
                            &ierr, (ftnlen)1);
                    i__2 = *n;
                    i__1 = chunk;
                    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
                        i__3 = *n - i__ + 1;
                        blk = min(i__3, chunk);
                        dgemm_((char *)"N", (char *)"N", m, &blk, m, &c_b84, &work[ivt], &ldwkvt,
                               &a[i__ * a_dim1 + 1], lda, &c_b63, &work[il], m, (ftnlen)1,
                               (ftnlen)1);
                        dlacpy_((char *)"F", m, &blk, &work[il], m, &a[i__ * a_dim1 + 1], lda, (ftnlen)1);
                    }
                }
            } else if (wntqs) {
                dlaset_((char *)"F", m, n, &c_b63, &c_b63, &vt[vt_offset], ldvt, (ftnlen)1);
                dbdsdc_((char *)"L", (char *)"I", m, &s[1], &work[ie], &u[u_offset], ldu, &vt[vt_offset], ldvt, dum,
                        idum, &work[nwork], &iwork[1], info, (ftnlen)1, (ftnlen)1);
                i__1 = *lwork - nwork + 1;
                dormbr_((char *)"Q", (char *)"L", (char *)"N", m, m, n, &a[a_offset], lda, &work[itauq], &u[u_offset], ldu,
                        &work[nwork], &i__1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                i__1 = *lwork - nwork + 1;
                dormbr_((char *)"P", (char *)"R", (char *)"T", m, n, m, &a[a_offset], lda, &work[itaup], &vt[vt_offset],
                        ldvt, &work[nwork], &i__1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            } else if (wntqa) {
                dlaset_((char *)"F", n, n, &c_b63, &c_b63, &vt[vt_offset], ldvt, (ftnlen)1);
                dbdsdc_((char *)"L", (char *)"I", m, &s[1], &work[ie], &u[u_offset], ldu, &vt[vt_offset], ldvt, dum,
                        idum, &work[nwork], &iwork[1], info, (ftnlen)1, (ftnlen)1);
                if (*n > *m) {
                    i__1 = *n - *m;
                    i__2 = *n - *m;
                    dlaset_((char *)"F", &i__1, &i__2, &c_b63, &c_b84, &vt[*m + 1 + (*m + 1) * vt_dim1],
                            ldvt, (ftnlen)1);
                }
                i__1 = *lwork - nwork + 1;
                dormbr_((char *)"Q", (char *)"L", (char *)"N", m, m, n, &a[a_offset], lda, &work[itauq], &u[u_offset], ldu,
                        &work[nwork], &i__1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                i__1 = *lwork - nwork + 1;
                dormbr_((char *)"P", (char *)"R", (char *)"T", n, n, m, &a[a_offset], lda, &work[itaup], &vt[vt_offset],
                        ldvt, &work[nwork], &i__1, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            }
        }
    }
    if (iscl == 1) {
        if (anrm > bignum) {
            dlascl_((char *)"G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &minmn, &ierr,
                    (ftnlen)1);
        }
        if (anrm < smlnum) {
            dlascl_((char *)"G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &minmn, &ierr,
                    (ftnlen)1);
        }
    }
    work[1] = droundup_lwork__(&maxwrk);
    return 0;
}
#ifdef __cplusplus
}
#endif
