#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__6 = 6;
static integer c__0 = 0;
static integer c__2 = 2;
static integer c_n1 = -1;
static doublereal c_b57 = 0.;
static integer c__1 = 1;
static doublereal c_b79 = 1.;
int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n, doublereal *a, integer *lda,
            doublereal *s, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt,
            doublereal *work, integer *lwork, integer *info, ftnlen jobu_len, ftnlen jobvt_len)
{
    address a__1[2];
    integer a_dim1, a_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1[2], i__2, i__3, i__4;
    char ch__1[2];
    int s_lmp_cat(char *, char **, integer *, integer *, ftnlen);
    double sqrt(doublereal);
    integer i__, ie, ir, iu, blk, ncu;
    doublereal dum[1], eps;
    integer nru, iscl;
    doublereal anrm;
    integer ierr, itau, ncvt, nrvt, lwork_dgebrd__, lwork_dgelqf__, lwork_dgeqrf__;
    extern int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
                      integer *, doublereal *, integer *, doublereal *, doublereal *, integer *,
                      ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer chunk, minmn, wrkbl, itaup, itauq, mnthr, iwork;
    logical wntua, wntva, wntun, wntuo, wntvn, wntvo, wntus, wntvs;
    extern int dgebrd_(integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
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
        dbdsqr_(char *, integer *, integer *, integer *, integer *, doublereal *, doublereal *,
                doublereal *, integer *, doublereal *, integer *, doublereal *, integer *,
                doublereal *, integer *, ftnlen),
        dorgbr_(char *, integer *, integer *, integer *, doublereal *, integer *, doublereal *,
                doublereal *, integer *, integer *, ftnlen);
    doublereal bignum;
    extern int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    extern int dormbr_(char *, char *, char *, integer *, integer *, integer *, doublereal *,
                       integer *, doublereal *, doublereal *, integer *, doublereal *, integer *,
                       integer *, ftnlen, ftnlen, ftnlen),
        dorglq_(integer *, integer *, integer *, doublereal *, integer *, doublereal *,
                doublereal *, integer *, integer *),
        dorgqr_(integer *, integer *, integer *, doublereal *, integer *, doublereal *,
                doublereal *, integer *, integer *);
    integer ldwrkr, minwrk, ldwrku, maxwrk;
    doublereal smlnum;
    logical lquery, wntuas, wntvas;
    integer lwork_dorgbr_p__, lwork_dorgbr_q__, lwork_dorglq_m__, lwork_dorglq_n__,
        lwork_dorgqr_m__, lwork_dorgqr_n__;
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
    *info = 0;
    minmn = min(*m, *n);
    wntua = lsame_(jobu, (char *)"A", (ftnlen)1, (ftnlen)1);
    wntus = lsame_(jobu, (char *)"S", (ftnlen)1, (ftnlen)1);
    wntuas = wntua || wntus;
    wntuo = lsame_(jobu, (char *)"O", (ftnlen)1, (ftnlen)1);
    wntun = lsame_(jobu, (char *)"N", (ftnlen)1, (ftnlen)1);
    wntva = lsame_(jobvt, (char *)"A", (ftnlen)1, (ftnlen)1);
    wntvs = lsame_(jobvt, (char *)"S", (ftnlen)1, (ftnlen)1);
    wntvas = wntva || wntvs;
    wntvo = lsame_(jobvt, (char *)"O", (ftnlen)1, (ftnlen)1);
    wntvn = lsame_(jobvt, (char *)"N", (ftnlen)1, (ftnlen)1);
    lquery = *lwork == -1;
    if (!(wntua || wntus || wntuo || wntun)) {
        *info = -1;
    } else if (!(wntva || wntvs || wntvo || wntvn) || wntvo && wntuo) {
        *info = -2;
    } else if (*m < 0) {
        *info = -3;
    } else if (*n < 0) {
        *info = -4;
    } else if (*lda < max(1, *m)) {
        *info = -6;
    } else if (*ldu < 1 || wntuas && *ldu < *m) {
        *info = -9;
    } else if (*ldvt < 1 || wntva && *ldvt < *n || wntvs && *ldvt < minmn) {
        *info = -11;
    }
    if (*info == 0) {
        minwrk = 1;
        maxwrk = 1;
        if (*m >= *n && minmn > 0) {
            i__1[0] = 1, a__1[0] = jobu;
            i__1[1] = 1, a__1[1] = jobvt;
            s_lmp_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
            mnthr = ilaenv_(&c__6, (char *)"DGESVD", ch__1, m, n, &c__0, &c__0, (ftnlen)6, (ftnlen)2);
            bdspac = *n * 5;
            dgeqrf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
            lwork_dgeqrf__ = (integer)dum[0];
            dorgqr_(m, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
            lwork_dorgqr_n__ = (integer)dum[0];
            dorgqr_(m, m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
            lwork_dorgqr_m__ = (integer)dum[0];
            dgebrd_(n, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &c_n1, &ierr);
            lwork_dgebrd__ = (integer)dum[0];
            dorgbr_((char *)"P", n, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr, (ftnlen)1);
            lwork_dorgbr_p__ = (integer)dum[0];
            dorgbr_((char *)"Q", n, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr, (ftnlen)1);
            lwork_dorgbr_q__ = (integer)dum[0];
            if (*m >= mnthr) {
                if (wntun) {
                    maxwrk = *n + lwork_dgeqrf__;
                    i__2 = maxwrk, i__3 = *n * 3 + lwork_dgebrd__;
                    maxwrk = max(i__2, i__3);
                    if (wntvo || wntvas) {
                        i__2 = maxwrk, i__3 = *n * 3 + lwork_dorgbr_p__;
                        maxwrk = max(i__2, i__3);
                    }
                    maxwrk = max(maxwrk, bdspac);
                    i__2 = *n << 2;
                    minwrk = max(i__2, bdspac);
                } else if (wntuo && wntvn) {
                    wrkbl = *n + lwork_dgeqrf__;
                    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_n__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2, i__3);
                    wrkbl = max(wrkbl, bdspac);
                    i__2 = *n * *n + wrkbl, i__3 = *n * *n + *m * *n + *n;
                    maxwrk = max(i__2, i__3);
                    i__2 = *n * 3 + *m;
                    minwrk = max(i__2, bdspac);
                } else if (wntuo && wntvas) {
                    wrkbl = *n + lwork_dgeqrf__;
                    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_n__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2, i__3);
                    wrkbl = max(wrkbl, bdspac);
                    i__2 = *n * *n + wrkbl, i__3 = *n * *n + *m * *n + *n;
                    maxwrk = max(i__2, i__3);
                    i__2 = *n * 3 + *m;
                    minwrk = max(i__2, bdspac);
                } else if (wntus && wntvn) {
                    wrkbl = *n + lwork_dgeqrf__;
                    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_n__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2, i__3);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = *n * *n + wrkbl;
                    i__2 = *n * 3 + *m;
                    minwrk = max(i__2, bdspac);
                } else if (wntus && wntvo) {
                    wrkbl = *n + lwork_dgeqrf__;
                    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_n__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2, i__3);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = (*n << 1) * *n + wrkbl;
                    i__2 = *n * 3 + *m;
                    minwrk = max(i__2, bdspac);
                } else if (wntus && wntvas) {
                    wrkbl = *n + lwork_dgeqrf__;
                    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_n__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2, i__3);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = *n * *n + wrkbl;
                    i__2 = *n * 3 + *m;
                    minwrk = max(i__2, bdspac);
                } else if (wntua && wntvn) {
                    wrkbl = *n + lwork_dgeqrf__;
                    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_m__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2, i__3);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = *n * *n + wrkbl;
                    i__2 = *n * 3 + *m;
                    minwrk = max(i__2, bdspac);
                } else if (wntua && wntvo) {
                    wrkbl = *n + lwork_dgeqrf__;
                    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_m__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2, i__3);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = (*n << 1) * *n + wrkbl;
                    i__2 = *n * 3 + *m;
                    minwrk = max(i__2, bdspac);
                } else if (wntua && wntvas) {
                    wrkbl = *n + lwork_dgeqrf__;
                    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_m__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2, i__3);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = *n * *n + wrkbl;
                    i__2 = *n * 3 + *m;
                    minwrk = max(i__2, bdspac);
                }
            } else {
                dgebrd_(m, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &c_n1, &ierr);
                lwork_dgebrd__ = (integer)dum[0];
                maxwrk = *n * 3 + lwork_dgebrd__;
                if (wntus || wntuo) {
                    dorgbr_((char *)"Q", m, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr, (ftnlen)1);
                    lwork_dorgbr_q__ = (integer)dum[0];
                    i__2 = maxwrk, i__3 = *n * 3 + lwork_dorgbr_q__;
                    maxwrk = max(i__2, i__3);
                }
                if (wntua) {
                    dorgbr_((char *)"Q", m, m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr, (ftnlen)1);
                    lwork_dorgbr_q__ = (integer)dum[0];
                    i__2 = maxwrk, i__3 = *n * 3 + lwork_dorgbr_q__;
                    maxwrk = max(i__2, i__3);
                }
                if (!wntvn) {
                    i__2 = maxwrk, i__3 = *n * 3 + lwork_dorgbr_p__;
                    maxwrk = max(i__2, i__3);
                }
                maxwrk = max(maxwrk, bdspac);
                i__2 = *n * 3 + *m;
                minwrk = max(i__2, bdspac);
            }
        } else if (minmn > 0) {
            i__1[0] = 1, a__1[0] = jobu;
            i__1[1] = 1, a__1[1] = jobvt;
            s_lmp_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
            mnthr = ilaenv_(&c__6, (char *)"DGESVD", ch__1, m, n, &c__0, &c__0, (ftnlen)6, (ftnlen)2);
            bdspac = *m * 5;
            dgelqf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
            lwork_dgelqf__ = (integer)dum[0];
            dorglq_(n, n, m, dum, n, dum, dum, &c_n1, &ierr);
            lwork_dorglq_n__ = (integer)dum[0];
            dorglq_(m, n, m, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
            lwork_dorglq_m__ = (integer)dum[0];
            dgebrd_(m, m, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &c_n1, &ierr);
            lwork_dgebrd__ = (integer)dum[0];
            dorgbr_((char *)"P", m, m, m, &a[a_offset], n, dum, dum, &c_n1, &ierr, (ftnlen)1);
            lwork_dorgbr_p__ = (integer)dum[0];
            dorgbr_((char *)"Q", m, m, m, &a[a_offset], n, dum, dum, &c_n1, &ierr, (ftnlen)1);
            lwork_dorgbr_q__ = (integer)dum[0];
            if (*n >= mnthr) {
                if (wntvn) {
                    maxwrk = *m + lwork_dgelqf__;
                    i__2 = maxwrk, i__3 = *m * 3 + lwork_dgebrd__;
                    maxwrk = max(i__2, i__3);
                    if (wntuo || wntuas) {
                        i__2 = maxwrk, i__3 = *m * 3 + lwork_dorgbr_q__;
                        maxwrk = max(i__2, i__3);
                    }
                    maxwrk = max(maxwrk, bdspac);
                    i__2 = *m << 2;
                    minwrk = max(i__2, bdspac);
                } else if (wntvo && wntun) {
                    wrkbl = *m + lwork_dgelqf__;
                    i__2 = wrkbl, i__3 = *m + lwork_dorglq_m__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2, i__3);
                    wrkbl = max(wrkbl, bdspac);
                    i__2 = *m * *m + wrkbl, i__3 = *m * *m + *m * *n + *m;
                    maxwrk = max(i__2, i__3);
                    i__2 = *m * 3 + *n;
                    minwrk = max(i__2, bdspac);
                } else if (wntvo && wntuas) {
                    wrkbl = *m + lwork_dgelqf__;
                    i__2 = wrkbl, i__3 = *m + lwork_dorglq_m__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2, i__3);
                    wrkbl = max(wrkbl, bdspac);
                    i__2 = *m * *m + wrkbl, i__3 = *m * *m + *m * *n + *m;
                    maxwrk = max(i__2, i__3);
                    i__2 = *m * 3 + *n;
                    minwrk = max(i__2, bdspac);
                } else if (wntvs && wntun) {
                    wrkbl = *m + lwork_dgelqf__;
                    i__2 = wrkbl, i__3 = *m + lwork_dorglq_m__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2, i__3);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = *m * *m + wrkbl;
                    i__2 = *m * 3 + *n;
                    minwrk = max(i__2, bdspac);
                } else if (wntvs && wntuo) {
                    wrkbl = *m + lwork_dgelqf__;
                    i__2 = wrkbl, i__3 = *m + lwork_dorglq_m__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2, i__3);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = (*m << 1) * *m + wrkbl;
                    i__2 = *m * 3 + *n;
                    minwrk = max(i__2, bdspac);
                } else if (wntvs && wntuas) {
                    wrkbl = *m + lwork_dgelqf__;
                    i__2 = wrkbl, i__3 = *m + lwork_dorglq_m__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2, i__3);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = *m * *m + wrkbl;
                    i__2 = *m * 3 + *n;
                    minwrk = max(i__2, bdspac);
                } else if (wntva && wntun) {
                    wrkbl = *m + lwork_dgelqf__;
                    i__2 = wrkbl, i__3 = *m + lwork_dorglq_n__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2, i__3);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = *m * *m + wrkbl;
                    i__2 = *m * 3 + *n;
                    minwrk = max(i__2, bdspac);
                } else if (wntva && wntuo) {
                    wrkbl = *m + lwork_dgelqf__;
                    i__2 = wrkbl, i__3 = *m + lwork_dorglq_n__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2, i__3);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = (*m << 1) * *m + wrkbl;
                    i__2 = *m * 3 + *n;
                    minwrk = max(i__2, bdspac);
                } else if (wntva && wntuas) {
                    wrkbl = *m + lwork_dgelqf__;
                    i__2 = wrkbl, i__3 = *m + lwork_dorglq_n__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2, i__3);
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2, i__3);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = *m * *m + wrkbl;
                    i__2 = *m * 3 + *n;
                    minwrk = max(i__2, bdspac);
                }
            } else {
                dgebrd_(m, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &c_n1, &ierr);
                lwork_dgebrd__ = (integer)dum[0];
                maxwrk = *m * 3 + lwork_dgebrd__;
                if (wntvs || wntvo) {
                    dorgbr_((char *)"P", m, n, m, &a[a_offset], n, dum, dum, &c_n1, &ierr, (ftnlen)1);
                    lwork_dorgbr_p__ = (integer)dum[0];
                    i__2 = maxwrk, i__3 = *m * 3 + lwork_dorgbr_p__;
                    maxwrk = max(i__2, i__3);
                }
                if (wntva) {
                    dorgbr_((char *)"P", n, n, m, &a[a_offset], n, dum, dum, &c_n1, &ierr, (ftnlen)1);
                    lwork_dorgbr_p__ = (integer)dum[0];
                    i__2 = maxwrk, i__3 = *m * 3 + lwork_dorgbr_p__;
                    maxwrk = max(i__2, i__3);
                }
                if (!wntun) {
                    i__2 = maxwrk, i__3 = *m * 3 + lwork_dorgbr_q__;
                    maxwrk = max(i__2, i__3);
                }
                maxwrk = max(maxwrk, bdspac);
                i__2 = *m * 3 + *n;
                minwrk = max(i__2, bdspac);
            }
        }
        maxwrk = max(maxwrk, minwrk);
        work[1] = (doublereal)maxwrk;
        if (*lwork < minwrk && !lquery) {
            *info = -13;
        }
    }
    if (*info != 0) {
        i__2 = -(*info);
        xerbla_((char *)"DGESVD", &i__2, (ftnlen)6);
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
            if (wntun) {
                itau = 1;
                iwork = itau + *n;
                i__2 = *lwork - iwork + 1;
                dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                if (*n > 1) {
                    i__2 = *n - 1;
                    i__3 = *n - 1;
                    dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &a[a_dim1 + 2], lda, (ftnlen)1);
                }
                ie = 1;
                itauq = ie + *n;
                itaup = itauq + *n;
                iwork = itaup + *n;
                i__2 = *lwork - iwork + 1;
                dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &work[itaup],
                        &work[iwork], &i__2, &ierr);
                ncvt = 0;
                if (wntvo || wntvas) {
                    i__2 = *lwork - iwork + 1;
                    dorgbr_((char *)"P", n, n, n, &a[a_offset], lda, &work[itaup], &work[iwork], &i__2,
                            &ierr, (ftnlen)1);
                    ncvt = *n;
                }
                iwork = ie + *n;
                dbdsqr_((char *)"U", n, &ncvt, &c__0, &c__0, &s[1], &work[ie], &a[a_offset], lda, dum,
                        &c__1, dum, &c__1, &work[iwork], info, (ftnlen)1);
                if (wntvas) {
                    dlacpy_((char *)"F", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (ftnlen)1);
                }
            } else if (wntuo && wntvn) {
                i__2 = *n << 2;
                if (*lwork >= *n * *n + max(i__2, bdspac)) {
                    ir = 1;
                    i__2 = wrkbl, i__3 = *lda * *n + *n;
                    if (*lwork >= max(i__2, i__3) + *lda * *n) {
                        ldwrku = *lda;
                        ldwrkr = *lda;
                    } else {
                        i__2 = wrkbl, i__3 = *lda * *n + *n;
                        if (*lwork >= max(i__2, i__3) + *n * *n) {
                            ldwrku = *lda;
                            ldwrkr = *n;
                        } else {
                            ldwrku = (*lwork - *n * *n - *n) / *n;
                            ldwrkr = *n;
                        }
                    }
                    itau = ir + ldwrkr * *n;
                    iwork = itau + *n;
                    i__2 = *lwork - iwork + 1;
                    dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                    dlacpy_((char *)"U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, (ftnlen)1);
                    i__2 = *n - 1;
                    i__3 = *n - 1;
                    dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 1], &ldwrkr, (ftnlen)1);
                    i__2 = *lwork - iwork + 1;
                    dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                    ie = itau;
                    itauq = ie + *n;
                    itaup = itauq + *n;
                    iwork = itaup + *n;
                    i__2 = *lwork - iwork + 1;
                    dgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &work[itauq], &work[itaup],
                            &work[iwork], &i__2, &ierr);
                    i__2 = *lwork - iwork + 1;
                    dorgbr_((char *)"Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &work[iwork], &i__2,
                            &ierr, (ftnlen)1);
                    iwork = ie + *n;
                    dbdsqr_((char *)"U", n, &c__0, n, &c__0, &s[1], &work[ie], dum, &c__1, &work[ir],
                            &ldwrkr, dum, &c__1, &work[iwork], info, (ftnlen)1);
                    iu = ie + *n;
                    i__2 = *m;
                    i__3 = ldwrku;
                    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) {
                        i__4 = *m - i__ + 1;
                        chunk = min(i__4, ldwrku);
                        dgemm_((char *)"N", (char *)"N", &chunk, n, n, &c_b79, &a[i__ + a_dim1], lda, &work[ir],
                               &ldwrkr, &c_b57, &work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
                        dlacpy_((char *)"F", &chunk, n, &work[iu], &ldwrku, &a[i__ + a_dim1], lda,
                                (ftnlen)1);
                    }
                } else {
                    ie = 1;
                    itauq = ie + *n;
                    itaup = itauq + *n;
                    iwork = itaup + *n;
                    i__3 = *lwork - iwork + 1;
                    dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &work[itaup],
                            &work[iwork], &i__3, &ierr);
                    i__3 = *lwork - iwork + 1;
                    dorgbr_((char *)"Q", m, n, n, &a[a_offset], lda, &work[itauq], &work[iwork], &i__3,
                            &ierr, (ftnlen)1);
                    iwork = ie + *n;
                    dbdsqr_((char *)"U", n, &c__0, m, &c__0, &s[1], &work[ie], dum, &c__1, &a[a_offset],
                            lda, dum, &c__1, &work[iwork], info, (ftnlen)1);
                }
            } else if (wntuo && wntvas) {
                i__3 = *n << 2;
                if (*lwork >= *n * *n + max(i__3, bdspac)) {
                    ir = 1;
                    i__3 = wrkbl, i__2 = *lda * *n + *n;
                    if (*lwork >= max(i__3, i__2) + *lda * *n) {
                        ldwrku = *lda;
                        ldwrkr = *lda;
                    } else {
                        i__3 = wrkbl, i__2 = *lda * *n + *n;
                        if (*lwork >= max(i__3, i__2) + *n * *n) {
                            ldwrku = *lda;
                            ldwrkr = *n;
                        } else {
                            ldwrku = (*lwork - *n * *n - *n) / *n;
                            ldwrkr = *n;
                        }
                    }
                    itau = ir + ldwrkr * *n;
                    iwork = itau + *n;
                    i__3 = *lwork - iwork + 1;
                    dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__3, &ierr);
                    dlacpy_((char *)"U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (ftnlen)1);
                    if (*n > 1) {
                        i__3 = *n - 1;
                        i__2 = *n - 1;
                        dlaset_((char *)"L", &i__3, &i__2, &c_b57, &c_b57, &vt[vt_dim1 + 2], ldvt,
                                (ftnlen)1);
                    }
                    i__3 = *lwork - iwork + 1;
                    dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__3, &ierr);
                    ie = itau;
                    itauq = ie + *n;
                    itaup = itauq + *n;
                    iwork = itaup + *n;
                    i__3 = *lwork - iwork + 1;
                    dgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie], &work[itauq],
                            &work[itaup], &work[iwork], &i__3, &ierr);
                    dlacpy_((char *)"L", n, n, &vt[vt_offset], ldvt, &work[ir], &ldwrkr, (ftnlen)1);
                    i__3 = *lwork - iwork + 1;
                    dorgbr_((char *)"Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &work[iwork], &i__3,
                            &ierr, (ftnlen)1);
                    i__3 = *lwork - iwork + 1;
                    dorgbr_((char *)"P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &work[iwork], &i__3,
                            &ierr, (ftnlen)1);
                    iwork = ie + *n;
                    dbdsqr_((char *)"U", n, n, n, &c__0, &s[1], &work[ie], &vt[vt_offset], ldvt, &work[ir],
                            &ldwrkr, dum, &c__1, &work[iwork], info, (ftnlen)1);
                    iu = ie + *n;
                    i__3 = *m;
                    i__2 = ldwrku;
                    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ += i__2) {
                        i__4 = *m - i__ + 1;
                        chunk = min(i__4, ldwrku);
                        dgemm_((char *)"N", (char *)"N", &chunk, n, n, &c_b79, &a[i__ + a_dim1], lda, &work[ir],
                               &ldwrkr, &c_b57, &work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
                        dlacpy_((char *)"F", &chunk, n, &work[iu], &ldwrku, &a[i__ + a_dim1], lda,
                                (ftnlen)1);
                    }
                } else {
                    itau = 1;
                    iwork = itau + *n;
                    i__2 = *lwork - iwork + 1;
                    dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                    dlacpy_((char *)"U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (ftnlen)1);
                    if (*n > 1) {
                        i__2 = *n - 1;
                        i__3 = *n - 1;
                        dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &vt[vt_dim1 + 2], ldvt,
                                (ftnlen)1);
                    }
                    i__2 = *lwork - iwork + 1;
                    dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                    ie = itau;
                    itauq = ie + *n;
                    itaup = itauq + *n;
                    iwork = itaup + *n;
                    i__2 = *lwork - iwork + 1;
                    dgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie], &work[itauq],
                            &work[itaup], &work[iwork], &i__2, &ierr);
                    i__2 = *lwork - iwork + 1;
                    dormbr_((char *)"Q", (char *)"R", (char *)"N", m, n, n, &vt[vt_offset], ldvt, &work[itauq],
                            &a[a_offset], lda, &work[iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1,
                            (ftnlen)1);
                    i__2 = *lwork - iwork + 1;
                    dorgbr_((char *)"P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &work[iwork], &i__2,
                            &ierr, (ftnlen)1);
                    iwork = ie + *n;
                    dbdsqr_((char *)"U", n, n, m, &c__0, &s[1], &work[ie], &vt[vt_offset], ldvt,
                            &a[a_offset], lda, dum, &c__1, &work[iwork], info, (ftnlen)1);
                }
            } else if (wntus) {
                if (wntvn) {
                    i__2 = *n << 2;
                    if (*lwork >= *n * *n + max(i__2, bdspac)) {
                        ir = 1;
                        if (*lwork >= wrkbl + *lda * *n) {
                            ldwrkr = *lda;
                        } else {
                            ldwrkr = *n;
                        }
                        itau = ir + ldwrkr * *n;
                        iwork = itau + *n;
                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, (ftnlen)1);
                        i__2 = *n - 1;
                        i__3 = *n - 1;
                        dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 1], &ldwrkr,
                                (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        iwork = ie + *n;
                        dbdsqr_((char *)"U", n, &c__0, n, &c__0, &s[1], &work[ie], dum, &c__1, &work[ir],
                                &ldwrkr, dum, &c__1, &work[iwork], info, (ftnlen)1);
                        dgemm_((char *)"N", (char *)"N", m, n, n, &c_b79, &a[a_offset], lda, &work[ir], &ldwrkr,
                               &c_b57, &u[u_offset], ldu, (ftnlen)1, (ftnlen)1);
                    } else {
                        itau = 1;
                        iwork = itau + *n;
                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, n, n, &u[u_offset], ldu, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;
                        if (*n > 1) {
                            i__2 = *n - 1;
                            i__3 = *n - 1;
                            dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &a[a_dim1 + 2], lda,
                                    (ftnlen)1);
                        }
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"Q", (char *)"R", (char *)"N", m, n, n, &a[a_offset], lda, &work[itauq],
                                &u[u_offset], ldu, &work[iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1,
                                (ftnlen)1);
                        iwork = ie + *n;
                        dbdsqr_((char *)"U", n, &c__0, m, &c__0, &s[1], &work[ie], dum, &c__1, &u[u_offset],
                                ldu, dum, &c__1, &work[iwork], info, (ftnlen)1);
                    }
                } else if (wntvo) {
                    i__2 = *n << 2;
                    if (*lwork >= (*n << 1) * *n + max(i__2, bdspac)) {
                        iu = 1;
                        if (*lwork >= wrkbl + (*lda << 1) * *n) {
                            ldwrku = *lda;
                            ir = iu + ldwrku * *n;
                            ldwrkr = *lda;
                        } else if (*lwork >= wrkbl + (*lda + *n) * *n) {
                            ldwrku = *lda;
                            ir = iu + ldwrku * *n;
                            ldwrkr = *n;
                        } else {
                            ldwrku = *n;
                            ir = iu + ldwrku * *n;
                            ldwrkr = *n;
                        }
                        itau = ir + ldwrkr * *n;
                        iwork = itau + *n;
                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", n, n, &a[a_offset], lda, &work[iu], &ldwrku, (ftnlen)1);
                        i__2 = *n - 1;
                        i__3 = *n - 1;
                        dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 1], &ldwrku,
                                (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", n, n, &work[iu], &ldwrku, &work[ir], &ldwrkr, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", n, n, n, &work[iu], &ldwrku, &work[itauq], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", n, n, n, &work[ir], &ldwrkr, &work[itaup], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        iwork = ie + *n;
                        dbdsqr_((char *)"U", n, n, n, &c__0, &s[1], &work[ie], &work[ir], &ldwrkr,
                                &work[iu], &ldwrku, dum, &c__1, &work[iwork], info, (ftnlen)1);
                        dgemm_((char *)"N", (char *)"N", m, n, n, &c_b79, &a[a_offset], lda, &work[iu], &ldwrku,
                               &c_b57, &u[u_offset], ldu, (ftnlen)1, (ftnlen)1);
                        dlacpy_((char *)"F", n, n, &work[ir], &ldwrkr, &a[a_offset], lda, (ftnlen)1);
                    } else {
                        itau = 1;
                        iwork = itau + *n;
                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, n, n, &u[u_offset], ldu, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;
                        if (*n > 1) {
                            i__2 = *n - 1;
                            i__3 = *n - 1;
                            dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &a[a_dim1 + 2], lda,
                                    (ftnlen)1);
                        }
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"Q", (char *)"R", (char *)"N", m, n, n, &a[a_offset], lda, &work[itauq],
                                &u[u_offset], ldu, &work[iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1,
                                (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", n, n, n, &a[a_offset], lda, &work[itaup], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        iwork = ie + *n;
                        dbdsqr_((char *)"U", n, n, m, &c__0, &s[1], &work[ie], &a[a_offset], lda,
                                &u[u_offset], ldu, dum, &c__1, &work[iwork], info, (ftnlen)1);
                    }
                } else if (wntvas) {
                    i__2 = *n << 2;
                    if (*lwork >= *n * *n + max(i__2, bdspac)) {
                        iu = 1;
                        if (*lwork >= wrkbl + *lda * *n) {
                            ldwrku = *lda;
                        } else {
                            ldwrku = *n;
                        }
                        itau = iu + ldwrku * *n;
                        iwork = itau + *n;
                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", n, n, &a[a_offset], lda, &work[iu], &ldwrku, (ftnlen)1);
                        i__2 = *n - 1;
                        i__3 = *n - 1;
                        dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 1], &ldwrku,
                                (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", n, n, &work[iu], &ldwrku, &vt[vt_offset], ldvt, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", n, n, n, &work[iu], &ldwrku, &work[itauq], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &work[iwork],
                                &i__2, &ierr, (ftnlen)1);
                        iwork = ie + *n;
                        dbdsqr_((char *)"U", n, n, n, &c__0, &s[1], &work[ie], &vt[vt_offset], ldvt,
                                &work[iu], &ldwrku, dum, &c__1, &work[iwork], info, (ftnlen)1);
                        dgemm_((char *)"N", (char *)"N", m, n, n, &c_b79, &a[a_offset], lda, &work[iu], &ldwrku,
                               &c_b57, &u[u_offset], ldu, (ftnlen)1, (ftnlen)1);
                    } else {
                        itau = 1;
                        iwork = itau + *n;
                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, n, n, &u[u_offset], ldu, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        dlacpy_((char *)"U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (ftnlen)1);
                        if (*n > 1) {
                            i__2 = *n - 1;
                            i__3 = *n - 1;
                            dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &vt[vt_dim1 + 2], ldvt,
                                    (ftnlen)1);
                        }
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"Q", (char *)"R", (char *)"N", m, n, n, &vt[vt_offset], ldvt, &work[itauq],
                                &u[u_offset], ldu, &work[iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1,
                                (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &work[iwork],
                                &i__2, &ierr, (ftnlen)1);
                        iwork = ie + *n;
                        dbdsqr_((char *)"U", n, n, m, &c__0, &s[1], &work[ie], &vt[vt_offset], ldvt,
                                &u[u_offset], ldu, dum, &c__1, &work[iwork], info, (ftnlen)1);
                    }
                }
            } else if (wntua) {
                if (wntvn) {
                    i__2 = *n + *m, i__3 = *n << 2, i__2 = max(i__2, i__3);
                    if (*lwork >= *n * *n + max(i__2, bdspac)) {
                        ir = 1;
                        if (*lwork >= wrkbl + *lda * *n) {
                            ldwrkr = *lda;
                        } else {
                            ldwrkr = *n;
                        }
                        itau = ir + ldwrkr * *n;
                        iwork = itau + *n;
                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (ftnlen)1);
                        dlacpy_((char *)"U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr, (ftnlen)1);
                        i__2 = *n - 1;
                        i__3 = *n - 1;
                        dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 1], &ldwrkr,
                                (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        iwork = ie + *n;
                        dbdsqr_((char *)"U", n, &c__0, n, &c__0, &s[1], &work[ie], dum, &c__1, &work[ir],
                                &ldwrkr, dum, &c__1, &work[iwork], info, (ftnlen)1);
                        dgemm_((char *)"N", (char *)"N", m, n, n, &c_b79, &u[u_offset], ldu, &work[ir], &ldwrkr,
                               &c_b57, &a[a_offset], lda, (ftnlen)1, (ftnlen)1);
                        dlacpy_((char *)"F", m, n, &a[a_offset], lda, &u[u_offset], ldu, (ftnlen)1);
                    } else {
                        itau = 1;
                        iwork = itau + *n;
                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;
                        if (*n > 1) {
                            i__2 = *n - 1;
                            i__3 = *n - 1;
                            dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &a[a_dim1 + 2], lda,
                                    (ftnlen)1);
                        }
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"Q", (char *)"R", (char *)"N", m, n, n, &a[a_offset], lda, &work[itauq],
                                &u[u_offset], ldu, &work[iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1,
                                (ftnlen)1);
                        iwork = ie + *n;
                        dbdsqr_((char *)"U", n, &c__0, m, &c__0, &s[1], &work[ie], dum, &c__1, &u[u_offset],
                                ldu, dum, &c__1, &work[iwork], info, (ftnlen)1);
                    }
                } else if (wntvo) {
                    i__2 = *n + *m, i__3 = *n << 2, i__2 = max(i__2, i__3);
                    if (*lwork >= (*n << 1) * *n + max(i__2, bdspac)) {
                        iu = 1;
                        if (*lwork >= wrkbl + (*lda << 1) * *n) {
                            ldwrku = *lda;
                            ir = iu + ldwrku * *n;
                            ldwrkr = *lda;
                        } else if (*lwork >= wrkbl + (*lda + *n) * *n) {
                            ldwrku = *lda;
                            ir = iu + ldwrku * *n;
                            ldwrkr = *n;
                        } else {
                            ldwrku = *n;
                            ir = iu + ldwrku * *n;
                            ldwrkr = *n;
                        }
                        itau = ir + ldwrkr * *n;
                        iwork = itau + *n;
                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        dlacpy_((char *)"U", n, n, &a[a_offset], lda, &work[iu], &ldwrku, (ftnlen)1);
                        i__2 = *n - 1;
                        i__3 = *n - 1;
                        dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 1], &ldwrku,
                                (ftnlen)1);
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", n, n, &work[iu], &ldwrku, &work[ir], &ldwrkr, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", n, n, n, &work[iu], &ldwrku, &work[itauq], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", n, n, n, &work[ir], &ldwrkr, &work[itaup], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        iwork = ie + *n;
                        dbdsqr_((char *)"U", n, n, n, &c__0, &s[1], &work[ie], &work[ir], &ldwrkr,
                                &work[iu], &ldwrku, dum, &c__1, &work[iwork], info, (ftnlen)1);
                        dgemm_((char *)"N", (char *)"N", m, n, n, &c_b79, &u[u_offset], ldu, &work[iu], &ldwrku,
                               &c_b57, &a[a_offset], lda, (ftnlen)1, (ftnlen)1);
                        dlacpy_((char *)"F", m, n, &a[a_offset], lda, &u[u_offset], ldu, (ftnlen)1);
                        dlacpy_((char *)"F", n, n, &work[ir], &ldwrkr, &a[a_offset], lda, (ftnlen)1);
                    } else {
                        itau = 1;
                        iwork = itau + *n;
                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;
                        if (*n > 1) {
                            i__2 = *n - 1;
                            i__3 = *n - 1;
                            dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &a[a_dim1 + 2], lda,
                                    (ftnlen)1);
                        }
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"Q", (char *)"R", (char *)"N", m, n, n, &a[a_offset], lda, &work[itauq],
                                &u[u_offset], ldu, &work[iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1,
                                (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", n, n, n, &a[a_offset], lda, &work[itaup], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        iwork = ie + *n;
                        dbdsqr_((char *)"U", n, n, m, &c__0, &s[1], &work[ie], &a[a_offset], lda,
                                &u[u_offset], ldu, dum, &c__1, &work[iwork], info, (ftnlen)1);
                    }
                } else if (wntvas) {
                    i__2 = *n + *m, i__3 = *n << 2, i__2 = max(i__2, i__3);
                    if (*lwork >= *n * *n + max(i__2, bdspac)) {
                        iu = 1;
                        if (*lwork >= wrkbl + *lda * *n) {
                            ldwrku = *lda;
                        } else {
                            ldwrku = *n;
                        }
                        itau = iu + ldwrku * *n;
                        iwork = itau + *n;
                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        dlacpy_((char *)"U", n, n, &a[a_offset], lda, &work[iu], &ldwrku, (ftnlen)1);
                        i__2 = *n - 1;
                        i__3 = *n - 1;
                        dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &work[iu + 1], &ldwrku,
                                (ftnlen)1);
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", n, n, &work[iu], &ldwrku, &vt[vt_offset], ldvt, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", n, n, n, &work[iu], &ldwrku, &work[itauq], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &work[iwork],
                                &i__2, &ierr, (ftnlen)1);
                        iwork = ie + *n;
                        dbdsqr_((char *)"U", n, n, n, &c__0, &s[1], &work[ie], &vt[vt_offset], ldvt,
                                &work[iu], &ldwrku, dum, &c__1, &work[iwork], info, (ftnlen)1);
                        dgemm_((char *)"N", (char *)"N", m, n, n, &c_b79, &u[u_offset], ldu, &work[iu], &ldwrku,
                               &c_b57, &a[a_offset], lda, (ftnlen)1, (ftnlen)1);
                        dlacpy_((char *)"F", m, n, &a[a_offset], lda, &u[u_offset], ldu, (ftnlen)1);
                    } else {
                        itau = 1;
                        iwork = itau + *n;
                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        dlacpy_((char *)"U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (ftnlen)1);
                        if (*n > 1) {
                            i__2 = *n - 1;
                            i__3 = *n - 1;
                            dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &vt[vt_dim1 + 2], ldvt,
                                    (ftnlen)1);
                        }
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"Q", (char *)"R", (char *)"N", m, n, n, &vt[vt_offset], ldvt, &work[itauq],
                                &u[u_offset], ldu, &work[iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1,
                                (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &work[iwork],
                                &i__2, &ierr, (ftnlen)1);
                        iwork = ie + *n;
                        dbdsqr_((char *)"U", n, n, m, &c__0, &s[1], &work[ie], &vt[vt_offset], ldvt,
                                &u[u_offset], ldu, dum, &c__1, &work[iwork], info, (ftnlen)1);
                    }
                }
            }
        } else {
            ie = 1;
            itauq = ie + *n;
            itaup = itauq + *n;
            iwork = itaup + *n;
            i__2 = *lwork - iwork + 1;
            dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &work[itaup],
                    &work[iwork], &i__2, &ierr);
            if (wntuas) {
                dlacpy_((char *)"L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (ftnlen)1);
                if (wntus) {
                    ncu = *n;
                }
                if (wntua) {
                    ncu = *m;
                }
                i__2 = *lwork - iwork + 1;
                dorgbr_((char *)"Q", m, &ncu, n, &u[u_offset], ldu, &work[itauq], &work[iwork], &i__2,
                        &ierr, (ftnlen)1);
            }
            if (wntvas) {
                dlacpy_((char *)"U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (ftnlen)1);
                i__2 = *lwork - iwork + 1;
                dorgbr_((char *)"P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &work[iwork], &i__2,
                        &ierr, (ftnlen)1);
            }
            if (wntuo) {
                i__2 = *lwork - iwork + 1;
                dorgbr_((char *)"Q", m, n, n, &a[a_offset], lda, &work[itauq], &work[iwork], &i__2, &ierr,
                        (ftnlen)1);
            }
            if (wntvo) {
                i__2 = *lwork - iwork + 1;
                dorgbr_((char *)"P", n, n, n, &a[a_offset], lda, &work[itaup], &work[iwork], &i__2, &ierr,
                        (ftnlen)1);
            }
            iwork = ie + *n;
            if (wntuas || wntuo) {
                nru = *m;
            }
            if (wntun) {
                nru = 0;
            }
            if (wntvas || wntvo) {
                ncvt = *n;
            }
            if (wntvn) {
                ncvt = 0;
            }
            if (!wntuo && !wntvo) {
                dbdsqr_((char *)"U", n, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[vt_offset], ldvt,
                        &u[u_offset], ldu, dum, &c__1, &work[iwork], info, (ftnlen)1);
            } else if (!wntuo && wntvo) {
                dbdsqr_((char *)"U", n, &ncvt, &nru, &c__0, &s[1], &work[ie], &a[a_offset], lda,
                        &u[u_offset], ldu, dum, &c__1, &work[iwork], info, (ftnlen)1);
            } else {
                dbdsqr_((char *)"U", n, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[vt_offset], ldvt,
                        &a[a_offset], lda, dum, &c__1, &work[iwork], info, (ftnlen)1);
            }
        }
    } else {
        if (*n >= mnthr) {
            if (wntvn) {
                itau = 1;
                iwork = itau + *m;
                i__2 = *lwork - iwork + 1;
                dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                i__2 = *m - 1;
                i__3 = *m - 1;
                dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 << 1) + 1], lda, (ftnlen)1);
                ie = 1;
                itauq = ie + *m;
                itaup = itauq + *m;
                iwork = itaup + *m;
                i__2 = *lwork - iwork + 1;
                dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &work[itaup],
                        &work[iwork], &i__2, &ierr);
                if (wntuo || wntuas) {
                    i__2 = *lwork - iwork + 1;
                    dorgbr_((char *)"Q", m, m, m, &a[a_offset], lda, &work[itauq], &work[iwork], &i__2,
                            &ierr, (ftnlen)1);
                }
                iwork = ie + *m;
                nru = 0;
                if (wntuo || wntuas) {
                    nru = *m;
                }
                dbdsqr_((char *)"U", m, &c__0, &nru, &c__0, &s[1], &work[ie], dum, &c__1, &a[a_offset], lda,
                        dum, &c__1, &work[iwork], info, (ftnlen)1);
                if (wntuas) {
                    dlacpy_((char *)"F", m, m, &a[a_offset], lda, &u[u_offset], ldu, (ftnlen)1);
                }
            } else if (wntvo && wntun) {
                i__2 = *m << 2;
                if (*lwork >= *m * *m + max(i__2, bdspac)) {
                    ir = 1;
                    i__2 = wrkbl, i__3 = *lda * *n + *m;
                    if (*lwork >= max(i__2, i__3) + *lda * *m) {
                        ldwrku = *lda;
                        chunk = *n;
                        ldwrkr = *lda;
                    } else {
                        i__2 = wrkbl, i__3 = *lda * *n + *m;
                        if (*lwork >= max(i__2, i__3) + *m * *m) {
                            ldwrku = *lda;
                            chunk = *n;
                            ldwrkr = *m;
                        } else {
                            ldwrku = *m;
                            chunk = (*lwork - *m * *m - *m) / *m;
                            ldwrkr = *m;
                        }
                    }
                    itau = ir + ldwrkr * *m;
                    iwork = itau + *m;
                    i__2 = *lwork - iwork + 1;
                    dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                    dlacpy_((char *)"L", m, m, &a[a_offset], lda, &work[ir], &ldwrkr, (ftnlen)1);
                    i__2 = *m - 1;
                    i__3 = *m - 1;
                    dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &work[ir + ldwrkr], &ldwrkr,
                            (ftnlen)1);
                    i__2 = *lwork - iwork + 1;
                    dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                    ie = itau;
                    itauq = ie + *m;
                    itaup = itauq + *m;
                    iwork = itaup + *m;
                    i__2 = *lwork - iwork + 1;
                    dgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &work[ie], &work[itauq], &work[itaup],
                            &work[iwork], &i__2, &ierr);
                    i__2 = *lwork - iwork + 1;
                    dorgbr_((char *)"P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &work[iwork], &i__2,
                            &ierr, (ftnlen)1);
                    iwork = ie + *m;
                    dbdsqr_((char *)"U", m, m, &c__0, &c__0, &s[1], &work[ie], &work[ir], &ldwrkr, dum,
                            &c__1, dum, &c__1, &work[iwork], info, (ftnlen)1);
                    iu = ie + *m;
                    i__2 = *n;
                    i__3 = chunk;
                    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) {
                        i__4 = *n - i__ + 1;
                        blk = min(i__4, chunk);
                        dgemm_((char *)"N", (char *)"N", m, &blk, m, &c_b79, &work[ir], &ldwrkr,
                               &a[i__ * a_dim1 + 1], lda, &c_b57, &work[iu], &ldwrku, (ftnlen)1,
                               (ftnlen)1);
                        dlacpy_((char *)"F", m, &blk, &work[iu], &ldwrku, &a[i__ * a_dim1 + 1], lda,
                                (ftnlen)1);
                    }
                } else {
                    ie = 1;
                    itauq = ie + *m;
                    itaup = itauq + *m;
                    iwork = itaup + *m;
                    i__3 = *lwork - iwork + 1;
                    dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &work[itaup],
                            &work[iwork], &i__3, &ierr);
                    i__3 = *lwork - iwork + 1;
                    dorgbr_((char *)"P", m, n, m, &a[a_offset], lda, &work[itaup], &work[iwork], &i__3,
                            &ierr, (ftnlen)1);
                    iwork = ie + *m;
                    dbdsqr_((char *)"L", m, n, &c__0, &c__0, &s[1], &work[ie], &a[a_offset], lda, dum,
                            &c__1, dum, &c__1, &work[iwork], info, (ftnlen)1);
                }
            } else if (wntvo && wntuas) {
                i__3 = *m << 2;
                if (*lwork >= *m * *m + max(i__3, bdspac)) {
                    ir = 1;
                    i__3 = wrkbl, i__2 = *lda * *n + *m;
                    if (*lwork >= max(i__3, i__2) + *lda * *m) {
                        ldwrku = *lda;
                        chunk = *n;
                        ldwrkr = *lda;
                    } else {
                        i__3 = wrkbl, i__2 = *lda * *n + *m;
                        if (*lwork >= max(i__3, i__2) + *m * *m) {
                            ldwrku = *lda;
                            chunk = *n;
                            ldwrkr = *m;
                        } else {
                            ldwrku = *m;
                            chunk = (*lwork - *m * *m - *m) / *m;
                            ldwrkr = *m;
                        }
                    }
                    itau = ir + ldwrkr * *m;
                    iwork = itau + *m;
                    i__3 = *lwork - iwork + 1;
                    dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__3, &ierr);
                    dlacpy_((char *)"L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (ftnlen)1);
                    i__3 = *m - 1;
                    i__2 = *m - 1;
                    dlaset_((char *)"U", &i__3, &i__2, &c_b57, &c_b57, &u[(u_dim1 << 1) + 1], ldu,
                            (ftnlen)1);
                    i__3 = *lwork - iwork + 1;
                    dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[iwork], &i__3, &ierr);
                    ie = itau;
                    itauq = ie + *m;
                    itaup = itauq + *m;
                    iwork = itaup + *m;
                    i__3 = *lwork - iwork + 1;
                    dgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &work[itauq], &work[itaup],
                            &work[iwork], &i__3, &ierr);
                    dlacpy_((char *)"U", m, m, &u[u_offset], ldu, &work[ir], &ldwrkr, (ftnlen)1);
                    i__3 = *lwork - iwork + 1;
                    dorgbr_((char *)"P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &work[iwork], &i__3,
                            &ierr, (ftnlen)1);
                    i__3 = *lwork - iwork + 1;
                    dorgbr_((char *)"Q", m, m, m, &u[u_offset], ldu, &work[itauq], &work[iwork], &i__3,
                            &ierr, (ftnlen)1);
                    iwork = ie + *m;
                    dbdsqr_((char *)"U", m, m, m, &c__0, &s[1], &work[ie], &work[ir], &ldwrkr, &u[u_offset],
                            ldu, dum, &c__1, &work[iwork], info, (ftnlen)1);
                    iu = ie + *m;
                    i__3 = *n;
                    i__2 = chunk;
                    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ += i__2) {
                        i__4 = *n - i__ + 1;
                        blk = min(i__4, chunk);
                        dgemm_((char *)"N", (char *)"N", m, &blk, m, &c_b79, &work[ir], &ldwrkr,
                               &a[i__ * a_dim1 + 1], lda, &c_b57, &work[iu], &ldwrku, (ftnlen)1,
                               (ftnlen)1);
                        dlacpy_((char *)"F", m, &blk, &work[iu], &ldwrku, &a[i__ * a_dim1 + 1], lda,
                                (ftnlen)1);
                    }
                } else {
                    itau = 1;
                    iwork = itau + *m;
                    i__2 = *lwork - iwork + 1;
                    dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                    dlacpy_((char *)"L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (ftnlen)1);
                    i__2 = *m - 1;
                    i__3 = *m - 1;
                    dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &u[(u_dim1 << 1) + 1], ldu,
                            (ftnlen)1);
                    i__2 = *lwork - iwork + 1;
                    dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                    ie = itau;
                    itauq = ie + *m;
                    itaup = itauq + *m;
                    iwork = itaup + *m;
                    i__2 = *lwork - iwork + 1;
                    dgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &work[itauq], &work[itaup],
                            &work[iwork], &i__2, &ierr);
                    i__2 = *lwork - iwork + 1;
                    dormbr_((char *)"P", (char *)"L", (char *)"T", m, n, m, &u[u_offset], ldu, &work[itaup], &a[a_offset],
                            lda, &work[iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                    i__2 = *lwork - iwork + 1;
                    dorgbr_((char *)"Q", m, m, m, &u[u_offset], ldu, &work[itauq], &work[iwork], &i__2,
                            &ierr, (ftnlen)1);
                    iwork = ie + *m;
                    dbdsqr_((char *)"U", m, n, m, &c__0, &s[1], &work[ie], &a[a_offset], lda, &u[u_offset],
                            ldu, dum, &c__1, &work[iwork], info, (ftnlen)1);
                }
            } else if (wntvs) {
                if (wntun) {
                    i__2 = *m << 2;
                    if (*lwork >= *m * *m + max(i__2, bdspac)) {
                        ir = 1;
                        if (*lwork >= wrkbl + *lda * *m) {
                            ldwrkr = *lda;
                        } else {
                            ldwrkr = *m;
                        }
                        itau = ir + ldwrkr * *m;
                        iwork = itau + *m;
                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, m, &a[a_offset], lda, &work[ir], &ldwrkr, (ftnlen)1);
                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &work[ir + ldwrkr], &ldwrkr,
                                (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        iwork = ie + *m;
                        dbdsqr_((char *)"U", m, m, &c__0, &c__0, &s[1], &work[ie], &work[ir], &ldwrkr, dum,
                                &c__1, dum, &c__1, &work[iwork], info, (ftnlen)1);
                        dgemm_((char *)"N", (char *)"N", m, n, m, &c_b79, &work[ir], &ldwrkr, &a[a_offset], lda,
                               &c_b57, &vt[vt_offset], ldvt, (ftnlen)1, (ftnlen)1);
                    } else {
                        itau = 1;
                        iwork = itau + *m;
                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;
                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 << 1) + 1], lda,
                                (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"P", (char *)"L", (char *)"T", m, n, m, &a[a_offset], lda, &work[itaup],
                                &vt[vt_offset], ldvt, &work[iwork], &i__2, &ierr, (ftnlen)1,
                                (ftnlen)1, (ftnlen)1);
                        iwork = ie + *m;
                        dbdsqr_((char *)"U", m, n, &c__0, &c__0, &s[1], &work[ie], &vt[vt_offset], ldvt,
                                dum, &c__1, dum, &c__1, &work[iwork], info, (ftnlen)1);
                    }
                } else if (wntuo) {
                    i__2 = *m << 2;
                    if (*lwork >= (*m << 1) * *m + max(i__2, bdspac)) {
                        iu = 1;
                        if (*lwork >= wrkbl + (*lda << 1) * *m) {
                            ldwrku = *lda;
                            ir = iu + ldwrku * *m;
                            ldwrkr = *lda;
                        } else if (*lwork >= wrkbl + (*lda + *m) * *m) {
                            ldwrku = *lda;
                            ir = iu + ldwrku * *m;
                            ldwrkr = *m;
                        } else {
                            ldwrku = *m;
                            ir = iu + ldwrku * *m;
                            ldwrkr = *m;
                        }
                        itau = ir + ldwrkr * *m;
                        iwork = itau + *m;
                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, m, &a[a_offset], lda, &work[iu], &ldwrku, (ftnlen)1);
                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &work[iu + ldwrku], &ldwrku,
                                (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, m, &work[iu], &ldwrku, &work[ir], &ldwrkr, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", m, m, m, &work[iu], &ldwrku, &work[itaup], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", m, m, m, &work[ir], &ldwrkr, &work[itauq], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        iwork = ie + *m;
                        dbdsqr_((char *)"U", m, m, m, &c__0, &s[1], &work[ie], &work[iu], &ldwrku,
                                &work[ir], &ldwrkr, dum, &c__1, &work[iwork], info, (ftnlen)1);
                        dgemm_((char *)"N", (char *)"N", m, n, m, &c_b79, &work[iu], &ldwrku, &a[a_offset], lda,
                               &c_b57, &vt[vt_offset], ldvt, (ftnlen)1, (ftnlen)1);
                        dlacpy_((char *)"F", m, m, &work[ir], &ldwrkr, &a[a_offset], lda, (ftnlen)1);
                    } else {
                        itau = 1;
                        iwork = itau + *m;
                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;
                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 << 1) + 1], lda,
                                (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"P", (char *)"L", (char *)"T", m, n, m, &a[a_offset], lda, &work[itaup],
                                &vt[vt_offset], ldvt, &work[iwork], &i__2, &ierr, (ftnlen)1,
                                (ftnlen)1, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", m, m, m, &a[a_offset], lda, &work[itauq], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        iwork = ie + *m;
                        dbdsqr_((char *)"U", m, n, m, &c__0, &s[1], &work[ie], &vt[vt_offset], ldvt,
                                &a[a_offset], lda, dum, &c__1, &work[iwork], info, (ftnlen)1);
                    }
                } else if (wntuas) {
                    i__2 = *m << 2;
                    if (*lwork >= *m * *m + max(i__2, bdspac)) {
                        iu = 1;
                        if (*lwork >= wrkbl + *lda * *m) {
                            ldwrku = *lda;
                        } else {
                            ldwrku = *m;
                        }
                        itau = iu + ldwrku * *m;
                        iwork = itau + *m;
                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, m, &a[a_offset], lda, &work[iu], &ldwrku, (ftnlen)1);
                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &work[iu + ldwrku], &ldwrku,
                                (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, m, &work[iu], &ldwrku, &u[u_offset], ldu, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", m, m, m, &work[iu], &ldwrku, &work[itaup], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", m, m, m, &u[u_offset], ldu, &work[itauq], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        iwork = ie + *m;
                        dbdsqr_((char *)"U", m, m, m, &c__0, &s[1], &work[ie], &work[iu], &ldwrku,
                                &u[u_offset], ldu, dum, &c__1, &work[iwork], info, (ftnlen)1);
                        dgemm_((char *)"N", (char *)"N", m, n, m, &c_b79, &work[iu], &ldwrku, &a[a_offset], lda,
                               &c_b57, &vt[vt_offset], ldvt, (ftnlen)1, (ftnlen)1);
                    } else {
                        itau = 1;
                        iwork = itau + *m;
                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        dlacpy_((char *)"L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (ftnlen)1);
                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &u[(u_dim1 << 1) + 1], ldu,
                                (ftnlen)1);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"P", (char *)"L", (char *)"T", m, n, m, &u[u_offset], ldu, &work[itaup],
                                &vt[vt_offset], ldvt, &work[iwork], &i__2, &ierr, (ftnlen)1,
                                (ftnlen)1, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", m, m, m, &u[u_offset], ldu, &work[itauq], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        iwork = ie + *m;
                        dbdsqr_((char *)"U", m, n, m, &c__0, &s[1], &work[ie], &vt[vt_offset], ldvt,
                                &u[u_offset], ldu, dum, &c__1, &work[iwork], info, (ftnlen)1);
                    }
                }
            } else if (wntva) {
                if (wntun) {
                    i__2 = *n + *m, i__3 = *m << 2, i__2 = max(i__2, i__3);
                    if (*lwork >= *m * *m + max(i__2, bdspac)) {
                        ir = 1;
                        if (*lwork >= wrkbl + *lda * *m) {
                            ldwrkr = *lda;
                        } else {
                            ldwrkr = *m;
                        }
                        itau = ir + ldwrkr * *m;
                        iwork = itau + *m;
                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (ftnlen)1);
                        dlacpy_((char *)"L", m, m, &a[a_offset], lda, &work[ir], &ldwrkr, (ftnlen)1);
                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &work[ir + ldwrkr], &ldwrkr,
                                (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        iwork = ie + *m;
                        dbdsqr_((char *)"U", m, m, &c__0, &c__0, &s[1], &work[ie], &work[ir], &ldwrkr, dum,
                                &c__1, dum, &c__1, &work[iwork], info, (ftnlen)1);
                        dgemm_((char *)"N", (char *)"N", m, n, m, &c_b79, &work[ir], &ldwrkr, &vt[vt_offset], ldvt,
                               &c_b57, &a[a_offset], lda, (ftnlen)1, (ftnlen)1);
                        dlacpy_((char *)"F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (ftnlen)1);
                    } else {
                        itau = 1;
                        iwork = itau + *m;
                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;
                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 << 1) + 1], lda,
                                (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"P", (char *)"L", (char *)"T", m, n, m, &a[a_offset], lda, &work[itaup],
                                &vt[vt_offset], ldvt, &work[iwork], &i__2, &ierr, (ftnlen)1,
                                (ftnlen)1, (ftnlen)1);
                        iwork = ie + *m;
                        dbdsqr_((char *)"U", m, n, &c__0, &c__0, &s[1], &work[ie], &vt[vt_offset], ldvt,
                                dum, &c__1, dum, &c__1, &work[iwork], info, (ftnlen)1);
                    }
                } else if (wntuo) {
                    i__2 = *n + *m, i__3 = *m << 2, i__2 = max(i__2, i__3);
                    if (*lwork >= (*m << 1) * *m + max(i__2, bdspac)) {
                        iu = 1;
                        if (*lwork >= wrkbl + (*lda << 1) * *m) {
                            ldwrku = *lda;
                            ir = iu + ldwrku * *m;
                            ldwrkr = *lda;
                        } else if (*lwork >= wrkbl + (*lda + *m) * *m) {
                            ldwrku = *lda;
                            ir = iu + ldwrku * *m;
                            ldwrkr = *m;
                        } else {
                            ldwrku = *m;
                            ir = iu + ldwrku * *m;
                            ldwrkr = *m;
                        }
                        itau = ir + ldwrkr * *m;
                        iwork = itau + *m;
                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        dlacpy_((char *)"L", m, m, &a[a_offset], lda, &work[iu], &ldwrku, (ftnlen)1);
                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &work[iu + ldwrku], &ldwrku,
                                (ftnlen)1);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, m, &work[iu], &ldwrku, &work[ir], &ldwrkr, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", m, m, m, &work[iu], &ldwrku, &work[itaup], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", m, m, m, &work[ir], &ldwrkr, &work[itauq], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        iwork = ie + *m;
                        dbdsqr_((char *)"U", m, m, m, &c__0, &s[1], &work[ie], &work[iu], &ldwrku,
                                &work[ir], &ldwrkr, dum, &c__1, &work[iwork], info, (ftnlen)1);
                        dgemm_((char *)"N", (char *)"N", m, n, m, &c_b79, &work[iu], &ldwrku, &vt[vt_offset], ldvt,
                               &c_b57, &a[a_offset], lda, (ftnlen)1, (ftnlen)1);
                        dlacpy_((char *)"F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (ftnlen)1);
                        dlacpy_((char *)"F", m, m, &work[ir], &ldwrkr, &a[a_offset], lda, (ftnlen)1);
                    } else {
                        itau = 1;
                        iwork = itau + *m;
                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;
                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 << 1) + 1], lda,
                                (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"P", (char *)"L", (char *)"T", m, n, m, &a[a_offset], lda, &work[itaup],
                                &vt[vt_offset], ldvt, &work[iwork], &i__2, &ierr, (ftnlen)1,
                                (ftnlen)1, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", m, m, m, &a[a_offset], lda, &work[itauq], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        iwork = ie + *m;
                        dbdsqr_((char *)"U", m, n, m, &c__0, &s[1], &work[ie], &vt[vt_offset], ldvt,
                                &a[a_offset], lda, dum, &c__1, &work[iwork], info, (ftnlen)1);
                    }
                } else if (wntuas) {
                    i__2 = *n + *m, i__3 = *m << 2, i__2 = max(i__2, i__3);
                    if (*lwork >= *m * *m + max(i__2, bdspac)) {
                        iu = 1;
                        if (*lwork >= wrkbl + *lda * *m) {
                            ldwrku = *lda;
                        } else {
                            ldwrku = *m;
                        }
                        itau = iu + ldwrku * *m;
                        iwork = itau + *m;
                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        dlacpy_((char *)"L", m, m, &a[a_offset], lda, &work[iu], &ldwrku, (ftnlen)1);
                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &work[iu + ldwrku], &ldwrku,
                                (ftnlen)1);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, m, &work[iu], &ldwrku, &u[u_offset], ldu, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", m, m, m, &work[iu], &ldwrku, &work[itaup], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", m, m, m, &u[u_offset], ldu, &work[itauq], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        iwork = ie + *m;
                        dbdsqr_((char *)"U", m, m, m, &c__0, &s[1], &work[ie], &work[iu], &ldwrku,
                                &u[u_offset], ldu, dum, &c__1, &work[iwork], info, (ftnlen)1);
                        dgemm_((char *)"N", (char *)"N", m, n, m, &c_b79, &work[iu], &ldwrku, &vt[vt_offset], ldvt,
                               &c_b57, &a[a_offset], lda, (ftnlen)1, (ftnlen)1);
                        dlacpy_((char *)"F", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (ftnlen)1);
                    } else {
                        itau = 1;
                        iwork = itau + *m;
                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &work[iwork], &i__2,
                                &ierr);
                        dlacpy_((char *)"L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (ftnlen)1);
                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &u[(u_dim1 << 1) + 1], ldu,
                                (ftnlen)1);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;
                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &work[itauq],
                                &work[itaup], &work[iwork], &i__2, &ierr);
                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"P", (char *)"L", (char *)"T", m, n, m, &u[u_offset], ldu, &work[itaup],
                                &vt[vt_offset], ldvt, &work[iwork], &i__2, &ierr, (ftnlen)1,
                                (ftnlen)1, (ftnlen)1);
                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", m, m, m, &u[u_offset], ldu, &work[itauq], &work[iwork], &i__2,
                                &ierr, (ftnlen)1);
                        iwork = ie + *m;
                        dbdsqr_((char *)"U", m, n, m, &c__0, &s[1], &work[ie], &vt[vt_offset], ldvt,
                                &u[u_offset], ldu, dum, &c__1, &work[iwork], info, (ftnlen)1);
                    }
                }
            }
        } else {
            ie = 1;
            itauq = ie + *m;
            itaup = itauq + *m;
            iwork = itaup + *m;
            i__2 = *lwork - iwork + 1;
            dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &work[itaup],
                    &work[iwork], &i__2, &ierr);
            if (wntuas) {
                dlacpy_((char *)"L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (ftnlen)1);
                i__2 = *lwork - iwork + 1;
                dorgbr_((char *)"Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[iwork], &i__2, &ierr,
                        (ftnlen)1);
            }
            if (wntvas) {
                dlacpy_((char *)"U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (ftnlen)1);
                if (wntva) {
                    nrvt = *n;
                }
                if (wntvs) {
                    nrvt = *m;
                }
                i__2 = *lwork - iwork + 1;
                dorgbr_((char *)"P", &nrvt, n, m, &vt[vt_offset], ldvt, &work[itaup], &work[iwork], &i__2,
                        &ierr, (ftnlen)1);
            }
            if (wntuo) {
                i__2 = *lwork - iwork + 1;
                dorgbr_((char *)"Q", m, m, n, &a[a_offset], lda, &work[itauq], &work[iwork], &i__2, &ierr,
                        (ftnlen)1);
            }
            if (wntvo) {
                i__2 = *lwork - iwork + 1;
                dorgbr_((char *)"P", m, n, m, &a[a_offset], lda, &work[itaup], &work[iwork], &i__2, &ierr,
                        (ftnlen)1);
            }
            iwork = ie + *m;
            if (wntuas || wntuo) {
                nru = *m;
            }
            if (wntun) {
                nru = 0;
            }
            if (wntvas || wntvo) {
                ncvt = *n;
            }
            if (wntvn) {
                ncvt = 0;
            }
            if (!wntuo && !wntvo) {
                dbdsqr_((char *)"L", m, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[vt_offset], ldvt,
                        &u[u_offset], ldu, dum, &c__1, &work[iwork], info, (ftnlen)1);
            } else if (!wntuo && wntvo) {
                dbdsqr_((char *)"L", m, &ncvt, &nru, &c__0, &s[1], &work[ie], &a[a_offset], lda,
                        &u[u_offset], ldu, dum, &c__1, &work[iwork], info, (ftnlen)1);
            } else {
                dbdsqr_((char *)"L", m, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[vt_offset], ldvt,
                        &a[a_offset], lda, dum, &c__1, &work[iwork], info, (ftnlen)1);
            }
        }
    }
    if (*info != 0) {
        if (ie > 2) {
            i__2 = minmn - 1;
            for (i__ = 1; i__ <= i__2; ++i__) {
                work[i__ + 1] = work[i__ + ie - 1];
            }
        }
        if (ie < 2) {
            for (i__ = minmn - 1; i__ >= 1; --i__) {
                work[i__ + 1] = work[i__ + ie - 1];
            }
        }
    }
    if (iscl == 1) {
        if (anrm > bignum) {
            dlascl_((char *)"G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &minmn, &ierr,
                    (ftnlen)1);
        }
        if (*info != 0 && anrm > bignum) {
            i__2 = minmn - 1;
            dlascl_((char *)"G", &c__0, &c__0, &bignum, &anrm, &i__2, &c__1, &work[2], &minmn, &ierr,
                    (ftnlen)1);
        }
        if (anrm < smlnum) {
            dlascl_((char *)"G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &minmn, &ierr,
                    (ftnlen)1);
        }
        if (*info != 0 && anrm < smlnum) {
            i__2 = minmn - 1;
            dlascl_((char *)"G", &c__0, &c__0, &smlnum, &anrm, &i__2, &c__1, &work[2], &minmn, &ierr,
                    (ftnlen)1);
        }
    }
    work[1] = (doublereal)maxwrk;
    return 0;
}
#ifdef __cplusplus
}
#endif
