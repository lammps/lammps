#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__6 = 6;
static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b46 = 0.;
static integer c__1 = 1;
static doublereal c_b79 = 1.;
int dgelss_(integer *m, integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *b,
            integer *ldb, doublereal *s, doublereal *rcond, integer *rank, doublereal *work,
            integer *lwork, integer *info)
{
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;
    integer i__, bl, ie, il, mm;
    doublereal dum[1], eps, thr, anrm, bnrm;
    integer itau, lwork_dgebrd__, lwork_dgelqf__, lwork_dgeqrf__, lwork_dorgbr__, lwork_dormbr__,
        lwork_dormlq__, lwork_dormqr__;
    extern int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
                      integer *, doublereal *, integer *, doublereal *, doublereal *, integer *,
                      ftnlen, ftnlen);
    integer iascl, ibscl;
    extern int dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *,
                      doublereal *, integer *, doublereal *, doublereal *, integer *, ftnlen),
        drscl_(integer *, doublereal *, doublereal *, integer *);
    integer chunk;
    doublereal sfmin;
    integer minmn;
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer maxmn, itaup, itauq, mnthr, iwork;
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
        xerbla_(char *, integer *, ftnlen),
        dbdsqr_(char *, integer *, integer *, integer *, integer *, doublereal *, doublereal *,
                doublereal *, integer *, doublereal *, integer *, doublereal *, integer *,
                doublereal *, integer *, ftnlen),
        dorgbr_(char *, integer *, integer *, integer *, doublereal *, integer *, doublereal *,
                doublereal *, integer *, integer *, ftnlen);
    doublereal bignum;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    extern int dormbr_(char *, char *, char *, integer *, integer *, integer *, doublereal *,
                       integer *, doublereal *, doublereal *, integer *, doublereal *, integer *,
                       integer *, ftnlen, ftnlen, ftnlen),
        dormlq_(char *, char *, integer *, integer *, integer *, doublereal *, integer *,
                doublereal *, doublereal *, integer *, doublereal *, integer *, integer *, ftnlen,
                ftnlen);
    integer ldwork;
    extern int dormqr_(char *, char *, integer *, integer *, integer *, doublereal *, integer *,
                       doublereal *, doublereal *, integer *, doublereal *, integer *, integer *,
                       ftnlen, ftnlen);
    integer minwrk, maxwrk;
    doublereal smlnum;
    logical lquery;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --s;
    --work;
    *info = 0;
    minmn = min(*m, *n);
    maxmn = max(*m, *n);
    lquery = *lwork == -1;
    if (*m < 0) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*nrhs < 0) {
        *info = -3;
    } else if (*lda < max(1, *m)) {
        *info = -5;
    } else if (*ldb < max(1, maxmn)) {
        *info = -7;
    }
    if (*info == 0) {
        minwrk = 1;
        maxwrk = 1;
        if (minmn > 0) {
            mm = *m;
            mnthr = ilaenv_(&c__6, (char *)"DGELSS", (char *)" ", m, n, nrhs, &c_n1, (ftnlen)6, (ftnlen)1);
            if (*m >= *n && *m >= mnthr) {
                dgeqrf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, info);
                lwork_dgeqrf__ = (integer)dum[0];
                dormqr_((char *)"L", (char *)"T", m, nrhs, n, &a[a_offset], lda, dum, &b[b_offset], ldb, dum, &c_n1,
                        info, (ftnlen)1, (ftnlen)1);
                lwork_dormqr__ = (integer)dum[0];
                mm = *n;
                i__1 = maxwrk, i__2 = *n + lwork_dgeqrf__;
                maxwrk = max(i__1, i__2);
                i__1 = maxwrk, i__2 = *n + lwork_dormqr__;
                maxwrk = max(i__1, i__2);
            }
            if (*m >= *n) {
                i__1 = 1, i__2 = *n * 5;
                bdspac = max(i__1, i__2);
                dgebrd_(&mm, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &c_n1, info);
                lwork_dgebrd__ = (integer)dum[0];
                dormbr_((char *)"Q", (char *)"L", (char *)"T", &mm, nrhs, n, &a[a_offset], lda, dum, &b[b_offset], ldb, dum,
                        &c_n1, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                lwork_dormbr__ = (integer)dum[0];
                dorgbr_((char *)"P", n, n, n, &a[a_offset], lda, dum, dum, &c_n1, info, (ftnlen)1);
                lwork_dorgbr__ = (integer)dum[0];
                i__1 = maxwrk, i__2 = *n * 3 + lwork_dgebrd__;
                maxwrk = max(i__1, i__2);
                i__1 = maxwrk, i__2 = *n * 3 + lwork_dormbr__;
                maxwrk = max(i__1, i__2);
                i__1 = maxwrk, i__2 = *n * 3 + lwork_dorgbr__;
                maxwrk = max(i__1, i__2);
                maxwrk = max(maxwrk, bdspac);
                i__1 = maxwrk, i__2 = *n * *nrhs;
                maxwrk = max(i__1, i__2);
                i__1 = *n * 3 + mm, i__2 = *n * 3 + *nrhs, i__1 = max(i__1, i__2);
                minwrk = max(i__1, bdspac);
                maxwrk = max(minwrk, maxwrk);
            }
            if (*n > *m) {
                i__1 = 1, i__2 = *m * 5;
                bdspac = max(i__1, i__2);
                i__1 = *m * 3 + *nrhs, i__2 = *m * 3 + *n, i__1 = max(i__1, i__2);
                minwrk = max(i__1, bdspac);
                if (*n >= mnthr) {
                    dgelqf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, info);
                    lwork_dgelqf__ = (integer)dum[0];
                    dgebrd_(m, m, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &c_n1, info);
                    lwork_dgebrd__ = (integer)dum[0];
                    dormbr_((char *)"Q", (char *)"L", (char *)"T", m, nrhs, n, &a[a_offset], lda, dum, &b[b_offset], ldb,
                            dum, &c_n1, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                    lwork_dormbr__ = (integer)dum[0];
                    dorgbr_((char *)"P", m, m, m, &a[a_offset], lda, dum, dum, &c_n1, info, (ftnlen)1);
                    lwork_dorgbr__ = (integer)dum[0];
                    dormlq_((char *)"L", (char *)"T", n, nrhs, m, &a[a_offset], lda, dum, &b[b_offset], ldb, dum,
                            &c_n1, info, (ftnlen)1, (ftnlen)1);
                    lwork_dormlq__ = (integer)dum[0];
                    maxwrk = *m + lwork_dgelqf__;
                    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + lwork_dgebrd__;
                    maxwrk = max(i__1, i__2);
                    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + lwork_dormbr__;
                    maxwrk = max(i__1, i__2);
                    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + lwork_dorgbr__;
                    maxwrk = max(i__1, i__2);
                    i__1 = maxwrk, i__2 = *m * *m + *m + bdspac;
                    maxwrk = max(i__1, i__2);
                    if (*nrhs > 1) {
                        i__1 = maxwrk, i__2 = *m * *m + *m + *m * *nrhs;
                        maxwrk = max(i__1, i__2);
                    } else {
                        i__1 = maxwrk, i__2 = *m * *m + (*m << 1);
                        maxwrk = max(i__1, i__2);
                    }
                    i__1 = maxwrk, i__2 = *m + lwork_dormlq__;
                    maxwrk = max(i__1, i__2);
                } else {
                    dgebrd_(m, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &c_n1, info);
                    lwork_dgebrd__ = (integer)dum[0];
                    dormbr_((char *)"Q", (char *)"L", (char *)"T", m, nrhs, m, &a[a_offset], lda, dum, &b[b_offset], ldb,
                            dum, &c_n1, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                    lwork_dormbr__ = (integer)dum[0];
                    dorgbr_((char *)"P", m, n, m, &a[a_offset], lda, dum, dum, &c_n1, info, (ftnlen)1);
                    lwork_dorgbr__ = (integer)dum[0];
                    maxwrk = *m * 3 + lwork_dgebrd__;
                    i__1 = maxwrk, i__2 = *m * 3 + lwork_dormbr__;
                    maxwrk = max(i__1, i__2);
                    i__1 = maxwrk, i__2 = *m * 3 + lwork_dorgbr__;
                    maxwrk = max(i__1, i__2);
                    maxwrk = max(maxwrk, bdspac);
                    i__1 = maxwrk, i__2 = *n * *nrhs;
                    maxwrk = max(i__1, i__2);
                }
            }
            maxwrk = max(minwrk, maxwrk);
        }
        work[1] = (doublereal)maxwrk;
        if (*lwork < minwrk && !lquery) {
            *info = -12;
        }
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DGELSS", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    if (*m == 0 || *n == 0) {
        *rank = 0;
        return 0;
    }
    eps = dlamch_((char *)"P", (ftnlen)1);
    sfmin = dlamch_((char *)"S", (ftnlen)1);
    smlnum = sfmin / eps;
    bignum = 1. / smlnum;
    anrm = dlange_((char *)"M", m, n, &a[a_offset], lda, &work[1], (ftnlen)1);
    iascl = 0;
    if (anrm > 0. && anrm < smlnum) {
        dlascl_((char *)"G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, info, (ftnlen)1);
        iascl = 1;
    } else if (anrm > bignum) {
        dlascl_((char *)"G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, info, (ftnlen)1);
        iascl = 2;
    } else if (anrm == 0.) {
        i__1 = max(*m, *n);
        dlaset_((char *)"F", &i__1, nrhs, &c_b46, &c_b46, &b[b_offset], ldb, (ftnlen)1);
        dlaset_((char *)"F", &minmn, &c__1, &c_b46, &c_b46, &s[1], &minmn, (ftnlen)1);
        *rank = 0;
        goto L70;
    }
    bnrm = dlange_((char *)"M", m, nrhs, &b[b_offset], ldb, &work[1], (ftnlen)1);
    ibscl = 0;
    if (bnrm > 0. && bnrm < smlnum) {
        dlascl_((char *)"G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &b[b_offset], ldb, info, (ftnlen)1);
        ibscl = 1;
    } else if (bnrm > bignum) {
        dlascl_((char *)"G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &b[b_offset], ldb, info, (ftnlen)1);
        ibscl = 2;
    }
    if (*m >= *n) {
        mm = *m;
        if (*m >= mnthr) {
            mm = *n;
            itau = 1;
            iwork = itau + *n;
            i__1 = *lwork - iwork + 1;
            dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__1, info);
            i__1 = *lwork - iwork + 1;
            dormqr_((char *)"L", (char *)"T", m, nrhs, n, &a[a_offset], lda, &work[itau], &b[b_offset], ldb,
                    &work[iwork], &i__1, info, (ftnlen)1, (ftnlen)1);
            if (*n > 1) {
                i__1 = *n - 1;
                i__2 = *n - 1;
                dlaset_((char *)"L", &i__1, &i__2, &c_b46, &c_b46, &a[a_dim1 + 2], lda, (ftnlen)1);
            }
        }
        ie = 1;
        itauq = ie + *n;
        itaup = itauq + *n;
        iwork = itaup + *n;
        i__1 = *lwork - iwork + 1;
        dgebrd_(&mm, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &work[itaup],
                &work[iwork], &i__1, info);
        i__1 = *lwork - iwork + 1;
        dormbr_((char *)"Q", (char *)"L", (char *)"T", &mm, nrhs, n, &a[a_offset], lda, &work[itauq], &b[b_offset], ldb,
                &work[iwork], &i__1, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
        i__1 = *lwork - iwork + 1;
        dorgbr_((char *)"P", n, n, n, &a[a_offset], lda, &work[itaup], &work[iwork], &i__1, info,
                (ftnlen)1);
        iwork = ie + *n;
        dbdsqr_((char *)"U", n, n, &c__0, nrhs, &s[1], &work[ie], &a[a_offset], lda, dum, &c__1,
                &b[b_offset], ldb, &work[iwork], info, (ftnlen)1);
        if (*info != 0) {
            goto L70;
        }
        d__1 = *rcond * s[1];
        thr = max(d__1, sfmin);
        if (*rcond < 0.) {
            d__1 = eps * s[1];
            thr = max(d__1, sfmin);
        }
        *rank = 0;
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            if (s[i__] > thr) {
                drscl_(nrhs, &s[i__], &b[i__ + b_dim1], ldb);
                ++(*rank);
            } else {
                dlaset_((char *)"F", &c__1, nrhs, &c_b46, &c_b46, &b[i__ + b_dim1], ldb, (ftnlen)1);
            }
        }
        if (*lwork >= *ldb * *nrhs && *nrhs > 1) {
            dgemm_((char *)"T", (char *)"N", n, nrhs, n, &c_b79, &a[a_offset], lda, &b[b_offset], ldb, &c_b46,
                   &work[1], ldb, (ftnlen)1, (ftnlen)1);
            dlacpy_((char *)"G", n, nrhs, &work[1], ldb, &b[b_offset], ldb, (ftnlen)1);
        } else if (*nrhs > 1) {
            chunk = *lwork / *n;
            i__1 = *nrhs;
            i__2 = chunk;
            for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
                i__3 = *nrhs - i__ + 1;
                bl = min(i__3, chunk);
                dgemm_((char *)"T", (char *)"N", n, &bl, n, &c_b79, &a[a_offset], lda, &b[i__ * b_dim1 + 1], ldb,
                       &c_b46, &work[1], n, (ftnlen)1, (ftnlen)1);
                dlacpy_((char *)"G", n, &bl, &work[1], n, &b[i__ * b_dim1 + 1], ldb, (ftnlen)1);
            }
        } else if (*nrhs == 1) {
            dgemv_((char *)"T", n, n, &c_b79, &a[a_offset], lda, &b[b_offset], &c__1, &c_b46, &work[1],
                   &c__1, (ftnlen)1);
            dcopy_(n, &work[1], &c__1, &b[b_offset], &c__1);
        }
    } else {
        i__2 = *m, i__1 = (*m << 1) - 4, i__2 = max(i__2, i__1), i__2 = max(i__2, *nrhs),
        i__1 = *n - *m * 3;
        if (*n >= mnthr && *lwork >= (*m << 2) + *m * *m + max(i__2, i__1)) {
            ldwork = *m;
            i__3 = *m, i__4 = (*m << 1) - 4, i__3 = max(i__3, i__4), i__3 = max(i__3, *nrhs),
            i__4 = *n - *m * 3;
            i__2 = (*m << 2) + *m * *lda + max(i__3, i__4), i__1 = *m * *lda + *m + *m * *nrhs;
            if (*lwork >= max(i__2, i__1)) {
                ldwork = *lda;
            }
            itau = 1;
            iwork = *m + 1;
            i__2 = *lwork - iwork + 1;
            dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, info);
            il = iwork;
            dlacpy_((char *)"L", m, m, &a[a_offset], lda, &work[il], &ldwork, (ftnlen)1);
            i__2 = *m - 1;
            i__1 = *m - 1;
            dlaset_((char *)"U", &i__2, &i__1, &c_b46, &c_b46, &work[il + ldwork], &ldwork, (ftnlen)1);
            ie = il + ldwork * *m;
            itauq = ie + *m;
            itaup = itauq + *m;
            iwork = itaup + *m;
            i__2 = *lwork - iwork + 1;
            dgebrd_(m, m, &work[il], &ldwork, &s[1], &work[ie], &work[itauq], &work[itaup],
                    &work[iwork], &i__2, info);
            i__2 = *lwork - iwork + 1;
            dormbr_((char *)"Q", (char *)"L", (char *)"T", m, nrhs, m, &work[il], &ldwork, &work[itauq], &b[b_offset], ldb,
                    &work[iwork], &i__2, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            i__2 = *lwork - iwork + 1;
            dorgbr_((char *)"P", m, m, m, &work[il], &ldwork, &work[itaup], &work[iwork], &i__2, info,
                    (ftnlen)1);
            iwork = ie + *m;
            dbdsqr_((char *)"U", m, m, &c__0, nrhs, &s[1], &work[ie], &work[il], &ldwork, &a[a_offset], lda,
                    &b[b_offset], ldb, &work[iwork], info, (ftnlen)1);
            if (*info != 0) {
                goto L70;
            }
            d__1 = *rcond * s[1];
            thr = max(d__1, sfmin);
            if (*rcond < 0.) {
                d__1 = eps * s[1];
                thr = max(d__1, sfmin);
            }
            *rank = 0;
            i__2 = *m;
            for (i__ = 1; i__ <= i__2; ++i__) {
                if (s[i__] > thr) {
                    drscl_(nrhs, &s[i__], &b[i__ + b_dim1], ldb);
                    ++(*rank);
                } else {
                    dlaset_((char *)"F", &c__1, nrhs, &c_b46, &c_b46, &b[i__ + b_dim1], ldb, (ftnlen)1);
                }
            }
            iwork = ie;
            if (*lwork >= *ldb * *nrhs + iwork - 1 && *nrhs > 1) {
                dgemm_((char *)"T", (char *)"N", m, nrhs, m, &c_b79, &work[il], &ldwork, &b[b_offset], ldb, &c_b46,
                       &work[iwork], ldb, (ftnlen)1, (ftnlen)1);
                dlacpy_((char *)"G", m, nrhs, &work[iwork], ldb, &b[b_offset], ldb, (ftnlen)1);
            } else if (*nrhs > 1) {
                chunk = (*lwork - iwork + 1) / *m;
                i__2 = *nrhs;
                i__1 = chunk;
                for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
                    i__3 = *nrhs - i__ + 1;
                    bl = min(i__3, chunk);
                    dgemm_((char *)"T", (char *)"N", m, &bl, m, &c_b79, &work[il], &ldwork, &b[i__ * b_dim1 + 1],
                           ldb, &c_b46, &work[iwork], m, (ftnlen)1, (ftnlen)1);
                    dlacpy_((char *)"G", m, &bl, &work[iwork], m, &b[i__ * b_dim1 + 1], ldb, (ftnlen)1);
                }
            } else if (*nrhs == 1) {
                dgemv_((char *)"T", m, m, &c_b79, &work[il], &ldwork, &b[b_dim1 + 1], &c__1, &c_b46,
                       &work[iwork], &c__1, (ftnlen)1);
                dcopy_(m, &work[iwork], &c__1, &b[b_dim1 + 1], &c__1);
            }
            i__1 = *n - *m;
            dlaset_((char *)"F", &i__1, nrhs, &c_b46, &c_b46, &b[*m + 1 + b_dim1], ldb, (ftnlen)1);
            iwork = itau + *m;
            i__1 = *lwork - iwork + 1;
            dormlq_((char *)"L", (char *)"T", n, nrhs, m, &a[a_offset], lda, &work[itau], &b[b_offset], ldb,
                    &work[iwork], &i__1, info, (ftnlen)1, (ftnlen)1);
        } else {
            ie = 1;
            itauq = ie + *m;
            itaup = itauq + *m;
            iwork = itaup + *m;
            i__1 = *lwork - iwork + 1;
            dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &work[itaup],
                    &work[iwork], &i__1, info);
            i__1 = *lwork - iwork + 1;
            dormbr_((char *)"Q", (char *)"L", (char *)"T", m, nrhs, n, &a[a_offset], lda, &work[itauq], &b[b_offset], ldb,
                    &work[iwork], &i__1, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            i__1 = *lwork - iwork + 1;
            dorgbr_((char *)"P", m, n, m, &a[a_offset], lda, &work[itaup], &work[iwork], &i__1, info,
                    (ftnlen)1);
            iwork = ie + *m;
            dbdsqr_((char *)"L", m, n, &c__0, nrhs, &s[1], &work[ie], &a[a_offset], lda, dum, &c__1,
                    &b[b_offset], ldb, &work[iwork], info, (ftnlen)1);
            if (*info != 0) {
                goto L70;
            }
            d__1 = *rcond * s[1];
            thr = max(d__1, sfmin);
            if (*rcond < 0.) {
                d__1 = eps * s[1];
                thr = max(d__1, sfmin);
            }
            *rank = 0;
            i__1 = *m;
            for (i__ = 1; i__ <= i__1; ++i__) {
                if (s[i__] > thr) {
                    drscl_(nrhs, &s[i__], &b[i__ + b_dim1], ldb);
                    ++(*rank);
                } else {
                    dlaset_((char *)"F", &c__1, nrhs, &c_b46, &c_b46, &b[i__ + b_dim1], ldb, (ftnlen)1);
                }
            }
            if (*lwork >= *ldb * *nrhs && *nrhs > 1) {
                dgemm_((char *)"T", (char *)"N", n, nrhs, m, &c_b79, &a[a_offset], lda, &b[b_offset], ldb, &c_b46,
                       &work[1], ldb, (ftnlen)1, (ftnlen)1);
                dlacpy_((char *)"F", n, nrhs, &work[1], ldb, &b[b_offset], ldb, (ftnlen)1);
            } else if (*nrhs > 1) {
                chunk = *lwork / *n;
                i__1 = *nrhs;
                i__2 = chunk;
                for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
                    i__3 = *nrhs - i__ + 1;
                    bl = min(i__3, chunk);
                    dgemm_((char *)"T", (char *)"N", n, &bl, m, &c_b79, &a[a_offset], lda, &b[i__ * b_dim1 + 1],
                           ldb, &c_b46, &work[1], n, (ftnlen)1, (ftnlen)1);
                    dlacpy_((char *)"F", n, &bl, &work[1], n, &b[i__ * b_dim1 + 1], ldb, (ftnlen)1);
                }
            } else if (*nrhs == 1) {
                dgemv_((char *)"T", m, n, &c_b79, &a[a_offset], lda, &b[b_offset], &c__1, &c_b46, &work[1],
                       &c__1, (ftnlen)1);
                dcopy_(n, &work[1], &c__1, &b[b_offset], &c__1);
            }
        }
    }
    if (iascl == 1) {
        dlascl_((char *)"G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &b[b_offset], ldb, info, (ftnlen)1);
        dlascl_((char *)"G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &minmn, info, (ftnlen)1);
    } else if (iascl == 2) {
        dlascl_((char *)"G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &b[b_offset], ldb, info, (ftnlen)1);
        dlascl_((char *)"G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &minmn, info, (ftnlen)1);
    }
    if (ibscl == 1) {
        dlascl_((char *)"G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &b[b_offset], ldb, info, (ftnlen)1);
    } else if (ibscl == 2) {
        dlascl_((char *)"G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &b[b_offset], ldb, info, (ftnlen)1);
    }
L70:
    work[1] = (doublereal)maxwrk;
    return 0;
}
#ifdef __cplusplus
}
#endif
