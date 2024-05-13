#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__6 = 6;
static integer c_n1 = -1;
static integer c__9 = 9;
static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b82 = 0.;
int dgelsd_(integer *m, integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *b,
            integer *ldb, doublereal *s, doublereal *rcond, integer *rank, doublereal *work,
            integer *lwork, integer *iwork, integer *info)
{
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;
    double log(doublereal);
    integer ie, il, mm;
    doublereal eps, anrm, bnrm;
    integer itau, nlvl, iascl, ibscl;
    doublereal sfmin;
    integer minmn, maxmn, itaup, itauq, mnthr, nwork;
    extern int dgebrd_(integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                       doublereal *, doublereal *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen),
        dlange_(char *, integer *, integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern int dgelqf_(integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                       integer *, integer *),
        dlalsd_(char *, integer *, integer *, integer *, doublereal *, doublereal *, doublereal *,
                integer *, doublereal *, integer *, doublereal *, integer *, integer *, ftnlen),
        dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *,
                doublereal *, integer *, integer *, ftnlen),
        dgeqrf_(integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                integer *, integer *),
        dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *,
                ftnlen),
        dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                ftnlen),
        xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    doublereal bignum;
    extern int dormbr_(char *, char *, char *, integer *, integer *, integer *, doublereal *,
                       integer *, doublereal *, doublereal *, integer *, doublereal *, integer *,
                       integer *, ftnlen, ftnlen, ftnlen);
    integer wlalsd;
    extern int dormlq_(char *, char *, integer *, integer *, integer *, doublereal *, integer *,
                       doublereal *, doublereal *, integer *, doublereal *, integer *, integer *,
                       ftnlen, ftnlen);
    integer ldwork;
    extern int dormqr_(char *, char *, integer *, integer *, integer *, doublereal *, integer *,
                       doublereal *, doublereal *, integer *, doublereal *, integer *, integer *,
                       ftnlen, ftnlen);
    integer liwork, minwrk, maxwrk;
    doublereal smlnum;
    logical lquery;
    integer smlsiz;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --s;
    --work;
    --iwork;
    *info = 0;
    minmn = min(*m, *n);
    maxmn = max(*m, *n);
    mnthr = ilaenv_(&c__6, (char *)"DGELSD", (char *)" ", m, n, nrhs, &c_n1, (ftnlen)6, (ftnlen)1);
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
    smlsiz = ilaenv_(&c__9, (char *)"DGELSD", (char *)" ", &c__0, &c__0, &c__0, &c__0, (ftnlen)6, (ftnlen)1);
    minwrk = 1;
    liwork = 1;
    minmn = max(1, minmn);
    i__1 = (integer)(log((doublereal)minmn / (doublereal)(smlsiz + 1)) / log(2.)) + 1;
    nlvl = max(i__1, 0);
    if (*info == 0) {
        maxwrk = 0;
        liwork = minmn * 3 * nlvl + minmn * 11;
        mm = *m;
        if (*m >= *n && *m >= mnthr) {
            mm = *n;
            i__1 = maxwrk, i__2 = *n + *n * ilaenv_(&c__1, (char *)"DGEQRF", (char *)" ", m, n, &c_n1, &c_n1,
                                                    (ftnlen)6, (ftnlen)1);
            maxwrk = max(i__1, i__2);
            i__1 = maxwrk, i__2 = *n + *nrhs * ilaenv_(&c__1, (char *)"DORMQR", (char *)"LT", m, nrhs, n, &c_n1,
                                                       (ftnlen)6, (ftnlen)2);
            maxwrk = max(i__1, i__2);
        }
        if (*m >= *n) {
            i__1 = maxwrk, i__2 = *n * 3 + (mm + *n) * ilaenv_(&c__1, (char *)"DGEBRD", (char *)" ", &mm, n, &c_n1,
                                                               &c_n1, (ftnlen)6, (ftnlen)1);
            maxwrk = max(i__1, i__2);
            i__1 = maxwrk, i__2 = *n * 3 + *nrhs * ilaenv_(&c__1, (char *)"DORMBR", (char *)"QLT", &mm, nrhs, n,
                                                           &c_n1, (ftnlen)6, (ftnlen)3);
            maxwrk = max(i__1, i__2);
            i__1 = maxwrk, i__2 = *n * 3 + (*n - 1) * ilaenv_(&c__1, (char *)"DORMBR", (char *)"PLN", n, nrhs, n,
                                                              &c_n1, (ftnlen)6, (ftnlen)3);
            maxwrk = max(i__1, i__2);
            i__1 = smlsiz + 1;
            wlalsd = *n * 9 + (*n << 1) * smlsiz + (*n << 3) * nlvl + *n * *nrhs + i__1 * i__1;
            i__1 = maxwrk, i__2 = *n * 3 + wlalsd;
            maxwrk = max(i__1, i__2);
            i__1 = *n * 3 + mm, i__2 = *n * 3 + *nrhs, i__1 = max(i__1, i__2),
            i__2 = *n * 3 + wlalsd;
            minwrk = max(i__1, i__2);
        }
        if (*n > *m) {
            i__1 = smlsiz + 1;
            wlalsd = *m * 9 + (*m << 1) * smlsiz + (*m << 3) * nlvl + *m * *nrhs + i__1 * i__1;
            if (*n >= mnthr) {
                maxwrk = *m + *m * ilaenv_(&c__1, (char *)"DGELQF", (char *)" ", m, n, &c_n1, &c_n1, (ftnlen)6,
                                           (ftnlen)1);
                i__1 = maxwrk, i__2 = *m * *m + (*m << 2) +
                                      (*m << 1) * ilaenv_(&c__1, (char *)"DGEBRD", (char *)" ", m, m, &c_n1, &c_n1,
                                                          (ftnlen)6, (ftnlen)1);
                maxwrk = max(i__1, i__2);
                i__1 = maxwrk, i__2 = *m * *m + (*m << 2) +
                                      *nrhs * ilaenv_(&c__1, (char *)"DORMBR", (char *)"QLT", m, nrhs, m, &c_n1,
                                                      (ftnlen)6, (ftnlen)3);
                maxwrk = max(i__1, i__2);
                i__1 = maxwrk, i__2 = *m * *m + (*m << 2) +
                                      (*m - 1) * ilaenv_(&c__1, (char *)"DORMBR", (char *)"PLN", m, nrhs, m, &c_n1,
                                                         (ftnlen)6, (ftnlen)3);
                maxwrk = max(i__1, i__2);
                if (*nrhs > 1) {
                    i__1 = maxwrk, i__2 = *m * *m + *m + *m * *nrhs;
                    maxwrk = max(i__1, i__2);
                } else {
                    i__1 = maxwrk, i__2 = *m * *m + (*m << 1);
                    maxwrk = max(i__1, i__2);
                }
                i__1 = maxwrk, i__2 = *m + *nrhs * ilaenv_(&c__1, (char *)"DORMLQ", (char *)"LT", n, nrhs, m, &c_n1,
                                                           (ftnlen)6, (ftnlen)2);
                maxwrk = max(i__1, i__2);
                i__1 = maxwrk, i__2 = *m * *m + (*m << 2) + wlalsd;
                maxwrk = max(i__1, i__2);
                i__3 = *m, i__4 = (*m << 1) - 4, i__3 = max(i__3, i__4), i__3 = max(i__3, *nrhs),
                i__4 = *n - *m * 3;
                i__1 = maxwrk, i__2 = (*m << 2) + *m * *m + max(i__3, i__4);
                maxwrk = max(i__1, i__2);
            } else {
                maxwrk = *m * 3 + (*n + *m) * ilaenv_(&c__1, (char *)"DGEBRD", (char *)" ", m, n, &c_n1, &c_n1,
                                                      (ftnlen)6, (ftnlen)1);
                i__1 = maxwrk, i__2 = *m * 3 + *nrhs * ilaenv_(&c__1, (char *)"DORMBR", (char *)"QLT", m, nrhs, n,
                                                               &c_n1, (ftnlen)6, (ftnlen)3);
                maxwrk = max(i__1, i__2);
                i__1 = maxwrk, i__2 = *m * 3 + *m * ilaenv_(&c__1, (char *)"DORMBR", (char *)"PLN", n, nrhs, m,
                                                            &c_n1, (ftnlen)6, (ftnlen)3);
                maxwrk = max(i__1, i__2);
                i__1 = maxwrk, i__2 = *m * 3 + wlalsd;
                maxwrk = max(i__1, i__2);
            }
            i__1 = *m * 3 + *nrhs, i__2 = *m * 3 + *m, i__1 = max(i__1, i__2),
            i__2 = *m * 3 + wlalsd;
            minwrk = max(i__1, i__2);
        }
        minwrk = min(minwrk, maxwrk);
        work[1] = (doublereal)maxwrk;
        iwork[1] = liwork;
        if (*lwork < minwrk && !lquery) {
            *info = -12;
        }
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DGELSD", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        goto L10;
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
        dlaset_((char *)"F", &i__1, nrhs, &c_b82, &c_b82, &b[b_offset], ldb, (ftnlen)1);
        dlaset_((char *)"F", &minmn, &c__1, &c_b82, &c_b82, &s[1], &c__1, (ftnlen)1);
        *rank = 0;
        goto L10;
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
    if (*m < *n) {
        i__1 = *n - *m;
        dlaset_((char *)"F", &i__1, nrhs, &c_b82, &c_b82, &b[*m + 1 + b_dim1], ldb, (ftnlen)1);
    }
    if (*m >= *n) {
        mm = *m;
        if (*m >= mnthr) {
            mm = *n;
            itau = 1;
            nwork = itau + *n;
            i__1 = *lwork - nwork + 1;
            dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &i__1, info);
            i__1 = *lwork - nwork + 1;
            dormqr_((char *)"L", (char *)"T", m, nrhs, n, &a[a_offset], lda, &work[itau], &b[b_offset], ldb,
                    &work[nwork], &i__1, info, (ftnlen)1, (ftnlen)1);
            if (*n > 1) {
                i__1 = *n - 1;
                i__2 = *n - 1;
                dlaset_((char *)"L", &i__1, &i__2, &c_b82, &c_b82, &a[a_dim1 + 2], lda, (ftnlen)1);
            }
        }
        ie = 1;
        itauq = ie + *n;
        itaup = itauq + *n;
        nwork = itaup + *n;
        i__1 = *lwork - nwork + 1;
        dgebrd_(&mm, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &work[itaup],
                &work[nwork], &i__1, info);
        i__1 = *lwork - nwork + 1;
        dormbr_((char *)"Q", (char *)"L", (char *)"T", &mm, nrhs, n, &a[a_offset], lda, &work[itauq], &b[b_offset], ldb,
                &work[nwork], &i__1, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
        dlalsd_((char *)"U", &smlsiz, n, nrhs, &s[1], &work[ie], &b[b_offset], ldb, rcond, rank,
                &work[nwork], &iwork[1], info, (ftnlen)1);
        if (*info != 0) {
            goto L10;
        }
        i__1 = *lwork - nwork + 1;
        dormbr_((char *)"P", (char *)"L", (char *)"N", n, nrhs, n, &a[a_offset], lda, &work[itaup], &b[b_offset], ldb,
                &work[nwork], &i__1, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
    } else {
        i__1 = *m, i__2 = (*m << 1) - 4, i__1 = max(i__1, i__2), i__1 = max(i__1, *nrhs),
        i__2 = *n - *m * 3, i__1 = max(i__1, i__2);
        if (*n >= mnthr && *lwork >= (*m << 2) + *m * *m + max(i__1, wlalsd)) {
            ldwork = *m;
            i__3 = *m, i__4 = (*m << 1) - 4, i__3 = max(i__3, i__4), i__3 = max(i__3, *nrhs),
            i__4 = *n - *m * 3;
            i__1 = (*m << 2) + *m * *lda + max(i__3, i__4), i__2 = *m * *lda + *m + *m * *nrhs,
            i__1 = max(i__1, i__2), i__2 = (*m << 2) + *m * *lda + wlalsd;
            if (*lwork >= max(i__1, i__2)) {
                ldwork = *lda;
            }
            itau = 1;
            nwork = *m + 1;
            i__1 = *lwork - nwork + 1;
            dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[nwork], &i__1, info);
            il = nwork;
            dlacpy_((char *)"L", m, m, &a[a_offset], lda, &work[il], &ldwork, (ftnlen)1);
            i__1 = *m - 1;
            i__2 = *m - 1;
            dlaset_((char *)"U", &i__1, &i__2, &c_b82, &c_b82, &work[il + ldwork], &ldwork, (ftnlen)1);
            ie = il + ldwork * *m;
            itauq = ie + *m;
            itaup = itauq + *m;
            nwork = itaup + *m;
            i__1 = *lwork - nwork + 1;
            dgebrd_(m, m, &work[il], &ldwork, &s[1], &work[ie], &work[itauq], &work[itaup],
                    &work[nwork], &i__1, info);
            i__1 = *lwork - nwork + 1;
            dormbr_((char *)"Q", (char *)"L", (char *)"T", m, nrhs, m, &work[il], &ldwork, &work[itauq], &b[b_offset], ldb,
                    &work[nwork], &i__1, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            dlalsd_((char *)"U", &smlsiz, m, nrhs, &s[1], &work[ie], &b[b_offset], ldb, rcond, rank,
                    &work[nwork], &iwork[1], info, (ftnlen)1);
            if (*info != 0) {
                goto L10;
            }
            i__1 = *lwork - nwork + 1;
            dormbr_((char *)"P", (char *)"L", (char *)"N", m, nrhs, m, &work[il], &ldwork, &work[itaup], &b[b_offset], ldb,
                    &work[nwork], &i__1, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            i__1 = *n - *m;
            dlaset_((char *)"F", &i__1, nrhs, &c_b82, &c_b82, &b[*m + 1 + b_dim1], ldb, (ftnlen)1);
            nwork = itau + *m;
            i__1 = *lwork - nwork + 1;
            dormlq_((char *)"L", (char *)"T", n, nrhs, m, &a[a_offset], lda, &work[itau], &b[b_offset], ldb,
                    &work[nwork], &i__1, info, (ftnlen)1, (ftnlen)1);
        } else {
            ie = 1;
            itauq = ie + *m;
            itaup = itauq + *m;
            nwork = itaup + *m;
            i__1 = *lwork - nwork + 1;
            dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &work[itaup],
                    &work[nwork], &i__1, info);
            i__1 = *lwork - nwork + 1;
            dormbr_((char *)"Q", (char *)"L", (char *)"T", m, nrhs, n, &a[a_offset], lda, &work[itauq], &b[b_offset], ldb,
                    &work[nwork], &i__1, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            dlalsd_((char *)"L", &smlsiz, m, nrhs, &s[1], &work[ie], &b[b_offset], ldb, rcond, rank,
                    &work[nwork], &iwork[1], info, (ftnlen)1);
            if (*info != 0) {
                goto L10;
            }
            i__1 = *lwork - nwork + 1;
            dormbr_((char *)"P", (char *)"L", (char *)"N", n, nrhs, m, &a[a_offset], lda, &work[itaup], &b[b_offset], ldb,
                    &work[nwork], &i__1, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
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
L10:
    work[1] = (doublereal)maxwrk;
    iwork[1] = liwork;
    return 0;
}
#ifdef __cplusplus
}
#endif
