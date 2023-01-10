#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b17 = 1.;
int dsyev_(char *jobz, char *uplo, integer *n, doublereal *a, integer *lda, doublereal *w,
           doublereal *work, integer *lwork, integer *info, ftnlen jobz_len, ftnlen uplo_len)
{
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;
    double sqrt(doublereal);
    integer nb;
    doublereal eps;
    integer inde;
    doublereal anrm;
    integer imax;
    doublereal rmin, rmax;
    extern int dscal_(integer *, doublereal *, doublereal *, integer *);
    doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer iinfo;
    logical lower, wantz;
    extern doublereal dlamch_(char *, ftnlen);
    integer iscale;
    extern int dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *,
                       integer *, doublereal *, integer *, integer *, ftnlen);
    doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    extern int xerbla_(char *, integer *, ftnlen);
    doublereal bignum;
    integer indtau;
    extern int dsterf_(integer *, doublereal *, doublereal *, integer *);
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, integer *, doublereal *,
                              ftnlen, ftnlen);
    integer indwrk;
    extern int dorgtr_(char *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                       integer *, integer *, ftnlen),
        dsteqr_(char *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                doublereal *, integer *, ftnlen),
        dsytrd_(char *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                doublereal *, doublereal *, integer *, integer *, ftnlen);
    integer llwork;
    doublereal smlnum;
    integer lwkopt;
    logical lquery;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;
    --work;
    wantz = lsame_(jobz, (char *)"V", (ftnlen)1, (ftnlen)1);
    lower = lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1);
    lquery = *lwork == -1;
    *info = 0;
    if (!(wantz || lsame_(jobz, (char *)"N", (ftnlen)1, (ftnlen)1))) {
        *info = -1;
    } else if (!(lower || lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1))) {
        *info = -2;
    } else if (*n < 0) {
        *info = -3;
    } else if (*lda < max(1, *n)) {
        *info = -5;
    }
    if (*info == 0) {
        nb = ilaenv_(&c__1, (char *)"DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
        i__1 = 1, i__2 = (nb + 2) * *n;
        lwkopt = max(i__1, i__2);
        work[1] = (doublereal)lwkopt;
        i__1 = 1, i__2 = *n * 3 - 1;
        if (*lwork < max(i__1, i__2) && !lquery) {
            *info = -8;
        }
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DSYEV ", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    if (*n == 1) {
        w[1] = a[a_dim1 + 1];
        work[1] = 2.;
        if (wantz) {
            a[a_dim1 + 1] = 1.;
        }
        return 0;
    }
    safmin = dlamch_((char *)"Safe minimum", (ftnlen)12);
    eps = dlamch_((char *)"Precision", (ftnlen)9);
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = sqrt(smlnum);
    rmax = sqrt(bignum);
    anrm = dlansy_((char *)"M", uplo, n, &a[a_offset], lda, &work[1], (ftnlen)1, (ftnlen)1);
    iscale = 0;
    if (anrm > 0. && anrm < rmin) {
        iscale = 1;
        sigma = rmin / anrm;
    } else if (anrm > rmax) {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if (iscale == 1) {
        dlascl_(uplo, &c__0, &c__0, &c_b17, &sigma, n, n, &a[a_offset], lda, info, (ftnlen)1);
    }
    inde = 1;
    indtau = inde + *n;
    indwrk = indtau + *n;
    llwork = *lwork - indwrk + 1;
    dsytrd_(uplo, n, &a[a_offset], lda, &w[1], &work[inde], &work[indtau], &work[indwrk], &llwork,
            &iinfo, (ftnlen)1);
    if (!wantz) {
        dsterf_(n, &w[1], &work[inde], info);
    } else {
        dorgtr_(uplo, n, &a[a_offset], lda, &work[indtau], &work[indwrk], &llwork, &iinfo,
                (ftnlen)1);
        dsteqr_(jobz, n, &w[1], &work[inde], &a[a_offset], lda, &work[indtau], info, (ftnlen)1);
    }
    if (iscale == 1) {
        if (*info == 0) {
            imax = *n;
        } else {
            imax = *info - 1;
        }
        d__1 = 1. / sigma;
        dscal_(&imax, &d__1, &w[1], &c__1);
    }
    work[1] = (doublereal)lwkopt;
    return 0;
}
#ifdef __cplusplus
}
#endif
