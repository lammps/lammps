#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b18 = 1.;
int zheevd_(char *jobz, char *uplo, integer *n, doublecomplex *a, integer *lda, doublereal *w,
            doublecomplex *work, integer *lwork, doublereal *rwork, integer *lrwork, integer *iwork,
            integer *liwork, integer *info, ftnlen jobz_len, ftnlen uplo_len)
{
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;
    double sqrt(doublereal);
    doublereal eps;
    integer inde;
    doublereal anrm;
    integer imax;
    doublereal rmin, rmax;
    integer lopt;
    extern int dscal_(integer *, doublereal *, doublereal *, integer *);
    doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer iinfo, lwmin, liopt;
    logical lower;
    integer llrwk, lropt;
    logical wantz;
    integer indwk2, llwrk2;
    extern doublereal dlamch_(char *, ftnlen);
    integer iscale;
    doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    extern int xerbla_(char *, integer *, ftnlen);
    doublereal bignum;
    extern doublereal zlanhe_(char *, char *, integer *, doublecomplex *, integer *, doublereal *,
                              ftnlen, ftnlen);
    integer indtau;
    extern int dsterf_(integer *, doublereal *, doublereal *, integer *),
        zlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *,
                doublecomplex *, integer *, integer *, ftnlen),
        zstedc_(char *, integer *, doublereal *, doublereal *, doublecomplex *, integer *,
                doublecomplex *, integer *, doublereal *, integer *, integer *, integer *,
                integer *, ftnlen);
    integer indrwk, indwrk, liwmin;
    extern int zhetrd_(char *, integer *, doublecomplex *, integer *, doublereal *, doublereal *,
                       doublecomplex *, doublecomplex *, integer *, integer *, ftnlen),
        zlacpy_(char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *,
                integer *, ftnlen);
    integer lrwmin, llwork;
    doublereal smlnum;
    logical lquery;
    extern int zunmtr_(char *, char *, char *, integer *, integer *, doublecomplex *, integer *,
                       doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *,
                       integer *, ftnlen, ftnlen, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;
    --work;
    --rwork;
    --iwork;
    wantz = lsame_(jobz, (char *)"V", (ftnlen)1, (ftnlen)1);
    lower = lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1);
    lquery = *lwork == -1 || *lrwork == -1 || *liwork == -1;
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
        if (*n <= 1) {
            lwmin = 1;
            lrwmin = 1;
            liwmin = 1;
            lopt = lwmin;
            lropt = lrwmin;
            liopt = liwmin;
        } else {
            if (wantz) {
                lwmin = (*n << 1) + *n * *n;
                i__1 = *n;
                lrwmin = *n * 5 + 1 + (i__1 * i__1 << 1);
                liwmin = *n * 5 + 3;
            } else {
                lwmin = *n + 1;
                lrwmin = *n;
                liwmin = 1;
            }
            i__1 = lwmin, i__2 = *n + *n * ilaenv_(&c__1, (char *)"ZHETRD", uplo, n, &c_n1, &c_n1, &c_n1,
                                                   (ftnlen)6, (ftnlen)1);
            lopt = max(i__1, i__2);
            lropt = lrwmin;
            liopt = liwmin;
        }
        work[1].r = (doublereal)lopt, work[1].i = 0.;
        rwork[1] = (doublereal)lropt;
        iwork[1] = liopt;
        if (*lwork < lwmin && !lquery) {
            *info = -8;
        } else if (*lrwork < lrwmin && !lquery) {
            *info = -10;
        } else if (*liwork < liwmin && !lquery) {
            *info = -12;
        }
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"ZHEEVD", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    if (*n == 1) {
        i__1 = a_dim1 + 1;
        w[1] = a[i__1].r;
        if (wantz) {
            i__1 = a_dim1 + 1;
            a[i__1].r = 1., a[i__1].i = 0.;
        }
        return 0;
    }
    safmin = dlamch_((char *)"Safe minimum", (ftnlen)12);
    eps = dlamch_((char *)"Precision", (ftnlen)9);
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = sqrt(smlnum);
    rmax = sqrt(bignum);
    anrm = zlanhe_((char *)"M", uplo, n, &a[a_offset], lda, &rwork[1], (ftnlen)1, (ftnlen)1);
    iscale = 0;
    if (anrm > 0. && anrm < rmin) {
        iscale = 1;
        sigma = rmin / anrm;
    } else if (anrm > rmax) {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if (iscale == 1) {
        zlascl_(uplo, &c__0, &c__0, &c_b18, &sigma, n, n, &a[a_offset], lda, info, (ftnlen)1);
    }
    inde = 1;
    indtau = 1;
    indwrk = indtau + *n;
    indrwk = inde + *n;
    indwk2 = indwrk + *n * *n;
    llwork = *lwork - indwrk + 1;
    llwrk2 = *lwork - indwk2 + 1;
    llrwk = *lrwork - indrwk + 1;
    zhetrd_(uplo, n, &a[a_offset], lda, &w[1], &rwork[inde], &work[indtau], &work[indwrk], &llwork,
            &iinfo, (ftnlen)1);
    if (!wantz) {
        dsterf_(n, &w[1], &rwork[inde], info);
    } else {
        zstedc_((char *)"I", n, &w[1], &rwork[inde], &work[indwrk], n, &work[indwk2], &llwrk2,
                &rwork[indrwk], &llrwk, &iwork[1], liwork, info, (ftnlen)1);
        zunmtr_((char *)"L", uplo, (char *)"N", n, n, &a[a_offset], lda, &work[indtau], &work[indwrk], n,
                &work[indwk2], &llwrk2, &iinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1);
        zlacpy_((char *)"A", n, n, &work[indwrk], n, &a[a_offset], lda, (ftnlen)1);
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
    work[1].r = (doublereal)lopt, work[1].i = 0.;
    rwork[1] = (doublereal)lropt;
    iwork[1] = liopt;
    return 0;
}
#ifdef __cplusplus
}
#endif
