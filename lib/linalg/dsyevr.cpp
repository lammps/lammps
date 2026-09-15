#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__10 = 10;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c_n1 = -1;
int dsyevr_(char *jobz, char *range, char *uplo, integer *n, doublereal *a, integer *lda,
            doublereal *vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol,
            integer *m, doublereal *w, doublereal *z__, integer *ldz, integer *isuppz,
            doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer *info,
            ftnlen jobz_len, ftnlen range_len, ftnlen uplo_len)
{
    integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;
    double sqrt(doublereal);
    integer i__, j, nb, jj;
    doublereal eps, vll, vuu, tmp1;
    integer indd, inde;
    doublereal anrm;
    integer imax;
    doublereal rmin, rmax;
    integer inddd, indee;
    extern int dscal_(integer *, doublereal *, doublereal *, integer *);
    doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer iinfo;
    char order[1];
    integer indwk;
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *),
        dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer lwmin;
    logical lower, wantz;
    extern doublereal dlamch_(char *, ftnlen);
    logical alleig, indeig;
    integer iscale, ieeeok, indibl, indifl;
    logical valeig;
    doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    extern int xerbla_(char *, integer *, ftnlen);
    doublereal abstll, bignum;
    integer indtau, indisp;
    extern int dstein_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *,
                       integer *, doublereal *, integer *, doublereal *, integer *, integer *,
                       integer *),
        dsterf_(integer *, doublereal *, doublereal *, integer *);
    integer indiwo, indwkn;
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, integer *, doublereal *,
                              ftnlen, ftnlen);
    extern int dstebz_(char *, char *, integer *, doublereal *, doublereal *, integer *, integer *,
                       doublereal *, doublereal *, doublereal *, integer *, integer *, doublereal *,
                       integer *, integer *, doublereal *, integer *, integer *, ftnlen, ftnlen),
        dstemr_(char *, char *, integer *, doublereal *, doublereal *, doublereal *, doublereal *,
                integer *, integer *, integer *, doublereal *, doublereal *, integer *, integer *,
                integer *, logical *, doublereal *, integer *, integer *, integer *, integer *,
                ftnlen, ftnlen);
    integer liwmin;
    logical tryrac;
    extern int dormtr_(char *, char *, char *, integer *, integer *, doublereal *, integer *,
                       doublereal *, doublereal *, integer *, doublereal *, integer *, integer *,
                       ftnlen, ftnlen, ftnlen);
    integer llwrkn, llwork, nsplit;
    doublereal smlnum;
    extern int dsytrd_(char *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                       doublereal *, doublereal *, integer *, integer *, ftnlen);
    integer lwkopt;
    logical lquery;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --isuppz;
    --work;
    --iwork;
    ieeeok = ilaenv_(&c__10, (char *)"DSYEVR", (char *)"N", &c__1, &c__2, &c__3, &c__4, (ftnlen)6, (ftnlen)1);
    lower = lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1);
    wantz = lsame_(jobz, (char *)"V", (ftnlen)1, (ftnlen)1);
    alleig = lsame_(range, (char *)"A", (ftnlen)1, (ftnlen)1);
    valeig = lsame_(range, (char *)"V", (ftnlen)1, (ftnlen)1);
    indeig = lsame_(range, (char *)"I", (ftnlen)1, (ftnlen)1);
    lquery = *lwork == -1 || *liwork == -1;
    if (*n <= 1) {
        lwmin = 1;
        liwmin = 1;
    } else {
        lwmin = *n * 26;
        liwmin = *n * 10;
    }
    *info = 0;
    if (!(wantz || lsame_(jobz, (char *)"N", (ftnlen)1, (ftnlen)1))) {
        *info = -1;
    } else if (!(alleig || valeig || indeig)) {
        *info = -2;
    } else if (!(lower || lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1))) {
        *info = -3;
    } else if (*n < 0) {
        *info = -4;
    } else if (*lda < max(1, *n)) {
        *info = -6;
    } else {
        if (valeig) {
            if (*n > 0 && *vu <= *vl) {
                *info = -8;
            }
        } else if (indeig) {
            if (*il < 1 || *il > max(1, *n)) {
                *info = -9;
            } else if (*iu < min(*n, *il) || *iu > *n) {
                *info = -10;
            }
        }
    }
    if (*info == 0) {
        if (*ldz < 1 || wantz && *ldz < *n) {
            *info = -15;
        } else if (*lwork < lwmin && !lquery) {
            *info = -18;
        } else if (*liwork < liwmin && !lquery) {
            *info = -20;
        }
    }
    if (*info == 0) {
        nb = ilaenv_(&c__1, (char *)"DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
        i__1 = nb,
        i__2 = ilaenv_(&c__1, (char *)"DORMTR", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
        nb = max(i__1, i__2);
        i__1 = (nb + 1) * *n;
        lwkopt = max(i__1, lwmin);
        work[1] = (doublereal)lwkopt;
        iwork[1] = liwmin;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DSYEVR", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    *m = 0;
    if (*n == 0) {
        work[1] = 1.;
        return 0;
    }
    if (*n == 1) {
        work[1] = 1.;
        if (alleig || indeig) {
            *m = 1;
            w[1] = a[a_dim1 + 1];
        } else {
            if (*vl < a[a_dim1 + 1] && *vu >= a[a_dim1 + 1]) {
                *m = 1;
                w[1] = a[a_dim1 + 1];
            }
        }
        if (wantz) {
            z__[z_dim1 + 1] = 1.;
            isuppz[1] = 1;
            isuppz[2] = 1;
        }
        return 0;
    }
    safmin = dlamch_((char *)"Safe minimum", (ftnlen)12);
    eps = dlamch_((char *)"Precision", (ftnlen)9);
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = sqrt(smlnum);
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
    rmax = min(d__1, d__2);
    iscale = 0;
    abstll = *abstol;
    if (valeig) {
        vll = *vl;
        vuu = *vu;
    }
    anrm = dlansy_((char *)"M", uplo, n, &a[a_offset], lda, &work[1], (ftnlen)1, (ftnlen)1);
    if (anrm > 0. && anrm < rmin) {
        iscale = 1;
        sigma = rmin / anrm;
    } else if (anrm > rmax) {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if (iscale == 1) {
        if (lower) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *n - j + 1;
                dscal_(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                dscal_(&j, &sigma, &a[j * a_dim1 + 1], &c__1);
            }
        }
        if (*abstol > 0.) {
            abstll = *abstol * sigma;
        }
        if (valeig) {
            vll = *vl * sigma;
            vuu = *vu * sigma;
        }
    }
    indtau = 1;
    indd = indtau + *n;
    inde = indd + *n;
    inddd = inde + *n;
    indee = inddd + *n;
    indwk = indee + *n;
    llwork = *lwork - indwk + 1;
    indibl = 1;
    indisp = indibl + *n;
    indifl = indisp + *n;
    indiwo = indifl + *n;
    dsytrd_(uplo, n, &a[a_offset], lda, &work[indd], &work[inde], &work[indtau], &work[indwk],
            &llwork, &iinfo, (ftnlen)1);
    if ((alleig || indeig && *il == 1 && *iu == *n) && ieeeok == 1) {
        if (!wantz) {
            dcopy_(n, &work[indd], &c__1, &w[1], &c__1);
            i__1 = *n - 1;
            dcopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
            dsterf_(n, &w[1], &work[indee], info);
        } else {
            i__1 = *n - 1;
            dcopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
            dcopy_(n, &work[indd], &c__1, &work[inddd], &c__1);
            if (*abstol <= *n * 2. * eps) {
                tryrac = TRUE_;
            } else {
                tryrac = FALSE_;
            }
            dstemr_(jobz, (char *)"A", n, &work[inddd], &work[indee], vl, vu, il, iu, m, &w[1],
                    &z__[z_offset], ldz, n, &isuppz[1], &tryrac, &work[indwk], lwork, &iwork[1],
                    liwork, info, (ftnlen)1, (ftnlen)1);
            if (wantz && *info == 0) {
                indwkn = inde;
                llwrkn = *lwork - indwkn + 1;
                dormtr_((char *)"L", uplo, (char *)"N", n, m, &a[a_offset], lda, &work[indtau], &z__[z_offset], ldz,
                        &work[indwkn], &llwrkn, &iinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            }
        }
        if (*info == 0) {
            *m = *n;
            goto L30;
        }
        *info = 0;
    }
    if (wantz) {
        *(unsigned char *)order = 'B';
    } else {
        *(unsigned char *)order = 'E';
    }
    dstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &work[indd], &work[inde], m, &nsplit,
            &w[1], &iwork[indibl], &iwork[indisp], &work[indwk], &iwork[indiwo], info, (ftnlen)1,
            (ftnlen)1);
    if (wantz) {
        dstein_(n, &work[indd], &work[inde], m, &w[1], &iwork[indibl], &iwork[indisp],
                &z__[z_offset], ldz, &work[indwk], &iwork[indiwo], &iwork[indifl], info);
        indwkn = inde;
        llwrkn = *lwork - indwkn + 1;
        dormtr_((char *)"L", uplo, (char *)"N", n, m, &a[a_offset], lda, &work[indtau], &z__[z_offset], ldz,
                &work[indwkn], &llwrkn, &iinfo, (ftnlen)1, (ftnlen)1, (ftnlen)1);
    }
L30:
    if (iscale == 1) {
        if (*info == 0) {
            imax = *m;
        } else {
            imax = *info - 1;
        }
        d__1 = 1. / sigma;
        dscal_(&imax, &d__1, &w[1], &c__1);
    }
    if (wantz) {
        i__1 = *m - 1;
        for (j = 1; j <= i__1; ++j) {
            i__ = 0;
            tmp1 = w[j];
            i__2 = *m;
            for (jj = j + 1; jj <= i__2; ++jj) {
                if (w[jj] < tmp1) {
                    i__ = jj;
                    tmp1 = w[jj];
                }
            }
            if (i__ != 0) {
                w[i__] = w[j];
                w[j] = tmp1;
                dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1], &c__1);
            }
        }
    }
    work[1] = (doublereal)lwkopt;
    iwork[1] = liwmin;
    return 0;
}
#ifdef __cplusplus
}
#endif
