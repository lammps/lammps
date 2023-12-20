#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__9 = 9;
static integer c__0 = 0;
static integer c__2 = 2;
static doublereal c_b17 = 0.;
static doublereal c_b18 = 1.;
static integer c__1 = 1;
int zstedc_(char *compz, integer *n, doublereal *d__, doublereal *e, doublecomplex *z__,
            integer *ldz, doublecomplex *work, integer *lwork, doublereal *rwork, integer *lrwork,
            integer *iwork, integer *liwork, integer *info, ftnlen compz_len)
{
    integer z_dim1, z_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    double log(doublereal);
    integer pow_lmp_ii(integer *, integer *);
    double sqrt(doublereal);
    integer i__, j, k, m;
    doublereal p;
    integer ii, ll, lgn;
    doublereal eps, tiny;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer lwmin, start;
    extern int zswap_(integer *, doublecomplex *, integer *, doublecomplex *, integer *),
        zlaed0_(integer *, integer *, doublereal *, doublereal *, doublecomplex *, integer *,
                doublecomplex *, integer *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern int dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *,
                       integer *, doublereal *, integer *, integer *, ftnlen),
        dstedc_(char *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                doublereal *, integer *, integer *, integer *, integer *, ftnlen),
        dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                ftnlen),
        xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    integer finish;
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *, ftnlen);
    extern int dsterf_(integer *, doublereal *, doublereal *, integer *),
        zlacrm_(integer *, integer *, doublecomplex *, integer *, doublereal *, integer *,
                doublecomplex *, integer *, doublereal *);
    integer liwmin, icompz;
    extern int dsteqr_(char *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                       doublereal *, integer *, ftnlen),
        zlacpy_(char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *,
                integer *, ftnlen);
    doublereal orgnrm;
    integer lrwmin;
    logical lquery;
    integer smlsiz;
    extern int zsteqr_(char *, integer *, doublereal *, doublereal *, doublecomplex *, integer *,
                       doublereal *, integer *, ftnlen);
    --d__;
    --e;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --rwork;
    --iwork;
    *info = 0;
    lquery = *lwork == -1 || *lrwork == -1 || *liwork == -1;
    if (lsame_(compz, (char *)"N", (ftnlen)1, (ftnlen)1)) {
        icompz = 0;
    } else if (lsame_(compz, (char *)"V", (ftnlen)1, (ftnlen)1)) {
        icompz = 1;
    } else if (lsame_(compz, (char *)"I", (ftnlen)1, (ftnlen)1)) {
        icompz = 2;
    } else {
        icompz = -1;
    }
    if (icompz < 0) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*ldz < 1 || icompz > 0 && *ldz < max(1, *n)) {
        *info = -6;
    }
    if (*info == 0) {
        smlsiz = ilaenv_(&c__9, (char *)"ZSTEDC", (char *)" ", &c__0, &c__0, &c__0, &c__0, (ftnlen)6, (ftnlen)1);
        if (*n <= 1 || icompz == 0) {
            lwmin = 1;
            liwmin = 1;
            lrwmin = 1;
        } else if (*n <= smlsiz) {
            lwmin = 1;
            liwmin = 1;
            lrwmin = *n - 1 << 1;
        } else if (icompz == 1) {
            lgn = (integer)(log((doublereal)(*n)) / log(2.));
            if (pow_lmp_ii(&c__2, &lgn) < *n) {
                ++lgn;
            }
            if (pow_lmp_ii(&c__2, &lgn) < *n) {
                ++lgn;
            }
            lwmin = *n * *n;
            i__1 = *n;
            lrwmin = *n * 3 + 1 + (*n << 1) * lgn + (i__1 * i__1 << 2);
            liwmin = *n * 6 + 6 + *n * 5 * lgn;
        } else if (icompz == 2) {
            lwmin = 1;
            i__1 = *n;
            lrwmin = (*n << 2) + 1 + (i__1 * i__1 << 1);
            liwmin = *n * 5 + 3;
        }
        work[1].r = (doublereal)lwmin, work[1].i = 0.;
        rwork[1] = (doublereal)lrwmin;
        iwork[1] = liwmin;
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
        xerbla_((char *)"ZSTEDC", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    if (*n == 1) {
        if (icompz != 0) {
            i__1 = z_dim1 + 1;
            z__[i__1].r = 1., z__[i__1].i = 0.;
        }
        return 0;
    }
    if (icompz == 0) {
        dsterf_(n, &d__[1], &e[1], info);
        goto L70;
    }
    if (*n <= smlsiz) {
        zsteqr_(compz, n, &d__[1], &e[1], &z__[z_offset], ldz, &rwork[1], info, (ftnlen)1);
    } else {
        if (icompz == 2) {
            dlaset_((char *)"Full", n, n, &c_b17, &c_b18, &rwork[1], n, (ftnlen)4);
            ll = *n * *n + 1;
            i__1 = *lrwork - ll + 1;
            dstedc_((char *)"I", n, &d__[1], &e[1], &rwork[1], n, &rwork[ll], &i__1, &iwork[1], liwork,
                    info, (ftnlen)1);
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *n;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    i__3 = i__ + j * z_dim1;
                    i__4 = (j - 1) * *n + i__;
                    z__[i__3].r = rwork[i__4], z__[i__3].i = 0.;
                }
            }
            goto L70;
        }
        orgnrm = dlanst_((char *)"M", n, &d__[1], &e[1], (ftnlen)1);
        if (orgnrm == 0.) {
            goto L70;
        }
        eps = dlamch_((char *)"Epsilon", (ftnlen)7);
        start = 1;
    L30:
        if (start <= *n) {
            finish = start;
        L40:
            if (finish < *n) {
                tiny = eps * sqrt((d__1 = d__[finish], abs(d__1))) *
                       sqrt((d__2 = d__[finish + 1], abs(d__2)));
                if ((d__1 = e[finish], abs(d__1)) > tiny) {
                    ++finish;
                    goto L40;
                }
            }
            m = finish - start + 1;
            if (m > smlsiz) {
                orgnrm = dlanst_((char *)"M", &m, &d__[start], &e[start], (ftnlen)1);
                dlascl_((char *)"G", &c__0, &c__0, &orgnrm, &c_b18, &m, &c__1, &d__[start], &m, info,
                        (ftnlen)1);
                i__1 = m - 1;
                i__2 = m - 1;
                dlascl_((char *)"G", &c__0, &c__0, &orgnrm, &c_b18, &i__1, &c__1, &e[start], &i__2, info,
                        (ftnlen)1);
                zlaed0_(n, &m, &d__[start], &e[start], &z__[start * z_dim1 + 1], ldz, &work[1], n,
                        &rwork[1], &iwork[1], info);
                if (*info > 0) {
                    *info = (*info / (m + 1) + start - 1) * (*n + 1) + *info % (m + 1) + start - 1;
                    goto L70;
                }
                dlascl_((char *)"G", &c__0, &c__0, &c_b18, &orgnrm, &m, &c__1, &d__[start], &m, info,
                        (ftnlen)1);
            } else {
                dsteqr_((char *)"I", &m, &d__[start], &e[start], &rwork[1], &m, &rwork[m * m + 1], info,
                        (ftnlen)1);
                zlacrm_(n, &m, &z__[start * z_dim1 + 1], ldz, &rwork[1], &m, &work[1], n,
                        &rwork[m * m + 1]);
                zlacpy_((char *)"A", n, &m, &work[1], n, &z__[start * z_dim1 + 1], ldz, (ftnlen)1);
                if (*info > 0) {
                    *info = start * (*n + 1) + finish;
                    goto L70;
                }
            }
            start = finish + 1;
            goto L30;
        }
        i__1 = *n;
        for (ii = 2; ii <= i__1; ++ii) {
            i__ = ii - 1;
            k = i__;
            p = d__[i__];
            i__2 = *n;
            for (j = ii; j <= i__2; ++j) {
                if (d__[j] < p) {
                    k = j;
                    p = d__[j];
                }
            }
            if (k != i__) {
                d__[k] = d__[i__];
                d__[i__] = p;
                zswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[k * z_dim1 + 1], &c__1);
            }
        }
    }
L70:
    work[1].r = (doublereal)lwmin, work[1].i = 0.;
    rwork[1] = (doublereal)lrwmin;
    iwork[1] = liwmin;
    return 0;
}
#ifdef __cplusplus
}
#endif
