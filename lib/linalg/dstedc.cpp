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
int dstedc_(char *compz, integer *n, doublereal *d__, doublereal *e, doublereal *z__, integer *ldz,
            doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer *info,
            ftnlen compz_len)
{
    integer z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;
    double log(doublereal);
    integer pow_lmp_ii(integer *, integer *);
    double sqrt(doublereal);
    integer i__, j, k, m;
    doublereal p;
    integer ii, lgn;
    doublereal eps, tiny;
    extern int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
                      integer *, doublereal *, integer *, doublereal *, doublereal *, integer *,
                      ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer lwmin;
    extern int dlaed0_(integer *, integer *, integer *, doublereal *, doublereal *, doublereal *,
                       integer *, doublereal *, integer *, doublereal *, integer *, integer *);
    integer start;
    extern doublereal dlamch_(char *, ftnlen);
    extern int dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *,
                       integer *, doublereal *, integer *, integer *, ftnlen),
        dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *,
                ftnlen),
        dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    extern int xerbla_(char *, integer *, ftnlen);
    integer finish;
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *, ftnlen);
    extern int dsterf_(integer *, doublereal *, doublereal *, integer *),
        dlasrt_(char *, integer *, doublereal *, integer *, ftnlen);
    integer liwmin, icompz;
    extern int dsteqr_(char *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                       doublereal *, integer *, ftnlen);
    doublereal orgnrm;
    logical lquery;
    integer smlsiz, storez, strtrw;
    --d__;
    --e;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --iwork;
    *info = 0;
    lquery = *lwork == -1 || *liwork == -1;
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
        smlsiz = ilaenv_(&c__9, (char *)"DSTEDC", (char *)" ", &c__0, &c__0, &c__0, &c__0, (ftnlen)6, (ftnlen)1);
        if (*n <= 1 || icompz == 0) {
            liwmin = 1;
            lwmin = 1;
        } else if (*n <= smlsiz) {
            liwmin = 1;
            lwmin = *n - 1 << 1;
        } else {
            lgn = (integer)(log((doublereal)(*n)) / log(2.));
            if (pow_lmp_ii(&c__2, &lgn) < *n) {
                ++lgn;
            }
            if (pow_lmp_ii(&c__2, &lgn) < *n) {
                ++lgn;
            }
            if (icompz == 1) {
                i__1 = *n;
                lwmin = *n * 3 + 1 + (*n << 1) * lgn + (i__1 * i__1 << 2);
                liwmin = *n * 6 + 6 + *n * 5 * lgn;
            } else if (icompz == 2) {
                i__1 = *n;
                lwmin = (*n << 2) + 1 + i__1 * i__1;
                liwmin = *n * 5 + 3;
            }
        }
        work[1] = (doublereal)lwmin;
        iwork[1] = liwmin;
        if (*lwork < lwmin && !lquery) {
            *info = -8;
        } else if (*liwork < liwmin && !lquery) {
            *info = -10;
        }
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DSTEDC", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    if (*n == 1) {
        if (icompz != 0) {
            z__[z_dim1 + 1] = 1.;
        }
        return 0;
    }
    if (icompz == 0) {
        dsterf_(n, &d__[1], &e[1], info);
        goto L50;
    }
    if (*n <= smlsiz) {
        dsteqr_(compz, n, &d__[1], &e[1], &z__[z_offset], ldz, &work[1], info, (ftnlen)1);
    } else {
        if (icompz == 1) {
            storez = *n * *n + 1;
        } else {
            storez = 1;
        }
        if (icompz == 2) {
            dlaset_((char *)"Full", n, n, &c_b17, &c_b18, &z__[z_offset], ldz, (ftnlen)4);
        }
        orgnrm = dlanst_((char *)"M", n, &d__[1], &e[1], (ftnlen)1);
        if (orgnrm == 0.) {
            goto L50;
        }
        eps = dlamch_((char *)"Epsilon", (ftnlen)7);
        start = 1;
    L10:
        if (start <= *n) {
            finish = start;
        L20:
            if (finish < *n) {
                tiny = eps * sqrt((d__1 = d__[finish], abs(d__1))) *
                       sqrt((d__2 = d__[finish + 1], abs(d__2)));
                if ((d__1 = e[finish], abs(d__1)) > tiny) {
                    ++finish;
                    goto L20;
                }
            }
            m = finish - start + 1;
            if (m == 1) {
                start = finish + 1;
                goto L10;
            }
            if (m > smlsiz) {
                orgnrm = dlanst_((char *)"M", &m, &d__[start], &e[start], (ftnlen)1);
                dlascl_((char *)"G", &c__0, &c__0, &orgnrm, &c_b18, &m, &c__1, &d__[start], &m, info,
                        (ftnlen)1);
                i__1 = m - 1;
                i__2 = m - 1;
                dlascl_((char *)"G", &c__0, &c__0, &orgnrm, &c_b18, &i__1, &c__1, &e[start], &i__2, info,
                        (ftnlen)1);
                if (icompz == 1) {
                    strtrw = 1;
                } else {
                    strtrw = start;
                }
                dlaed0_(&icompz, n, &m, &d__[start], &e[start], &z__[strtrw + start * z_dim1], ldz,
                        &work[1], n, &work[storez], &iwork[1], info);
                if (*info != 0) {
                    *info = (*info / (m + 1) + start - 1) * (*n + 1) + *info % (m + 1) + start - 1;
                    goto L50;
                }
                dlascl_((char *)"G", &c__0, &c__0, &c_b18, &orgnrm, &m, &c__1, &d__[start], &m, info,
                        (ftnlen)1);
            } else {
                if (icompz == 1) {
                    dsteqr_((char *)"I", &m, &d__[start], &e[start], &work[1], &m, &work[m * m + 1], info,
                            (ftnlen)1);
                    dlacpy_((char *)"A", n, &m, &z__[start * z_dim1 + 1], ldz, &work[storez], n, (ftnlen)1);
                    dgemm_((char *)"N", (char *)"N", n, &m, &m, &c_b18, &work[storez], n, &work[1], &m, &c_b17,
                           &z__[start * z_dim1 + 1], ldz, (ftnlen)1, (ftnlen)1);
                } else if (icompz == 2) {
                    dsteqr_((char *)"I", &m, &d__[start], &e[start], &z__[start + start * z_dim1], ldz,
                            &work[1], info, (ftnlen)1);
                } else {
                    dsterf_(&m, &d__[start], &e[start], info);
                }
                if (*info != 0) {
                    *info = start * (*n + 1) + finish;
                    goto L50;
                }
            }
            start = finish + 1;
            goto L10;
        }
        if (icompz == 0) {
            dlasrt_((char *)"I", n, &d__[1], info, (ftnlen)1);
        } else {
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
                    dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[k * z_dim1 + 1], &c__1);
                }
            }
        }
    }
L50:
    work[1] = (doublereal)lwmin;
    iwork[1] = liwmin;
    return 0;
}
#ifdef __cplusplus
}
#endif
