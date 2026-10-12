#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static doublereal c_b18 = .001;
int dstemr_(char *jobz, char *range, integer *n, doublereal *d__, doublereal *e, doublereal *vl,
            doublereal *vu, integer *il, integer *iu, integer *m, doublereal *w, doublereal *z__,
            integer *ldz, integer *nzc, integer *isuppz, logical *tryrac, doublereal *work,
            integer *lwork, integer *iwork, integer *liwork, integer *info, ftnlen jobz_len,
            ftnlen range_len)
{
    integer z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;
    double sqrt(doublereal);
    integer i__, j;
    doublereal r1, r2;
    integer jj;
    doublereal cs;
    integer in;
    doublereal sn, wl, wu;
    integer iil, iiu;
    doublereal eps, tmp;
    integer indd, iend, jblk, wend;
    doublereal rmin, rmax;
    integer itmp;
    doublereal tnrm;
    extern int dlae2_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
    integer inde2, itmp2;
    doublereal rtol1, rtol2;
    extern int dscal_(integer *, doublereal *, doublereal *, integer *);
    doublereal scale;
    integer indgp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer iinfo, iindw, ilast;
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *),
        dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer lwmin;
    logical wantz;
    extern int dlaev2_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                       doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    logical alleig;
    integer ibegin;
    logical indeig;
    integer iindbl;
    logical valeig;
    extern int dlarrc_(char *, integer *, doublereal *, doublereal *, doublereal *, doublereal *,
                       doublereal *, integer *, integer *, integer *, integer *, ftnlen),
        dlarre_(char *, integer *, doublereal *, doublereal *, integer *, integer *, doublereal *,
                doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, integer *,
                integer *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                integer *, doublereal *, doublereal *, doublereal *, integer *, integer *, ftnlen);
    integer wbegin;
    doublereal safmin;
    extern int dlarrj_(integer *, doublereal *, doublereal *, integer *, integer *, doublereal *,
                       integer *, doublereal *, doublereal *, doublereal *, integer *, doublereal *,
                       doublereal *, integer *),
        xerbla_(char *, integer *, ftnlen);
    doublereal bignum;
    integer inderr, iindwk, indgrs, offset;
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *, ftnlen);
    extern int dlarrr_(integer *, doublereal *, doublereal *, integer *),
        dlarrv_(integer *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                integer *, integer *, integer *, integer *, doublereal *, doublereal *,
                doublereal *, doublereal *, doublereal *, doublereal *, integer *, integer *,
                doublereal *, doublereal *, integer *, integer *, doublereal *, integer *,
                integer *),
        dlasrt_(char *, integer *, doublereal *, integer *, ftnlen);
    doublereal thresh;
    integer iinspl, ifirst, indwrk, liwmin, nzcmin;
    doublereal pivmin;
    integer nsplit;
    doublereal smlnum;
    logical lquery, zquery, laeswap;
    --d__;
    --e;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --isuppz;
    --work;
    --iwork;
    wantz = lsame_(jobz, (char *)"V", (ftnlen)1, (ftnlen)1);
    alleig = lsame_(range, (char *)"A", (ftnlen)1, (ftnlen)1);
    valeig = lsame_(range, (char *)"V", (ftnlen)1, (ftnlen)1);
    indeig = lsame_(range, (char *)"I", (ftnlen)1, (ftnlen)1);
    lquery = *lwork == -1 || *liwork == -1;
    zquery = *nzc == -1;
    laeswap = FALSE_;
    if (wantz) {
        lwmin = *n * 18;
        liwmin = *n * 10;
    } else {
        lwmin = *n * 12;
        liwmin = *n << 3;
    }
    wl = 0.;
    wu = 0.;
    iil = 0;
    iiu = 0;
    nsplit = 0;
    if (valeig) {
        wl = *vl;
        wu = *vu;
    } else if (indeig) {
        iil = *il;
        iiu = *iu;
    }
    *info = 0;
    if (!(wantz || lsame_(jobz, (char *)"N", (ftnlen)1, (ftnlen)1))) {
        *info = -1;
    } else if (!(alleig || valeig || indeig)) {
        *info = -2;
    } else if (*n < 0) {
        *info = -3;
    } else if (valeig && *n > 0 && wu <= wl) {
        *info = -7;
    } else if (indeig && (iil < 1 || iil > *n)) {
        *info = -8;
    } else if (indeig && (iiu < iil || iiu > *n)) {
        *info = -9;
    } else if (*ldz < 1 || wantz && *ldz < *n) {
        *info = -13;
    } else if (*lwork < lwmin && !lquery) {
        *info = -17;
    } else if (*liwork < liwmin && !lquery) {
        *info = -19;
    }
    safmin = dlamch_((char *)"Safe minimum", (ftnlen)12);
    eps = dlamch_((char *)"Precision", (ftnlen)9);
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = sqrt(smlnum);
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
    rmax = min(d__1, d__2);
    if (*info == 0) {
        work[1] = (doublereal)lwmin;
        iwork[1] = liwmin;
        if (wantz && alleig) {
            nzcmin = *n;
        } else if (wantz && valeig) {
            dlarrc_((char *)"T", n, vl, vu, &d__[1], &e[1], &safmin, &nzcmin, &itmp, &itmp2, info,
                    (ftnlen)1);
        } else if (wantz && indeig) {
            nzcmin = iiu - iil + 1;
        } else {
            nzcmin = 0;
        }
        if (zquery && *info == 0) {
            z__[z_dim1 + 1] = (doublereal)nzcmin;
        } else if (*nzc < nzcmin && !zquery) {
            *info = -14;
        }
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DSTEMR", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery || zquery) {
        return 0;
    }
    *m = 0;
    if (*n == 0) {
        return 0;
    }
    if (*n == 1) {
        if (alleig || indeig) {
            *m = 1;
            w[1] = d__[1];
        } else {
            if (wl < d__[1] && wu >= d__[1]) {
                *m = 1;
                w[1] = d__[1];
            }
        }
        if (wantz && !zquery) {
            z__[z_dim1 + 1] = 1.;
            isuppz[1] = 1;
            isuppz[2] = 1;
        }
        return 0;
    }
    if (*n == 2) {
        if (!wantz) {
            dlae2_(&d__[1], &e[1], &d__[2], &r1, &r2);
        } else if (wantz && !zquery) {
            dlaev2_(&d__[1], &e[1], &d__[2], &r1, &r2, &cs, &sn);
        }
        if (r1 < r2) {
            e[2] = r1;
            r1 = r2;
            r2 = e[2];
            laeswap = TRUE_;
        }
        if (alleig || valeig && r2 > wl && r2 <= wu || indeig && iil == 1) {
            ++(*m);
            w[*m] = r2;
            if (wantz && !zquery) {
                if (laeswap) {
                    z__[*m * z_dim1 + 1] = cs;
                    z__[*m * z_dim1 + 2] = sn;
                } else {
                    z__[*m * z_dim1 + 1] = -sn;
                    z__[*m * z_dim1 + 2] = cs;
                }
                if (sn != 0.) {
                    if (cs != 0.) {
                        isuppz[(*m << 1) - 1] = 1;
                        isuppz[*m * 2] = 2;
                    } else {
                        isuppz[(*m << 1) - 1] = 1;
                        isuppz[*m * 2] = 1;
                    }
                } else {
                    isuppz[(*m << 1) - 1] = 2;
                    isuppz[*m * 2] = 2;
                }
            }
        }
        if (alleig || valeig && r1 > wl && r1 <= wu || indeig && iiu == 2) {
            ++(*m);
            w[*m] = r1;
            if (wantz && !zquery) {
                if (laeswap) {
                    z__[*m * z_dim1 + 1] = -sn;
                    z__[*m * z_dim1 + 2] = cs;
                } else {
                    z__[*m * z_dim1 + 1] = cs;
                    z__[*m * z_dim1 + 2] = sn;
                }
                if (sn != 0.) {
                    if (cs != 0.) {
                        isuppz[(*m << 1) - 1] = 1;
                        isuppz[*m * 2] = 2;
                    } else {
                        isuppz[(*m << 1) - 1] = 1;
                        isuppz[*m * 2] = 1;
                    }
                } else {
                    isuppz[(*m << 1) - 1] = 2;
                    isuppz[*m * 2] = 2;
                }
            }
        }
    } else {
        indgrs = 1;
        inderr = (*n << 1) + 1;
        indgp = *n * 3 + 1;
        indd = (*n << 2) + 1;
        inde2 = *n * 5 + 1;
        indwrk = *n * 6 + 1;
        iinspl = 1;
        iindbl = *n + 1;
        iindw = (*n << 1) + 1;
        iindwk = *n * 3 + 1;
        scale = 1.;
        tnrm = dlanst_((char *)"M", n, &d__[1], &e[1], (ftnlen)1);
        if (tnrm > 0. && tnrm < rmin) {
            scale = rmin / tnrm;
        } else if (tnrm > rmax) {
            scale = rmax / tnrm;
        }
        if (scale != 1.) {
            dscal_(n, &scale, &d__[1], &c__1);
            i__1 = *n - 1;
            dscal_(&i__1, &scale, &e[1], &c__1);
            tnrm *= scale;
            if (valeig) {
                wl *= scale;
                wu *= scale;
            }
        }
        if (*tryrac) {
            dlarrr_(n, &d__[1], &e[1], &iinfo);
        } else {
            iinfo = -1;
        }
        if (iinfo == 0) {
            thresh = eps;
        } else {
            thresh = -eps;
            *tryrac = FALSE_;
        }
        if (*tryrac) {
            dcopy_(n, &d__[1], &c__1, &work[indd], &c__1);
        }
        i__1 = *n - 1;
        for (j = 1; j <= i__1; ++j) {
            d__1 = e[j];
            work[inde2 + j - 1] = d__1 * d__1;
        }
        if (!wantz) {
            rtol1 = eps * 4.;
            rtol2 = eps * 4.;
        } else {
            rtol1 = sqrt(eps);
            d__1 = sqrt(eps) * .005, d__2 = eps * 4.;
            rtol2 = max(d__1, d__2);
        }
        dlarre_(range, n, &wl, &wu, &iil, &iiu, &d__[1], &e[1], &work[inde2], &rtol1, &rtol2,
                &thresh, &nsplit, &iwork[iinspl], m, &w[1], &work[inderr], &work[indgp],
                &iwork[iindbl], &iwork[iindw], &work[indgrs], &pivmin, &work[indwrk],
                &iwork[iindwk], &iinfo, (ftnlen)1);
        if (iinfo != 0) {
            *info = abs(iinfo) + 10;
            return 0;
        }
        if (wantz) {
            dlarrv_(n, &wl, &wu, &d__[1], &e[1], &pivmin, &iwork[iinspl], m, &c__1, m, &c_b18,
                    &rtol1, &rtol2, &w[1], &work[inderr], &work[indgp], &iwork[iindbl],
                    &iwork[iindw], &work[indgrs], &z__[z_offset], ldz, &isuppz[1], &work[indwrk],
                    &iwork[iindwk], &iinfo);
            if (iinfo != 0) {
                *info = abs(iinfo) + 20;
                return 0;
            }
        } else {
            i__1 = *m;
            for (j = 1; j <= i__1; ++j) {
                itmp = iwork[iindbl + j - 1];
                w[j] += e[iwork[iinspl + itmp - 1]];
            }
        }
        if (*tryrac) {
            ibegin = 1;
            wbegin = 1;
            i__1 = iwork[iindbl + *m - 1];
            for (jblk = 1; jblk <= i__1; ++jblk) {
                iend = iwork[iinspl + jblk - 1];
                in = iend - ibegin + 1;
                wend = wbegin - 1;
            L36:
                if (wend < *m) {
                    if (iwork[iindbl + wend] == jblk) {
                        ++wend;
                        goto L36;
                    }
                }
                if (wend < wbegin) {
                    ibegin = iend + 1;
                    goto L39;
                }
                offset = iwork[iindw + wbegin - 1] - 1;
                ifirst = iwork[iindw + wbegin - 1];
                ilast = iwork[iindw + wend - 1];
                rtol2 = eps * 4.;
                dlarrj_(&in, &work[indd + ibegin - 1], &work[inde2 + ibegin - 1], &ifirst, &ilast,
                        &rtol2, &offset, &w[wbegin], &work[inderr + wbegin - 1], &work[indwrk],
                        &iwork[iindwk], &pivmin, &tnrm, &iinfo);
                ibegin = iend + 1;
                wbegin = wend + 1;
            L39:;
            }
        }
        if (scale != 1.) {
            d__1 = 1. / scale;
            dscal_(m, &d__1, &w[1], &c__1);
        }
    }
    if (nsplit > 1 || *n == 2) {
        if (!wantz) {
            dlasrt_((char *)"I", m, &w[1], &iinfo, (ftnlen)1);
            if (iinfo != 0) {
                *info = 3;
                return 0;
            }
        } else {
            i__1 = *m - 1;
            for (j = 1; j <= i__1; ++j) {
                i__ = 0;
                tmp = w[j];
                i__2 = *m;
                for (jj = j + 1; jj <= i__2; ++jj) {
                    if (w[jj] < tmp) {
                        i__ = jj;
                        tmp = w[jj];
                    }
                }
                if (i__ != 0) {
                    w[i__] = w[j];
                    w[j] = tmp;
                    if (wantz) {
                        dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1], &c__1);
                        itmp = isuppz[(i__ << 1) - 1];
                        isuppz[(i__ << 1) - 1] = isuppz[(j << 1) - 1];
                        isuppz[(j << 1) - 1] = itmp;
                        itmp = isuppz[i__ * 2];
                        isuppz[i__ * 2] = isuppz[j * 2];
                        isuppz[j * 2] = itmp;
                    }
                }
            }
        }
    }
    work[1] = (doublereal)lwmin;
    iwork[1] = liwmin;
    return 0;
}
#ifdef __cplusplus
}
#endif
