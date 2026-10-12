#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
static integer c__0 = 0;
int dlarrd_(char *range, char *order, integer *n, doublereal *vl, doublereal *vu, integer *il,
            integer *iu, doublereal *gers, doublereal *reltol, doublereal *d__, doublereal *e,
            doublereal *e2, doublereal *pivmin, integer *nsplit, integer *isplit, integer *m,
            doublereal *w, doublereal *werr, doublereal *wl, doublereal *wu, integer *iblock,
            integer *indexw, doublereal *work, integer *iwork, integer *info, ftnlen range_len,
            ftnlen order_len)
{
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;
    double log(doublereal);
    integer i__, j, ib, ie, je, nb;
    doublereal gl;
    integer im, in;
    doublereal gu;
    integer iw, jee;
    doublereal eps;
    integer nwl;
    doublereal wlu, wul;
    integer nwu;
    doublereal tmp1, tmp2;
    integer iend, jblk, ioff, iout, itmp1, itmp2, jdisc;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer iinfo;
    doublereal atoli;
    integer iwoff, itmax;
    doublereal wkill, rtoli, uflow, tnorm;
    extern doublereal dlamch_(char *, ftnlen);
    integer ibegin;
    extern int dlaebz_(integer *, integer *, integer *, integer *, integer *, integer *,
                       doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                       doublereal *, integer *, doublereal *, doublereal *, integer *, integer *,
                       doublereal *, integer *, integer *);
    integer irange, idiscl, idumma[1];
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    integer idiscu;
    logical ncnvrg, toofew;
    --iwork;
    --work;
    --indexw;
    --iblock;
    --werr;
    --w;
    --isplit;
    --e2;
    --e;
    --d__;
    --gers;
    *info = 0;
    *m = 0;
    if (*n <= 0) {
        return 0;
    }
    if (lsame_(range, (char *)"A", (ftnlen)1, (ftnlen)1)) {
        irange = 1;
    } else if (lsame_(range, (char *)"V", (ftnlen)1, (ftnlen)1)) {
        irange = 2;
    } else if (lsame_(range, (char *)"I", (ftnlen)1, (ftnlen)1)) {
        irange = 3;
    } else {
        irange = 0;
    }
    if (irange <= 0) {
        *info = -1;
    } else if (!(lsame_(order, (char *)"B", (ftnlen)1, (ftnlen)1) ||
                 lsame_(order, (char *)"E", (ftnlen)1, (ftnlen)1))) {
        *info = -2;
    } else if (*n < 0) {
        *info = -3;
    } else if (irange == 2) {
        if (*vl >= *vu) {
            *info = -5;
        }
    } else if (irange == 3 && (*il < 1 || *il > max(1, *n))) {
        *info = -6;
    } else if (irange == 3 && (*iu < min(*n, *il) || *iu > *n)) {
        *info = -7;
    }
    if (*info != 0) {
        return 0;
    }
    ncnvrg = FALSE_;
    toofew = FALSE_;
    if (irange == 3 && *il == 1 && *iu == *n) {
        irange = 1;
    }
    eps = dlamch_((char *)"P", (ftnlen)1);
    uflow = dlamch_((char *)"U", (ftnlen)1);
    if (*n == 1) {
        if (irange == 1 || irange == 2 && d__[1] > *vl && d__[1] <= *vu ||
            irange == 3 && *il == 1 && *iu == 1) {
            *m = 1;
            w[1] = d__[1];
            werr[1] = 0.;
            iblock[1] = 1;
            indexw[1] = 1;
        }
        return 0;
    }
    nb = ilaenv_(&c__1, (char *)"DSTEBZ", (char *)" ", n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
    if (nb <= 1) {
        nb = 0;
    }
    gl = d__[1];
    gu = d__[1];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        d__1 = gl, d__2 = gers[(i__ << 1) - 1];
        gl = min(d__1, d__2);
        d__1 = gu, d__2 = gers[i__ * 2];
        gu = max(d__1, d__2);
    }
    d__1 = abs(gl), d__2 = abs(gu);
    tnorm = max(d__1, d__2);
    gl = gl - tnorm * 2. * eps * *n - *pivmin * 4.;
    gu = gu + tnorm * 2. * eps * *n + *pivmin * 4.;
    rtoli = *reltol;
    atoli = uflow * 4. + *pivmin * 4.;
    if (irange == 3) {
        itmax = (integer)((log(tnorm + *pivmin) - log(*pivmin)) / log(2.)) + 2;
        work[*n + 1] = gl;
        work[*n + 2] = gl;
        work[*n + 3] = gu;
        work[*n + 4] = gu;
        work[*n + 5] = gl;
        work[*n + 6] = gu;
        iwork[1] = -1;
        iwork[2] = -1;
        iwork[3] = *n + 1;
        iwork[4] = *n + 1;
        iwork[5] = *il - 1;
        iwork[6] = *iu;
        dlaebz_(&c__3, &itmax, n, &c__2, &c__2, &nb, &atoli, &rtoli, pivmin, &d__[1], &e[1], &e2[1],
                &iwork[5], &work[*n + 1], &work[*n + 5], &iout, &iwork[1], &w[1], &iblock[1],
                &iinfo);
        if (iinfo != 0) {
            *info = iinfo;
            return 0;
        }
        if (iwork[6] == *iu) {
            *wl = work[*n + 1];
            wlu = work[*n + 3];
            nwl = iwork[1];
            *wu = work[*n + 4];
            wul = work[*n + 2];
            nwu = iwork[4];
        } else {
            *wl = work[*n + 2];
            wlu = work[*n + 4];
            nwl = iwork[2];
            *wu = work[*n + 3];
            wul = work[*n + 1];
            nwu = iwork[3];
        }
        if (nwl < 0 || nwl >= *n || nwu < 1 || nwu > *n) {
            *info = 4;
            return 0;
        }
    } else if (irange == 2) {
        *wl = *vl;
        *wu = *vu;
    } else if (irange == 1) {
        *wl = gl;
        *wu = gu;
    }
    *m = 0;
    iend = 0;
    *info = 0;
    nwl = 0;
    nwu = 0;
    i__1 = *nsplit;
    for (jblk = 1; jblk <= i__1; ++jblk) {
        ioff = iend;
        ibegin = ioff + 1;
        iend = isplit[jblk];
        in = iend - ioff;
        if (in == 1) {
            if (*wl >= d__[ibegin] - *pivmin) {
                ++nwl;
            }
            if (*wu >= d__[ibegin] - *pivmin) {
                ++nwu;
            }
            if (irange == 1 || *wl < d__[ibegin] - *pivmin && *wu >= d__[ibegin] - *pivmin) {
                ++(*m);
                w[*m] = d__[ibegin];
                werr[*m] = 0.;
                iblock[*m] = jblk;
                indexw[*m] = 1;
            }
        } else {
            gu = d__[ibegin];
            gl = d__[ibegin];
            tmp1 = 0.;
            i__2 = iend;
            for (j = ibegin; j <= i__2; ++j) {
                d__1 = gl, d__2 = gers[(j << 1) - 1];
                gl = min(d__1, d__2);
                d__1 = gu, d__2 = gers[j * 2];
                gu = max(d__1, d__2);
            }
            gl = gl - tnorm * 2. * eps * in - *pivmin * 2.;
            gu = gu + tnorm * 2. * eps * in + *pivmin * 2.;
            if (irange > 1) {
                if (gu < *wl) {
                    nwl += in;
                    nwu += in;
                    goto L70;
                }
                gl = max(gl, *wl);
                gu = min(gu, *wu);
                if (gl >= gu) {
                    goto L70;
                }
            }
            work[*n + 1] = gl;
            work[*n + in + 1] = gu;
            dlaebz_(&c__1, &c__0, &in, &in, &c__1, &nb, &atoli, &rtoli, pivmin, &d__[ibegin],
                    &e[ibegin], &e2[ibegin], idumma, &work[*n + 1], &work[*n + (in << 1) + 1], &im,
                    &iwork[1], &w[*m + 1], &iblock[*m + 1], &iinfo);
            if (iinfo != 0) {
                *info = iinfo;
                return 0;
            }
            nwl += iwork[1];
            nwu += iwork[in + 1];
            iwoff = *m - iwork[1];
            itmax = (integer)((log(gu - gl + *pivmin) - log(*pivmin)) / log(2.)) + 2;
            dlaebz_(&c__2, &itmax, &in, &in, &c__1, &nb, &atoli, &rtoli, pivmin, &d__[ibegin],
                    &e[ibegin], &e2[ibegin], idumma, &work[*n + 1], &work[*n + (in << 1) + 1],
                    &iout, &iwork[1], &w[*m + 1], &iblock[*m + 1], &iinfo);
            if (iinfo != 0) {
                *info = iinfo;
                return 0;
            }
            i__2 = iout;
            for (j = 1; j <= i__2; ++j) {
                tmp1 = (work[j + *n] + work[j + in + *n]) * .5;
                tmp2 = (d__1 = work[j + *n] - work[j + in + *n], abs(d__1)) * .5;
                if (j > iout - iinfo) {
                    ncnvrg = TRUE_;
                    ib = -jblk;
                } else {
                    ib = jblk;
                }
                i__3 = iwork[j + in] + iwoff;
                for (je = iwork[j] + 1 + iwoff; je <= i__3; ++je) {
                    w[je] = tmp1;
                    werr[je] = tmp2;
                    indexw[je] = je - iwoff;
                    iblock[je] = ib;
                }
            }
            *m += im;
        }
    L70:;
    }
    if (irange == 3) {
        idiscl = *il - 1 - nwl;
        idiscu = nwu - *iu;
        if (idiscl > 0) {
            im = 0;
            i__1 = *m;
            for (je = 1; je <= i__1; ++je) {
                if (w[je] <= wlu && idiscl > 0) {
                    --idiscl;
                } else {
                    ++im;
                    w[im] = w[je];
                    werr[im] = werr[je];
                    indexw[im] = indexw[je];
                    iblock[im] = iblock[je];
                }
            }
            *m = im;
        }
        if (idiscu > 0) {
            im = *m + 1;
            for (je = *m; je >= 1; --je) {
                if (w[je] >= wul && idiscu > 0) {
                    --idiscu;
                } else {
                    --im;
                    w[im] = w[je];
                    werr[im] = werr[je];
                    indexw[im] = indexw[je];
                    iblock[im] = iblock[je];
                }
            }
            jee = 0;
            i__1 = *m;
            for (je = im; je <= i__1; ++je) {
                ++jee;
                w[jee] = w[je];
                werr[jee] = werr[je];
                indexw[jee] = indexw[je];
                iblock[jee] = iblock[je];
            }
            *m = *m - im + 1;
        }
        if (idiscl > 0 || idiscu > 0) {
            if (idiscl > 0) {
                wkill = *wu;
                i__1 = idiscl;
                for (jdisc = 1; jdisc <= i__1; ++jdisc) {
                    iw = 0;
                    i__2 = *m;
                    for (je = 1; je <= i__2; ++je) {
                        if (iblock[je] != 0 && (w[je] < wkill || iw == 0)) {
                            iw = je;
                            wkill = w[je];
                        }
                    }
                    iblock[iw] = 0;
                }
            }
            if (idiscu > 0) {
                wkill = *wl;
                i__1 = idiscu;
                for (jdisc = 1; jdisc <= i__1; ++jdisc) {
                    iw = 0;
                    i__2 = *m;
                    for (je = 1; je <= i__2; ++je) {
                        if (iblock[je] != 0 && (w[je] >= wkill || iw == 0)) {
                            iw = je;
                            wkill = w[je];
                        }
                    }
                    iblock[iw] = 0;
                }
            }
            im = 0;
            i__1 = *m;
            for (je = 1; je <= i__1; ++je) {
                if (iblock[je] != 0) {
                    ++im;
                    w[im] = w[je];
                    werr[im] = werr[je];
                    indexw[im] = indexw[je];
                    iblock[im] = iblock[je];
                }
            }
            *m = im;
        }
        if (idiscl < 0 || idiscu < 0) {
            toofew = TRUE_;
        }
    }
    if (irange == 1 && *m != *n || irange == 3 && *m != *iu - *il + 1) {
        toofew = TRUE_;
    }
    if (lsame_(order, (char *)"E", (ftnlen)1, (ftnlen)1) && *nsplit > 1) {
        i__1 = *m - 1;
        for (je = 1; je <= i__1; ++je) {
            ie = 0;
            tmp1 = w[je];
            i__2 = *m;
            for (j = je + 1; j <= i__2; ++j) {
                if (w[j] < tmp1) {
                    ie = j;
                    tmp1 = w[j];
                }
            }
            if (ie != 0) {
                tmp2 = werr[ie];
                itmp1 = iblock[ie];
                itmp2 = indexw[ie];
                w[ie] = w[je];
                werr[ie] = werr[je];
                iblock[ie] = iblock[je];
                indexw[ie] = indexw[je];
                w[je] = tmp1;
                werr[je] = tmp2;
                iblock[je] = itmp1;
                indexw[je] = itmp2;
            }
        }
    }
    *info = 0;
    if (ncnvrg) {
        ++(*info);
    }
    if (toofew) {
        *info += 2;
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
