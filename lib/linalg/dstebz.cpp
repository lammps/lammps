#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
static integer c__0 = 0;
int dstebz_(char *range, char *order, integer *n, doublereal *vl, doublereal *vu, integer *il,
            integer *iu, doublereal *abstol, doublereal *d__, doublereal *e, integer *m,
            integer *nsplit, doublereal *w, integer *iblock, integer *isplit, doublereal *work,
            integer *iwork, integer *info, ftnlen range_len, ftnlen order_len)
{
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, d__3, d__4, d__5;
    double sqrt(doublereal), log(doublereal);
    integer j, ib, jb, ie, je, nb;
    doublereal gl;
    integer im, in;
    doublereal gu;
    integer iw;
    doublereal wl, wu;
    integer nwl;
    doublereal ulp, wlu, wul;
    integer nwu;
    doublereal tmp1, tmp2;
    integer iend, ioff, iout, itmp1, jdisc;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer iinfo;
    doublereal atoli;
    integer iwoff;
    doublereal bnorm;
    integer itmax;
    doublereal wkill, rtoli, tnorm;
    extern doublereal dlamch_(char *, ftnlen);
    integer ibegin;
    extern int dlaebz_(integer *, integer *, integer *, integer *, integer *, integer *,
                       doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                       doublereal *, integer *, doublereal *, doublereal *, integer *, integer *,
                       doublereal *, integer *, integer *);
    integer irange, idiscl;
    doublereal safemn;
    integer idumma[1];
    extern int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    integer idiscu, iorder;
    logical ncnvrg;
    doublereal pivmin;
    logical toofew;
    --iwork;
    --work;
    --isplit;
    --iblock;
    --w;
    --e;
    --d__;
    *info = 0;
    if (lsame_(range, (char *)"A", (ftnlen)1, (ftnlen)1)) {
        irange = 1;
    } else if (lsame_(range, (char *)"V", (ftnlen)1, (ftnlen)1)) {
        irange = 2;
    } else if (lsame_(range, (char *)"I", (ftnlen)1, (ftnlen)1)) {
        irange = 3;
    } else {
        irange = 0;
    }
    if (lsame_(order, (char *)"B", (ftnlen)1, (ftnlen)1)) {
        iorder = 2;
    } else if (lsame_(order, (char *)"E", (ftnlen)1, (ftnlen)1)) {
        iorder = 1;
    } else {
        iorder = 0;
    }
    if (irange <= 0) {
        *info = -1;
    } else if (iorder <= 0) {
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
        i__1 = -(*info);
        xerbla_((char *)"DSTEBZ", &i__1, (ftnlen)6);
        return 0;
    }
    *info = 0;
    ncnvrg = FALSE_;
    toofew = FALSE_;
    *m = 0;
    if (*n == 0) {
        return 0;
    }
    if (irange == 3 && *il == 1 && *iu == *n) {
        irange = 1;
    }
    safemn = dlamch_((char *)"S", (ftnlen)1);
    ulp = dlamch_((char *)"P", (ftnlen)1);
    rtoli = ulp * 2.;
    nb = ilaenv_(&c__1, (char *)"DSTEBZ", (char *)" ", n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
    if (nb <= 1) {
        nb = 0;
    }
    if (*n == 1) {
        *nsplit = 1;
        isplit[1] = 1;
        if (irange == 2 && (*vl >= d__[1] || *vu < d__[1])) {
            *m = 0;
        } else {
            w[1] = d__[1];
            iblock[1] = 1;
            *m = 1;
        }
        return 0;
    }
    *nsplit = 1;
    work[*n] = 0.;
    pivmin = 1.;
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
        d__1 = e[j - 1];
        tmp1 = d__1 * d__1;
        d__2 = ulp;
        if ((d__1 = d__[j] * d__[j - 1], abs(d__1)) * (d__2 * d__2) + safemn > tmp1) {
            isplit[*nsplit] = j - 1;
            ++(*nsplit);
            work[j - 1] = 0.;
        } else {
            work[j - 1] = tmp1;
            pivmin = max(pivmin, tmp1);
        }
    }
    isplit[*nsplit] = *n;
    pivmin *= safemn;
    if (irange == 3) {
        gu = d__[1];
        gl = d__[1];
        tmp1 = 0.;
        i__1 = *n - 1;
        for (j = 1; j <= i__1; ++j) {
            tmp2 = sqrt(work[j]);
            d__1 = gu, d__2 = d__[j] + tmp1 + tmp2;
            gu = max(d__1, d__2);
            d__1 = gl, d__2 = d__[j] - tmp1 - tmp2;
            gl = min(d__1, d__2);
            tmp1 = tmp2;
        }
        d__1 = gu, d__2 = d__[*n] + tmp1;
        gu = max(d__1, d__2);
        d__1 = gl, d__2 = d__[*n] - tmp1;
        gl = min(d__1, d__2);
        d__1 = abs(gl), d__2 = abs(gu);
        tnorm = max(d__1, d__2);
        gl = gl - tnorm * 2.1 * ulp * *n - pivmin * 4.2000000000000002;
        gu = gu + tnorm * 2.1 * ulp * *n + pivmin * 2.1;
        itmax = (integer)((log(tnorm + pivmin) - log(pivmin)) / log(2.)) + 2;
        if (*abstol <= 0.) {
            atoli = ulp * tnorm;
        } else {
            atoli = *abstol;
        }
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
        dlaebz_(&c__3, &itmax, n, &c__2, &c__2, &nb, &atoli, &rtoli, &pivmin, &d__[1], &e[1],
                &work[1], &iwork[5], &work[*n + 1], &work[*n + 5], &iout, &iwork[1], &w[1],
                &iblock[1], &iinfo);
        if (iwork[6] == *iu) {
            wl = work[*n + 1];
            wlu = work[*n + 3];
            nwl = iwork[1];
            wu = work[*n + 4];
            wul = work[*n + 2];
            nwu = iwork[4];
        } else {
            wl = work[*n + 2];
            wlu = work[*n + 4];
            nwl = iwork[2];
            wu = work[*n + 3];
            wul = work[*n + 1];
            nwu = iwork[3];
        }
        if (nwl < 0 || nwl >= *n || nwu < 1 || nwu > *n) {
            *info = 4;
            return 0;
        }
    } else {
        d__3 = abs(d__[1]) + abs(e[1]),
        d__4 = (d__1 = d__[*n], abs(d__1)) + (d__2 = e[*n - 1], abs(d__2));
        tnorm = max(d__3, d__4);
        i__1 = *n - 1;
        for (j = 2; j <= i__1; ++j) {
            d__4 = tnorm, d__5 = (d__1 = d__[j], abs(d__1)) + (d__2 = e[j - 1], abs(d__2)) +
                                 (d__3 = e[j], abs(d__3));
            tnorm = max(d__4, d__5);
        }
        if (*abstol <= 0.) {
            atoli = ulp * tnorm;
        } else {
            atoli = *abstol;
        }
        if (irange == 2) {
            wl = *vl;
            wu = *vu;
        } else {
            wl = 0.;
            wu = 0.;
        }
    }
    *m = 0;
    iend = 0;
    *info = 0;
    nwl = 0;
    nwu = 0;
    i__1 = *nsplit;
    for (jb = 1; jb <= i__1; ++jb) {
        ioff = iend;
        ibegin = ioff + 1;
        iend = isplit[jb];
        in = iend - ioff;
        if (in == 1) {
            if (irange == 1 || wl >= d__[ibegin] - pivmin) {
                ++nwl;
            }
            if (irange == 1 || wu >= d__[ibegin] - pivmin) {
                ++nwu;
            }
            if (irange == 1 || wl < d__[ibegin] - pivmin && wu >= d__[ibegin] - pivmin) {
                ++(*m);
                w[*m] = d__[ibegin];
                iblock[*m] = jb;
            }
        } else {
            gu = d__[ibegin];
            gl = d__[ibegin];
            tmp1 = 0.;
            i__2 = iend - 1;
            for (j = ibegin; j <= i__2; ++j) {
                tmp2 = (d__1 = e[j], abs(d__1));
                d__1 = gu, d__2 = d__[j] + tmp1 + tmp2;
                gu = max(d__1, d__2);
                d__1 = gl, d__2 = d__[j] - tmp1 - tmp2;
                gl = min(d__1, d__2);
                tmp1 = tmp2;
            }
            d__1 = gu, d__2 = d__[iend] + tmp1;
            gu = max(d__1, d__2);
            d__1 = gl, d__2 = d__[iend] - tmp1;
            gl = min(d__1, d__2);
            d__1 = abs(gl), d__2 = abs(gu);
            bnorm = max(d__1, d__2);
            gl = gl - bnorm * 2.1 * ulp * in - pivmin * 2.1;
            gu = gu + bnorm * 2.1 * ulp * in + pivmin * 2.1;
            if (*abstol <= 0.) {
                d__1 = abs(gl), d__2 = abs(gu);
                atoli = ulp * max(d__1, d__2);
            } else {
                atoli = *abstol;
            }
            if (irange > 1) {
                if (gu < wl) {
                    nwl += in;
                    nwu += in;
                    goto L70;
                }
                gl = max(gl, wl);
                gu = min(gu, wu);
                if (gl >= gu) {
                    goto L70;
                }
            }
            work[*n + 1] = gl;
            work[*n + in + 1] = gu;
            dlaebz_(&c__1, &c__0, &in, &in, &c__1, &nb, &atoli, &rtoli, &pivmin, &d__[ibegin],
                    &e[ibegin], &work[ibegin], idumma, &work[*n + 1], &work[*n + (in << 1) + 1],
                    &im, &iwork[1], &w[*m + 1], &iblock[*m + 1], &iinfo);
            nwl += iwork[1];
            nwu += iwork[in + 1];
            iwoff = *m - iwork[1];
            itmax = (integer)((log(gu - gl + pivmin) - log(pivmin)) / log(2.)) + 2;
            dlaebz_(&c__2, &itmax, &in, &in, &c__1, &nb, &atoli, &rtoli, &pivmin, &d__[ibegin],
                    &e[ibegin], &work[ibegin], idumma, &work[*n + 1], &work[*n + (in << 1) + 1],
                    &iout, &iwork[1], &w[*m + 1], &iblock[*m + 1], &iinfo);
            i__2 = iout;
            for (j = 1; j <= i__2; ++j) {
                tmp1 = (work[j + *n] + work[j + in + *n]) * .5;
                if (j > iout - iinfo) {
                    ncnvrg = TRUE_;
                    ib = -jb;
                } else {
                    ib = jb;
                }
                i__3 = iwork[j + in] + iwoff;
                for (je = iwork[j] + 1 + iwoff; je <= i__3; ++je) {
                    w[je] = tmp1;
                    iblock[je] = ib;
                }
            }
            *m += im;
        }
    L70:;
    }
    if (irange == 3) {
        im = 0;
        idiscl = *il - 1 - nwl;
        idiscu = nwu - *iu;
        if (idiscl > 0 || idiscu > 0) {
            i__1 = *m;
            for (je = 1; je <= i__1; ++je) {
                if (w[je] <= wlu && idiscl > 0) {
                    --idiscl;
                } else if (w[je] >= wul && idiscu > 0) {
                    --idiscu;
                } else {
                    ++im;
                    w[im] = w[je];
                    iblock[im] = iblock[je];
                }
            }
            *m = im;
        }
        if (idiscl > 0 || idiscu > 0) {
            if (idiscl > 0) {
                wkill = wu;
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
                wkill = wl;
                i__1 = idiscu;
                for (jdisc = 1; jdisc <= i__1; ++jdisc) {
                    iw = 0;
                    i__2 = *m;
                    for (je = 1; je <= i__2; ++je) {
                        if (iblock[je] != 0 && (w[je] > wkill || iw == 0)) {
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
                    iblock[im] = iblock[je];
                }
            }
            *m = im;
        }
        if (idiscl < 0 || idiscu < 0) {
            toofew = TRUE_;
        }
    }
    if (iorder == 1 && *nsplit > 1) {
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
                itmp1 = iblock[ie];
                w[ie] = w[je];
                iblock[ie] = iblock[je];
                w[je] = tmp1;
                iblock[je] = itmp1;
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
