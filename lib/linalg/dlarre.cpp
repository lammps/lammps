#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c__2 = 2;
int dlarre_(char *range, integer *n, doublereal *vl, doublereal *vu, integer *il, integer *iu,
            doublereal *d__, doublereal *e, doublereal *e2, doublereal *rtol1, doublereal *rtol2,
            doublereal *spltol, integer *nsplit, integer *isplit, integer *m, doublereal *w,
            doublereal *werr, doublereal *wgap, integer *iblock, integer *indexw, doublereal *gers,
            doublereal *pivmin, doublereal *work, integer *iwork, integer *info, ftnlen range_len)
{
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;
    double sqrt(doublereal), log(doublereal);
    integer i__, j;
    doublereal s1, s2;
    integer mb;
    doublereal gl;
    integer in, mm;
    doublereal gu;
    integer cnt;
    doublereal eps, tau, tmp, rtl;
    integer cnt1, cnt2;
    doublereal tmp1, eabs;
    integer iend, jblk;
    doublereal eold;
    integer indl;
    doublereal dmax__, emax;
    integer wend, idum, indu;
    doublereal rtol;
    integer iseed[4];
    doublereal avgap, sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer iinfo;
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    logical norep;
    extern int dlasq2_(integer *, doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    integer ibegin;
    logical forceb;
    integer irange;
    doublereal sgndef;
    extern int dlarra_(integer *, doublereal *, doublereal *, doublereal *, doublereal *,
                       doublereal *, integer *, integer *, integer *),
        dlarrb_(integer *, doublereal *, doublereal *, integer *, integer *, doublereal *,
                doublereal *, integer *, doublereal *, doublereal *, doublereal *, doublereal *,
                integer *, doublereal *, doublereal *, integer *, integer *),
        dlarrc_(char *, integer *, doublereal *, doublereal *, doublereal *, doublereal *,
                doublereal *, integer *, integer *, integer *, integer *, ftnlen);
    integer wbegin;
    extern int dlarrd_(char *, char *, integer *, doublereal *, doublereal *, integer *, integer *,
                       doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                       doublereal *, integer *, integer *, integer *, doublereal *, doublereal *,
                       doublereal *, doublereal *, integer *, integer *, doublereal *, integer *,
                       integer *, ftnlen, ftnlen);
    doublereal safmin, spdiam;
    extern int dlarrk_(integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *,
                       doublereal *, doublereal *, doublereal *, doublereal *, integer *);
    logical usedqd;
    doublereal clwdth, isleft;
    extern int dlarnv_(integer *, integer *, integer *, doublereal *);
    doublereal isrght, bsrtol, dpivot;
    --iwork;
    --work;
    --gers;
    --indexw;
    --iblock;
    --wgap;
    --werr;
    --w;
    --isplit;
    --e2;
    --e;
    --d__;
    *info = 0;
    *nsplit = 0;
    *m = 0;
    if (*n <= 0) {
        return 0;
    }
    if (lsame_(range, (char *)"A", (ftnlen)1, (ftnlen)1)) {
        irange = 1;
    } else if (lsame_(range, (char *)"V", (ftnlen)1, (ftnlen)1)) {
        irange = 3;
    } else if (lsame_(range, (char *)"I", (ftnlen)1, (ftnlen)1)) {
        irange = 2;
    }
    safmin = dlamch_((char *)"S", (ftnlen)1);
    eps = dlamch_((char *)"P", (ftnlen)1);
    rtl = sqrt(eps);
    bsrtol = sqrt(eps);
    if (*n == 1) {
        if (irange == 1 || irange == 3 && d__[1] > *vl && d__[1] <= *vu ||
            irange == 2 && *il == 1 && *iu == 1) {
            *m = 1;
            w[1] = d__[1];
            werr[1] = 0.;
            wgap[1] = 0.;
            iblock[1] = 1;
            indexw[1] = 1;
            gers[1] = d__[1];
            gers[2] = d__[1];
        }
        e[1] = 0.;
        return 0;
    }
    gl = d__[1];
    gu = d__[1];
    eold = 0.;
    emax = 0.;
    e[*n] = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        werr[i__] = 0.;
        wgap[i__] = 0.;
        eabs = (d__1 = e[i__], abs(d__1));
        if (eabs >= emax) {
            emax = eabs;
        }
        tmp1 = eabs + eold;
        gers[(i__ << 1) - 1] = d__[i__] - tmp1;
        d__1 = gl, d__2 = gers[(i__ << 1) - 1];
        gl = min(d__1, d__2);
        gers[i__ * 2] = d__[i__] + tmp1;
        d__1 = gu, d__2 = gers[i__ * 2];
        gu = max(d__1, d__2);
        eold = eabs;
    }
    d__3 = emax;
    d__1 = 1., d__2 = d__3 * d__3;
    *pivmin = safmin * max(d__1, d__2);
    spdiam = gu - gl;
    dlarra_(n, &d__[1], &e[1], &e2[1], spltol, &spdiam, nsplit, &isplit[1], &iinfo);
    forceb = FALSE_;
    usedqd = irange == 1 && !forceb;
    if (irange == 1 && !forceb) {
        *vl = gl;
        *vu = gu;
    } else {
        dlarrd_(range, (char *)"B", n, vl, vu, il, iu, &gers[1], &bsrtol, &d__[1], &e[1], &e2[1], pivmin,
                nsplit, &isplit[1], &mm, &w[1], &werr[1], vl, vu, &iblock[1], &indexw[1], &work[1],
                &iwork[1], &iinfo, (ftnlen)1, (ftnlen)1);
        if (iinfo != 0) {
            *info = -1;
            return 0;
        }
        i__1 = *n;
        for (i__ = mm + 1; i__ <= i__1; ++i__) {
            w[i__] = 0.;
            werr[i__] = 0.;
            iblock[i__] = 0;
            indexw[i__] = 0;
        }
    }
    ibegin = 1;
    wbegin = 1;
    i__1 = *nsplit;
    for (jblk = 1; jblk <= i__1; ++jblk) {
        iend = isplit[jblk];
        in = iend - ibegin + 1;
        if (in == 1) {
            if (irange == 1 || irange == 3 && d__[ibegin] > *vl && d__[ibegin] <= *vu ||
                irange == 2 && iblock[wbegin] == jblk) {
                ++(*m);
                w[*m] = d__[ibegin];
                werr[*m] = 0.;
                wgap[*m] = 0.;
                iblock[*m] = jblk;
                indexw[*m] = 1;
                ++wbegin;
            }
            e[iend] = 0.;
            ibegin = iend + 1;
            goto L170;
        }
        e[iend] = 0.;
        gl = d__[ibegin];
        gu = d__[ibegin];
        i__2 = iend;
        for (i__ = ibegin; i__ <= i__2; ++i__) {
            d__1 = gers[(i__ << 1) - 1];
            gl = min(d__1, gl);
            d__1 = gers[i__ * 2];
            gu = max(d__1, gu);
        }
        spdiam = gu - gl;
        if (!(irange == 1 && !forceb)) {
            mb = 0;
            i__2 = mm;
            for (i__ = wbegin; i__ <= i__2; ++i__) {
                if (iblock[i__] == jblk) {
                    ++mb;
                } else {
                    goto L21;
                }
            }
        L21:
            if (mb == 0) {
                e[iend] = 0.;
                ibegin = iend + 1;
                goto L170;
            } else {
                usedqd = (doublereal)mb > in * .5 && !forceb;
                wend = wbegin + mb - 1;
                sigma = 0.;
                i__2 = wend - 1;
                for (i__ = wbegin; i__ <= i__2; ++i__) {
                    d__1 = 0., d__2 = w[i__ + 1] - werr[i__ + 1] - (w[i__] + werr[i__]);
                    wgap[i__] = max(d__1, d__2);
                }
                d__1 = 0., d__2 = *vu - sigma - (w[wend] + werr[wend]);
                wgap[wend] = max(d__1, d__2);
                indl = indexw[wbegin];
                indu = indexw[wend];
            }
        }
        if (irange == 1 && !forceb || usedqd) {
            dlarrk_(&in, &c__1, &gl, &gu, &d__[ibegin], &e2[ibegin], pivmin, &rtl, &tmp, &tmp1,
                    &iinfo);
            if (iinfo != 0) {
                *info = -1;
                return 0;
            }
            d__2 = gl, d__3 = tmp - tmp1 - eps * 100. * (d__1 = tmp - tmp1, abs(d__1));
            isleft = max(d__2, d__3);
            dlarrk_(&in, &in, &gl, &gu, &d__[ibegin], &e2[ibegin], pivmin, &rtl, &tmp, &tmp1,
                    &iinfo);
            if (iinfo != 0) {
                *info = -1;
                return 0;
            }
            d__2 = gu, d__3 = tmp + tmp1 + eps * 100. * (d__1 = tmp + tmp1, abs(d__1));
            isrght = min(d__2, d__3);
            spdiam = isrght - isleft;
        } else {
            d__2 = gl, d__3 = w[wbegin] - werr[wbegin] -
                              eps * 100. * (d__1 = w[wbegin] - werr[wbegin], abs(d__1));
            isleft = max(d__2, d__3);
            d__2 = gu,
            d__3 = w[wend] + werr[wend] + eps * 100. * (d__1 = w[wend] + werr[wend], abs(d__1));
            isrght = min(d__2, d__3);
        }
        if (irange == 1 && !forceb) {
            usedqd = TRUE_;
            indl = 1;
            indu = in;
            mb = in;
            wend = wbegin + mb - 1;
            s1 = isleft + spdiam * .25;
            s2 = isrght - spdiam * .25;
        } else {
            if (usedqd) {
                s1 = isleft + spdiam * .25;
                s2 = isrght - spdiam * .25;
            } else {
                tmp = min(isrght, *vu) - max(isleft, *vl);
                s1 = max(isleft, *vl) + tmp * .25;
                s2 = min(isrght, *vu) - tmp * .25;
            }
        }
        if (mb > 1) {
            dlarrc_((char *)"T", &in, &s1, &s2, &d__[ibegin], &e[ibegin], pivmin, &cnt, &cnt1, &cnt2,
                    &iinfo, (ftnlen)1);
        }
        if (mb == 1) {
            sigma = gl;
            sgndef = 1.;
        } else if (cnt1 - indl >= indu - cnt2) {
            if (irange == 1 && !forceb) {
                sigma = max(isleft, gl);
            } else if (usedqd) {
                sigma = isleft;
            } else {
                sigma = max(isleft, *vl);
            }
            sgndef = 1.;
        } else {
            if (irange == 1 && !forceb) {
                sigma = min(isrght, gu);
            } else if (usedqd) {
                sigma = isrght;
            } else {
                sigma = min(isrght, *vu);
            }
            sgndef = -1.;
        }
        if (usedqd) {
            tau = spdiam * eps * *n + *pivmin * 2.;
            d__1 = tau, d__2 = eps * 2. * abs(sigma);
            tau = max(d__1, d__2);
        } else {
            if (mb > 1) {
                clwdth = w[wend] + werr[wend] - w[wbegin] - werr[wbegin];
                avgap = (d__1 = clwdth / (doublereal)(wend - wbegin), abs(d__1));
                if (sgndef == 1.) {
                    d__1 = wgap[wbegin];
                    tau = max(d__1, avgap) * .5;
                    d__1 = tau, d__2 = werr[wbegin];
                    tau = max(d__1, d__2);
                } else {
                    d__1 = wgap[wend - 1];
                    tau = max(d__1, avgap) * .5;
                    d__1 = tau, d__2 = werr[wend];
                    tau = max(d__1, d__2);
                }
            } else {
                tau = werr[wbegin];
            }
        }
        for (idum = 1; idum <= 6; ++idum) {
            dpivot = d__[ibegin] - sigma;
            work[1] = dpivot;
            dmax__ = abs(work[1]);
            j = ibegin;
            i__2 = in - 1;
            for (i__ = 1; i__ <= i__2; ++i__) {
                work[(in << 1) + i__] = 1. / work[i__];
                tmp = e[j] * work[(in << 1) + i__];
                work[in + i__] = tmp;
                dpivot = d__[j + 1] - sigma - tmp * e[j];
                work[i__ + 1] = dpivot;
                d__1 = dmax__, d__2 = abs(dpivot);
                dmax__ = max(d__1, d__2);
                ++j;
            }
            if (dmax__ > spdiam * 64.) {
                norep = TRUE_;
            } else {
                norep = FALSE_;
            }
            if (usedqd && !norep) {
                i__2 = in;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    tmp = sgndef * work[i__];
                    if (tmp < 0.) {
                        norep = TRUE_;
                    }
                }
            }
            if (norep) {
                if (idum == 5) {
                    if (sgndef == 1.) {
                        sigma = gl - spdiam * 2. * eps * *n - *pivmin * 4.;
                    } else {
                        sigma = gu + spdiam * 2. * eps * *n + *pivmin * 4.;
                    }
                } else {
                    sigma -= sgndef * tau;
                    tau *= 2.;
                }
            } else {
                goto L83;
            }
        }
        *info = 2;
        return 0;
    L83:
        e[iend] = sigma;
        dcopy_(&in, &work[1], &c__1, &d__[ibegin], &c__1);
        i__2 = in - 1;
        dcopy_(&i__2, &work[in + 1], &c__1, &e[ibegin], &c__1);
        if (mb > 1) {
            for (i__ = 1; i__ <= 4; ++i__) {
                iseed[i__ - 1] = 1;
            }
            i__2 = (in << 1) - 1;
            dlarnv_(&c__2, iseed, &i__2, &work[1]);
            i__2 = in - 1;
            for (i__ = 1; i__ <= i__2; ++i__) {
                d__[ibegin + i__ - 1] *= eps * 8. * work[i__] + 1.;
                e[ibegin + i__ - 1] *= eps * 8. * work[in + i__] + 1.;
            }
            d__[iend] *= eps * 4. * work[in] + 1.;
        }
        if (!usedqd) {
            i__2 = wend;
            for (j = wbegin; j <= i__2; ++j) {
                w[j] -= sigma;
                werr[j] += (d__1 = w[j], abs(d__1)) * eps;
            }
            i__2 = iend - 1;
            for (i__ = ibegin; i__ <= i__2; ++i__) {
                d__1 = e[i__];
                work[i__] = d__[i__] * (d__1 * d__1);
            }
            i__2 = indl - 1;
            dlarrb_(&in, &d__[ibegin], &work[ibegin], &indl, &indu, rtol1, rtol2, &i__2, &w[wbegin],
                    &wgap[wbegin], &werr[wbegin], &work[(*n << 1) + 1], &iwork[1], pivmin, &spdiam,
                    &in, &iinfo);
            if (iinfo != 0) {
                *info = -4;
                return 0;
            }
            d__1 = 0., d__2 = *vu - sigma - (w[wend] + werr[wend]);
            wgap[wend] = max(d__1, d__2);
            i__2 = indu;
            for (i__ = indl; i__ <= i__2; ++i__) {
                ++(*m);
                iblock[*m] = jblk;
                indexw[*m] = i__;
            }
        } else {
            rtol = log((doublereal)in) * 4. * eps;
            j = ibegin;
            i__2 = in - 1;
            for (i__ = 1; i__ <= i__2; ++i__) {
                work[(i__ << 1) - 1] = (d__1 = d__[j], abs(d__1));
                work[i__ * 2] = e[j] * e[j] * work[(i__ << 1) - 1];
                ++j;
            }
            work[(in << 1) - 1] = (d__1 = d__[iend], abs(d__1));
            work[in * 2] = 0.;
            dlasq2_(&in, &work[1], &iinfo);
            if (iinfo != 0) {
                *info = -5;
                return 0;
            } else {
                i__2 = in;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    if (work[i__] < 0.) {
                        *info = -6;
                        return 0;
                    }
                }
            }
            if (sgndef > 0.) {
                i__2 = indu;
                for (i__ = indl; i__ <= i__2; ++i__) {
                    ++(*m);
                    w[*m] = work[in - i__ + 1];
                    iblock[*m] = jblk;
                    indexw[*m] = i__;
                }
            } else {
                i__2 = indu;
                for (i__ = indl; i__ <= i__2; ++i__) {
                    ++(*m);
                    w[*m] = -work[i__];
                    iblock[*m] = jblk;
                    indexw[*m] = i__;
                }
            }
            i__2 = *m;
            for (i__ = *m - mb + 1; i__ <= i__2; ++i__) {
                werr[i__] = rtol * (d__1 = w[i__], abs(d__1));
            }
            i__2 = *m - 1;
            for (i__ = *m - mb + 1; i__ <= i__2; ++i__) {
                d__1 = 0., d__2 = w[i__ + 1] - werr[i__ + 1] - (w[i__] + werr[i__]);
                wgap[i__] = max(d__1, d__2);
            }
            d__1 = 0., d__2 = *vu - sigma - (w[*m] + werr[*m]);
            wgap[*m] = max(d__1, d__2);
        }
        ibegin = iend + 1;
        wbegin = wend + 1;
    L170:;
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
