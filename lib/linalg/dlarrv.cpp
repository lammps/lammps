#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b5 = 0.;
static integer c__1 = 1;
static integer c__2 = 2;
int dlarrv_(integer *n, doublereal *vl, doublereal *vu, doublereal *d__, doublereal *l,
            doublereal *pivmin, integer *isplit, integer *m, integer *dol, integer *dou,
            doublereal *minrgp, doublereal *rtol1, doublereal *rtol2, doublereal *w,
            doublereal *werr, doublereal *wgap, integer *iblock, integer *indexw, doublereal *gers,
            doublereal *z__, integer *ldz, integer *isuppz, doublereal *work, integer *iwork,
            integer *info)
{
    integer z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;
    logical L__1;
    double log(doublereal);
    integer minwsize, i__, j, k, p, q, miniwsize, ii;
    doublereal gl;
    integer im, in;
    doublereal gu, gap, eps, tau, tol, tmp;
    integer zto;
    doublereal ztz;
    integer iend, jblk;
    doublereal lgap;
    integer done;
    doublereal rgap, left;
    integer wend, iter;
    doublereal bstw;
    integer itmp1;
    extern int dscal_(integer *, doublereal *, doublereal *, integer *);
    integer indld;
    doublereal fudge;
    integer idone;
    doublereal sigma;
    integer iinfo, iindr;
    doublereal resid;
    logical eskip;
    doublereal right;
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer nclus, zfrom;
    doublereal rqtol;
    integer iindc1, iindc2;
    extern int dlar1v_(integer *, integer *, integer *, doublereal *, doublereal *, doublereal *,
                       doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                       logical *, integer *, doublereal *, doublereal *, integer *, integer *,
                       doublereal *, doublereal *, doublereal *, doublereal *);
    logical stp2ii;
    doublereal lambda;
    extern doublereal dlamch_(char *, ftnlen);
    integer ibegin, indeig;
    logical needbs;
    integer indlld;
    doublereal sgndef, mingma;
    extern int dlarrb_(integer *, doublereal *, doublereal *, integer *, integer *, doublereal *,
                       doublereal *, integer *, doublereal *, doublereal *, doublereal *,
                       doublereal *, integer *, doublereal *, doublereal *, integer *, integer *);
    integer oldien, oldncl, wbegin;
    doublereal spdiam;
    integer negcnt;
    extern int dlarrf_(integer *, doublereal *, doublereal *, doublereal *, integer *, integer *,
                       doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                       doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                       doublereal *, integer *);
    integer oldcls;
    doublereal savgap;
    integer ndepth;
    doublereal ssigma;
    extern int dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *,
                       integer *, ftnlen);
    logical usedbs;
    integer iindwk, offset;
    doublereal gaptol;
    integer newcls, oldfst, indwrk, windex, oldlst;
    logical usedrq;
    integer newfst, newftt, parity, windmn, windpl, isupmn, newlst, zusedl;
    doublereal bstres;
    integer newsiz, zusedu, zusedw;
    doublereal nrminv, rqcorr;
    logical tryrqc;
    integer isupmx;
    --d__;
    --l;
    --isplit;
    --w;
    --werr;
    --wgap;
    --iblock;
    --indexw;
    --gers;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --isuppz;
    --work;
    --iwork;
    *info = 0;
    if (*n <= 0 || *m <= 0) {
        return 0;
    }
    indld = *n + 1;
    indlld = (*n << 1) + 1;
    indwrk = *n * 3 + 1;
    minwsize = *n * 12;
    i__1 = minwsize;
    for (i__ = 1; i__ <= i__1; ++i__) {
        work[i__] = 0.;
    }
    iindr = 0;
    iindc1 = *n;
    iindc2 = *n << 1;
    iindwk = *n * 3 + 1;
    miniwsize = *n * 7;
    i__1 = miniwsize;
    for (i__ = 1; i__ <= i__1; ++i__) {
        iwork[i__] = 0;
    }
    zusedl = 1;
    if (*dol > 1) {
        zusedl = *dol - 1;
    }
    zusedu = *m;
    if (*dou < *m) {
        zusedu = *dou + 1;
    }
    zusedw = zusedu - zusedl + 1;
    dlaset_((char *)"Full", n, &zusedw, &c_b5, &c_b5, &z__[zusedl * z_dim1 + 1], ldz, (ftnlen)4);
    eps = dlamch_((char *)"Precision", (ftnlen)9);
    rqtol = eps * 2.;
    tryrqc = TRUE_;
    if (*dol == 1 && *dou == *m) {
    } else {
        *rtol1 = eps * 4.;
        *rtol2 = eps * 4.;
    }
    done = 0;
    ibegin = 1;
    wbegin = 1;
    i__1 = iblock[*m];
    for (jblk = 1; jblk <= i__1; ++jblk) {
        iend = isplit[jblk];
        sigma = l[iend];
        wend = wbegin - 1;
    L15:
        if (wend < *m) {
            if (iblock[wend + 1] == jblk) {
                ++wend;
                goto L15;
            }
        }
        if (wend < wbegin) {
            ibegin = iend + 1;
            goto L170;
        } else if (wend < *dol || wbegin > *dou) {
            ibegin = iend + 1;
            wbegin = wend + 1;
            goto L170;
        }
        gl = gers[(ibegin << 1) - 1];
        gu = gers[ibegin * 2];
        i__2 = iend;
        for (i__ = ibegin + 1; i__ <= i__2; ++i__) {
            d__1 = gers[(i__ << 1) - 1];
            gl = min(d__1, gl);
            d__1 = gers[i__ * 2];
            gu = max(d__1, gu);
        }
        spdiam = gu - gl;
        oldien = ibegin - 1;
        in = iend - ibegin + 1;
        im = wend - wbegin + 1;
        if (ibegin == iend) {
            ++done;
            z__[ibegin + wbegin * z_dim1] = 1.;
            isuppz[(wbegin << 1) - 1] = ibegin;
            isuppz[wbegin * 2] = ibegin;
            w[wbegin] += sigma;
            work[wbegin] = w[wbegin];
            ibegin = iend + 1;
            ++wbegin;
            goto L170;
        }
        dcopy_(&im, &w[wbegin], &c__1, &work[wbegin], &c__1);
        i__2 = im;
        for (i__ = 1; i__ <= i__2; ++i__) {
            w[wbegin + i__ - 1] += sigma;
        }
        ndepth = 0;
        parity = 1;
        nclus = 1;
        iwork[iindc1 + 1] = 1;
        iwork[iindc1 + 2] = im;
        idone = 0;
    L40:
        if (idone < im) {
            if (ndepth > *m) {
                *info = -2;
                return 0;
            }
            oldncl = nclus;
            nclus = 0;
            parity = 1 - parity;
            if (parity == 0) {
                oldcls = iindc1;
                newcls = iindc2;
            } else {
                oldcls = iindc2;
                newcls = iindc1;
            }
            i__2 = oldncl;
            for (i__ = 1; i__ <= i__2; ++i__) {
                j = oldcls + (i__ << 1);
                oldfst = iwork[j - 1];
                oldlst = iwork[j];
                if (ndepth > 0) {
                    if (*dol == 1 && *dou == *m) {
                        j = wbegin + oldfst - 1;
                    } else {
                        if (wbegin + oldfst - 1 < *dol) {
                            j = *dol - 1;
                        } else if (wbegin + oldfst - 1 > *dou) {
                            j = *dou;
                        } else {
                            j = wbegin + oldfst - 1;
                        }
                    }
                    dcopy_(&in, &z__[ibegin + j * z_dim1], &c__1, &d__[ibegin], &c__1);
                    i__3 = in - 1;
                    dcopy_(&i__3, &z__[ibegin + (j + 1) * z_dim1], &c__1, &l[ibegin], &c__1);
                    sigma = z__[iend + (j + 1) * z_dim1];
                    dlaset_((char *)"Full", &in, &c__2, &c_b5, &c_b5, &z__[ibegin + j * z_dim1], ldz,
                            (ftnlen)4);
                }
                i__3 = iend - 1;
                for (j = ibegin; j <= i__3; ++j) {
                    tmp = d__[j] * l[j];
                    work[indld - 1 + j] = tmp;
                    work[indlld - 1 + j] = tmp * l[j];
                }
                if (ndepth > 0) {
                    p = indexw[wbegin - 1 + oldfst];
                    q = indexw[wbegin - 1 + oldlst];
                    offset = indexw[wbegin] - 1;
                    dlarrb_(&in, &d__[ibegin], &work[indlld + ibegin - 1], &p, &q, rtol1, rtol2,
                            &offset, &work[wbegin], &wgap[wbegin], &werr[wbegin], &work[indwrk],
                            &iwork[iindwk], pivmin, &spdiam, &in, &iinfo);
                    if (iinfo != 0) {
                        *info = -1;
                        return 0;
                    }
                    if (oldfst > 1) {
                        d__1 = wgap[wbegin + oldfst - 2],
                        d__2 = w[wbegin + oldfst - 1] - werr[wbegin + oldfst - 1] -
                               w[wbegin + oldfst - 2] - werr[wbegin + oldfst - 2];
                        wgap[wbegin + oldfst - 2] = max(d__1, d__2);
                    }
                    if (wbegin + oldlst - 1 < wend) {
                        d__1 = wgap[wbegin + oldlst - 1],
                        d__2 = w[wbegin + oldlst] - werr[wbegin + oldlst] - w[wbegin + oldlst - 1] -
                               werr[wbegin + oldlst - 1];
                        wgap[wbegin + oldlst - 1] = max(d__1, d__2);
                    }
                    i__3 = oldlst;
                    for (j = oldfst; j <= i__3; ++j) {
                        w[wbegin + j - 1] = work[wbegin + j - 1] + sigma;
                    }
                }
                newfst = oldfst;
                i__3 = oldlst;
                for (j = oldfst; j <= i__3; ++j) {
                    if (j == oldlst) {
                        newlst = j;
                    } else if (wgap[wbegin + j - 1] >=
                               *minrgp * (d__1 = work[wbegin + j - 1], abs(d__1))) {
                        newlst = j;
                    } else {
                        goto L140;
                    }
                    newsiz = newlst - newfst + 1;
                    if (*dol == 1 && *dou == *m) {
                        newftt = wbegin + newfst - 1;
                    } else {
                        if (wbegin + newfst - 1 < *dol) {
                            newftt = *dol - 1;
                        } else if (wbegin + newfst - 1 > *dou) {
                            newftt = *dou;
                        } else {
                            newftt = wbegin + newfst - 1;
                        }
                    }
                    if (newsiz > 1) {
                        if (newfst == 1) {
                            d__1 = 0., d__2 = w[wbegin] - werr[wbegin] - *vl;
                            lgap = max(d__1, d__2);
                        } else {
                            lgap = wgap[wbegin + newfst - 2];
                        }
                        rgap = wgap[wbegin + newlst - 1];
                        for (k = 1; k <= 2; ++k) {
                            if (k == 1) {
                                p = indexw[wbegin - 1 + newfst];
                            } else {
                                p = indexw[wbegin - 1 + newlst];
                            }
                            offset = indexw[wbegin] - 1;
                            dlarrb_(&in, &d__[ibegin], &work[indlld + ibegin - 1], &p, &p, &rqtol,
                                    &rqtol, &offset, &work[wbegin], &wgap[wbegin], &werr[wbegin],
                                    &work[indwrk], &iwork[iindwk], pivmin, &spdiam, &in, &iinfo);
                        }
                        if (wbegin + newlst - 1 < *dol || wbegin + newfst - 1 > *dou) {
                            idone = idone + newlst - newfst + 1;
                            goto L139;
                        }
                        dlarrf_(&in, &d__[ibegin], &l[ibegin], &work[indld + ibegin - 1], &newfst,
                                &newlst, &work[wbegin], &wgap[wbegin], &werr[wbegin], &spdiam,
                                &lgap, &rgap, pivmin, &tau, &z__[ibegin + newftt * z_dim1],
                                &z__[ibegin + (newftt + 1) * z_dim1], &work[indwrk], &iinfo);
                        if (iinfo == 0) {
                            ssigma = sigma + tau;
                            z__[iend + (newftt + 1) * z_dim1] = ssigma;
                            i__4 = newlst;
                            for (k = newfst; k <= i__4; ++k) {
                                fudge = eps * 3. * (d__1 = work[wbegin + k - 1], abs(d__1));
                                work[wbegin + k - 1] -= tau;
                                fudge += eps * 4. * (d__1 = work[wbegin + k - 1], abs(d__1));
                                werr[wbegin + k - 1] += fudge;
                            }
                            ++nclus;
                            k = newcls + (nclus << 1);
                            iwork[k - 1] = newfst;
                            iwork[k] = newlst;
                        } else {
                            *info = -2;
                            return 0;
                        }
                    } else {
                        iter = 0;
                        tol = log((doublereal)in) * 4. * eps;
                        k = newfst;
                        windex = wbegin + k - 1;
                        i__4 = windex - 1;
                        windmn = max(i__4, 1);
                        i__4 = windex + 1;
                        windpl = min(i__4, *m);
                        lambda = work[windex];
                        ++done;
                        if (windex < *dol || windex > *dou) {
                            eskip = TRUE_;
                            goto L125;
                        } else {
                            eskip = FALSE_;
                        }
                        left = work[windex] - werr[windex];
                        right = work[windex] + werr[windex];
                        indeig = indexw[windex];
                        if (k == 1) {
                            d__1 = abs(left), d__2 = abs(right);
                            lgap = eps * max(d__1, d__2);
                        } else {
                            lgap = wgap[windmn];
                        }
                        if (k == im) {
                            d__1 = abs(left), d__2 = abs(right);
                            rgap = eps * max(d__1, d__2);
                        } else {
                            rgap = wgap[windex];
                        }
                        gap = min(lgap, rgap);
                        if (k == 1 || k == im) {
                            gaptol = 0.;
                        } else {
                            gaptol = gap * eps;
                        }
                        isupmn = in;
                        isupmx = 1;
                        savgap = wgap[windex];
                        wgap[windex] = gap;
                        usedbs = FALSE_;
                        usedrq = FALSE_;
                        needbs = !tryrqc;
                    L120:
                        if (needbs) {
                            usedbs = TRUE_;
                            itmp1 = iwork[iindr + windex];
                            offset = indexw[wbegin] - 1;
                            d__1 = eps * 2.;
                            dlarrb_(&in, &d__[ibegin], &work[indlld + ibegin - 1], &indeig, &indeig,
                                    &c_b5, &d__1, &offset, &work[wbegin], &wgap[wbegin],
                                    &werr[wbegin], &work[indwrk], &iwork[iindwk], pivmin, &spdiam,
                                    &itmp1, &iinfo);
                            if (iinfo != 0) {
                                *info = -3;
                                return 0;
                            }
                            lambda = work[windex];
                            iwork[iindr + windex] = 0;
                        }
                        L__1 = !usedbs;
                        dlar1v_(&in, &c__1, &in, &lambda, &d__[ibegin], &l[ibegin],
                                &work[indld + ibegin - 1], &work[indlld + ibegin - 1], pivmin,
                                &gaptol, &z__[ibegin + windex * z_dim1], &L__1, &negcnt, &ztz,
                                &mingma, &iwork[iindr + windex], &isuppz[(windex << 1) - 1],
                                &nrminv, &resid, &rqcorr, &work[indwrk]);
                        if (iter == 0) {
                            bstres = resid;
                            bstw = lambda;
                        } else if (resid < bstres) {
                            bstres = resid;
                            bstw = lambda;
                        }
                        i__4 = isupmn, i__5 = isuppz[(windex << 1) - 1];
                        isupmn = min(i__4, i__5);
                        i__4 = isupmx, i__5 = isuppz[windex * 2];
                        isupmx = max(i__4, i__5);
                        ++iter;
                        if (resid > tol * gap && abs(rqcorr) > rqtol * abs(lambda) && !usedbs) {
                            if (indeig <= negcnt) {
                                sgndef = -1.;
                            } else {
                                sgndef = 1.;
                            }
                            if (rqcorr * sgndef >= 0. && lambda + rqcorr <= right &&
                                lambda + rqcorr >= left) {
                                usedrq = TRUE_;
                                if (sgndef == 1.) {
                                    left = lambda;
                                } else {
                                    right = lambda;
                                }
                                work[windex] = (right + left) * .5;
                                lambda += rqcorr;
                                werr[windex] = (right - left) * .5;
                            } else {
                                needbs = TRUE_;
                            }
                            if (right - left < rqtol * abs(lambda)) {
                                usedbs = TRUE_;
                                goto L120;
                            } else if (iter < 10) {
                                goto L120;
                            } else if (iter == 10) {
                                needbs = TRUE_;
                                goto L120;
                            } else {
                                *info = 5;
                                return 0;
                            }
                        } else {
                            stp2ii = FALSE_;
                            if (usedrq && usedbs && bstres <= resid) {
                                lambda = bstw;
                                stp2ii = TRUE_;
                            }
                            if (stp2ii) {
                                L__1 = !usedbs;
                                dlar1v_(&in, &c__1, &in, &lambda, &d__[ibegin], &l[ibegin],
                                        &work[indld + ibegin - 1], &work[indlld + ibegin - 1],
                                        pivmin, &gaptol, &z__[ibegin + windex * z_dim1], &L__1,
                                        &negcnt, &ztz, &mingma, &iwork[iindr + windex],
                                        &isuppz[(windex << 1) - 1], &nrminv, &resid, &rqcorr,
                                        &work[indwrk]);
                            }
                            work[windex] = lambda;
                        }
                        isuppz[(windex << 1) - 1] += oldien;
                        isuppz[windex * 2] += oldien;
                        zfrom = isuppz[(windex << 1) - 1];
                        zto = isuppz[windex * 2];
                        isupmn += oldien;
                        isupmx += oldien;
                        if (isupmn < zfrom) {
                            i__4 = zfrom - 1;
                            for (ii = isupmn; ii <= i__4; ++ii) {
                                z__[ii + windex * z_dim1] = 0.;
                            }
                        }
                        if (isupmx > zto) {
                            i__4 = isupmx;
                            for (ii = zto + 1; ii <= i__4; ++ii) {
                                z__[ii + windex * z_dim1] = 0.;
                            }
                        }
                        i__4 = zto - zfrom + 1;
                        dscal_(&i__4, &nrminv, &z__[zfrom + windex * z_dim1], &c__1);
                    L125:
                        w[windex] = lambda + sigma;
                        if (!eskip) {
                            if (k > 1) {
                                d__1 = wgap[windmn],
                                d__2 = w[windex] - werr[windex] - w[windmn] - werr[windmn];
                                wgap[windmn] = max(d__1, d__2);
                            }
                            if (windex < wend) {
                                d__1 = savgap,
                                d__2 = w[windpl] - werr[windpl] - w[windex] - werr[windex];
                                wgap[windex] = max(d__1, d__2);
                            }
                        }
                        ++idone;
                    }
                L139:
                    newfst = j + 1;
                L140:;
                }
            }
            ++ndepth;
            goto L40;
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
