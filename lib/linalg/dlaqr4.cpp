#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__13 = 13;
static integer c__15 = 15;
static integer c_n1 = -1;
static integer c__12 = 12;
static integer c__14 = 14;
static integer c__16 = 16;
static logical c_false = FALSE_;
static integer c__1 = 1;
static integer c__3 = 3;
int dlaqr4_(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, doublereal *h__,
            integer *ldh, doublereal *wr, doublereal *wi, integer *iloz, integer *ihiz,
            doublereal *z__, integer *ldz, doublereal *work, integer *lwork, integer *info)
{
    integer h_dim1, h_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;
    integer i__, k;
    doublereal aa, bb, cc, dd;
    integer ld;
    doublereal cs;
    integer nh, it, ks, kt;
    doublereal sn;
    integer ku, kv, ls, ns;
    doublereal ss;
    integer nw, inf, kdu, nho, nve, kwh, nsr, nwr, kwv, ndec, ndfl, kbot, nmin;
    doublereal swap;
    integer ktop;
    doublereal zdum[1];
    integer kacc22, itmax, nsmax, nwmax, kwtop;
    extern int dlaqr2_(logical *, logical *, integer *, integer *, integer *, integer *,
                       doublereal *, integer *, integer *, integer *, doublereal *, integer *,
                       integer *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                       integer *, doublereal *, integer *, integer *, doublereal *, integer *,
                       doublereal *, integer *),
        dlanv2_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                doublereal *, doublereal *, doublereal *, doublereal *),
        dlaqr5_(logical *, logical *, integer *, integer *, integer *, integer *, integer *,
                doublereal *, doublereal *, doublereal *, integer *, integer *, integer *,
                doublereal *, integer *, doublereal *, integer *, doublereal *, integer *,
                integer *, doublereal *, integer *, integer *, doublereal *, integer *);
    integer nibble;
    extern int dlahqr_(logical *, logical *, integer *, integer *, integer *, doublereal *,
                       integer *, doublereal *, doublereal *, integer *, integer *, doublereal *,
                       integer *, integer *),
        dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *,
                ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    char jbcmpz[2];
    integer nwupbd;
    logical sorted;
    integer lwkopt;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --wr;
    --wi;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    *info = 0;
    if (*n == 0) {
        work[1] = 1.;
        return 0;
    }
    if (*n <= 15) {
        lwkopt = 1;
        if (*lwork != -1) {
            dlahqr_(wantt, wantz, n, ilo, ihi, &h__[h_offset], ldh, &wr[1], &wi[1], iloz, ihiz,
                    &z__[z_offset], ldz, info);
        }
    } else {
        *info = 0;
        if (*wantt) {
            *(unsigned char *)jbcmpz = 'S';
        } else {
            *(unsigned char *)jbcmpz = 'E';
        }
        if (*wantz) {
            *(unsigned char *)&jbcmpz[1] = 'V';
        } else {
            *(unsigned char *)&jbcmpz[1] = 'N';
        }
        nwr = ilaenv_(&c__13, (char *)"DLAQR4", jbcmpz, n, ilo, ihi, lwork, (ftnlen)6, (ftnlen)2);
        nwr = max(2, nwr);
        i__1 = *ihi - *ilo + 1, i__2 = (*n - 1) / 3, i__1 = min(i__1, i__2);
        nwr = min(i__1, nwr);
        nsr = ilaenv_(&c__15, (char *)"DLAQR4", jbcmpz, n, ilo, ihi, lwork, (ftnlen)6, (ftnlen)2);
        i__1 = nsr, i__2 = (*n - 3) / 6, i__1 = min(i__1, i__2), i__2 = *ihi - *ilo;
        nsr = min(i__1, i__2);
        i__1 = 2, i__2 = nsr - nsr % 2;
        nsr = max(i__1, i__2);
        i__1 = nwr + 1;
        dlaqr2_(wantt, wantz, n, ilo, ihi, &i__1, &h__[h_offset], ldh, iloz, ihiz, &z__[z_offset],
                ldz, &ls, &ld, &wr[1], &wi[1], &h__[h_offset], ldh, n, &h__[h_offset], ldh, n,
                &h__[h_offset], ldh, &work[1], &c_n1);
        i__1 = nsr * 3 / 2, i__2 = (integer)work[1];
        lwkopt = max(i__1, i__2);
        if (*lwork == -1) {
            work[1] = (doublereal)lwkopt;
            return 0;
        }
        nmin = ilaenv_(&c__12, (char *)"DLAQR4", jbcmpz, n, ilo, ihi, lwork, (ftnlen)6, (ftnlen)2);
        nmin = max(15, nmin);
        nibble = ilaenv_(&c__14, (char *)"DLAQR4", jbcmpz, n, ilo, ihi, lwork, (ftnlen)6, (ftnlen)2);
        nibble = max(0, nibble);
        kacc22 = ilaenv_(&c__16, (char *)"DLAQR4", jbcmpz, n, ilo, ihi, lwork, (ftnlen)6, (ftnlen)2);
        kacc22 = max(0, kacc22);
        kacc22 = min(2, kacc22);
        i__1 = (*n - 1) / 3, i__2 = *lwork / 2;
        nwmax = min(i__1, i__2);
        nw = nwmax;
        i__1 = (*n - 3) / 6, i__2 = (*lwork << 1) / 3;
        nsmax = min(i__1, i__2);
        nsmax -= nsmax % 2;
        ndfl = 1;
        i__1 = 10, i__2 = *ihi - *ilo + 1;
        itmax = max(i__1, i__2) * 30;
        kbot = *ihi;
        i__1 = itmax;
        for (it = 1; it <= i__1; ++it) {
            if (kbot < *ilo) {
                goto L90;
            }
            i__2 = *ilo + 1;
            for (k = kbot; k >= i__2; --k) {
                if (h__[k + (k - 1) * h_dim1] == 0.) {
                    goto L20;
                }
            }
            k = *ilo;
        L20:
            ktop = k;
            nh = kbot - ktop + 1;
            nwupbd = min(nh, nwmax);
            if (ndfl < 5) {
                nw = min(nwupbd, nwr);
            } else {
                i__2 = nwupbd, i__3 = nw << 1;
                nw = min(i__2, i__3);
            }
            if (nw < nwmax) {
                if (nw >= nh - 1) {
                    nw = nh;
                } else {
                    kwtop = kbot - nw + 1;
                    if ((d__1 = h__[kwtop + (kwtop - 1) * h_dim1], abs(d__1)) >
                        (d__2 = h__[kwtop - 1 + (kwtop - 2) * h_dim1], abs(d__2))) {
                        ++nw;
                    }
                }
            }
            if (ndfl < 5) {
                ndec = -1;
            } else if (ndec >= 0 || nw >= nwupbd) {
                ++ndec;
                if (nw - ndec < 2) {
                    ndec = 0;
                }
                nw -= ndec;
            }
            kv = *n - nw + 1;
            kt = nw + 1;
            nho = *n - nw - 1 - kt + 1;
            kwv = nw + 2;
            nve = *n - nw - kwv + 1;
            dlaqr2_(wantt, wantz, n, &ktop, &kbot, &nw, &h__[h_offset], ldh, iloz, ihiz,
                    &z__[z_offset], ldz, &ls, &ld, &wr[1], &wi[1], &h__[kv + h_dim1], ldh, &nho,
                    &h__[kv + kt * h_dim1], ldh, &nve, &h__[kwv + h_dim1], ldh, &work[1], lwork);
            kbot -= ld;
            ks = kbot - ls + 1;
            if (ld == 0 || ld * 100 <= nw * nibble && kbot - ktop + 1 > min(nmin, nwmax)) {
                i__4 = 2, i__5 = kbot - ktop;
                i__2 = min(nsmax, nsr), i__3 = max(i__4, i__5);
                ns = min(i__2, i__3);
                ns -= ns % 2;
                if (ndfl % 6 == 0) {
                    ks = kbot - ns + 1;
                    i__3 = ks + 1, i__4 = ktop + 2;
                    i__2 = max(i__3, i__4);
                    for (i__ = kbot; i__ >= i__2; i__ += -2) {
                        ss = (d__1 = h__[i__ + (i__ - 1) * h_dim1], abs(d__1)) +
                             (d__2 = h__[i__ - 1 + (i__ - 2) * h_dim1], abs(d__2));
                        aa = ss * .75 + h__[i__ + i__ * h_dim1];
                        bb = ss;
                        cc = ss * -.4375;
                        dd = aa;
                        dlanv2_(&aa, &bb, &cc, &dd, &wr[i__ - 1], &wi[i__ - 1], &wr[i__], &wi[i__],
                                &cs, &sn);
                    }
                    if (ks == ktop) {
                        wr[ks + 1] = h__[ks + 1 + (ks + 1) * h_dim1];
                        wi[ks + 1] = 0.;
                        wr[ks] = wr[ks + 1];
                        wi[ks] = wi[ks + 1];
                    }
                } else {
                    if (kbot - ks + 1 <= ns / 2) {
                        ks = kbot - ns + 1;
                        kt = *n - ns + 1;
                        dlacpy_((char *)"A", &ns, &ns, &h__[ks + ks * h_dim1], ldh, &h__[kt + h_dim1], ldh,
                                (ftnlen)1);
                        dlahqr_(&c_false, &c_false, &ns, &c__1, &ns, &h__[kt + h_dim1], ldh,
                                &wr[ks], &wi[ks], &c__1, &c__1, zdum, &c__1, &inf);
                        ks += inf;
                        if (ks >= kbot) {
                            aa = h__[kbot - 1 + (kbot - 1) * h_dim1];
                            cc = h__[kbot + (kbot - 1) * h_dim1];
                            bb = h__[kbot - 1 + kbot * h_dim1];
                            dd = h__[kbot + kbot * h_dim1];
                            dlanv2_(&aa, &bb, &cc, &dd, &wr[kbot - 1], &wi[kbot - 1], &wr[kbot],
                                    &wi[kbot], &cs, &sn);
                            ks = kbot - 1;
                        }
                    }
                    if (kbot - ks + 1 > ns) {
                        sorted = FALSE_;
                        i__2 = ks + 1;
                        for (k = kbot; k >= i__2; --k) {
                            if (sorted) {
                                goto L60;
                            }
                            sorted = TRUE_;
                            i__3 = k - 1;
                            for (i__ = ks; i__ <= i__3; ++i__) {
                                if ((d__1 = wr[i__], abs(d__1)) + (d__2 = wi[i__], abs(d__2)) <
                                    (d__3 = wr[i__ + 1], abs(d__3)) +
                                        (d__4 = wi[i__ + 1], abs(d__4))) {
                                    sorted = FALSE_;
                                    swap = wr[i__];
                                    wr[i__] = wr[i__ + 1];
                                    wr[i__ + 1] = swap;
                                    swap = wi[i__];
                                    wi[i__] = wi[i__ + 1];
                                    wi[i__ + 1] = swap;
                                }
                            }
                        }
                    L60:;
                    }
                    i__2 = ks + 2;
                    for (i__ = kbot; i__ >= i__2; i__ += -2) {
                        if (wi[i__] != -wi[i__ - 1]) {
                            swap = wr[i__];
                            wr[i__] = wr[i__ - 1];
                            wr[i__ - 1] = wr[i__ - 2];
                            wr[i__ - 2] = swap;
                            swap = wi[i__];
                            wi[i__] = wi[i__ - 1];
                            wi[i__ - 1] = wi[i__ - 2];
                            wi[i__ - 2] = swap;
                        }
                    }
                }
                if (kbot - ks + 1 == 2) {
                    if (wi[kbot] == 0.) {
                        if ((d__1 = wr[kbot] - h__[kbot + kbot * h_dim1], abs(d__1)) <
                            (d__2 = wr[kbot - 1] - h__[kbot + kbot * h_dim1], abs(d__2))) {
                            wr[kbot - 1] = wr[kbot];
                        } else {
                            wr[kbot] = wr[kbot - 1];
                        }
                    }
                }
                i__2 = ns, i__3 = kbot - ks + 1;
                ns = min(i__2, i__3);
                ns -= ns % 2;
                ks = kbot - ns + 1;
                kdu = ns << 1;
                ku = *n - kdu + 1;
                kwh = kdu + 1;
                nho = *n - kdu - 3 - (kdu + 1) + 1;
                kwv = kdu + 4;
                nve = *n - kdu - kwv + 1;
                dlaqr5_(wantt, wantz, &kacc22, n, &ktop, &kbot, &ns, &wr[ks], &wi[ks],
                        &h__[h_offset], ldh, iloz, ihiz, &z__[z_offset], ldz, &work[1], &c__3,
                        &h__[ku + h_dim1], ldh, &nve, &h__[kwv + h_dim1], ldh, &nho,
                        &h__[ku + kwh * h_dim1], ldh);
            }
            if (ld > 0) {
                ndfl = 1;
            } else {
                ++ndfl;
            }
        }
        *info = kbot;
    L90:;
    }
    work[1] = (doublereal)lwkopt;
    return 0;
}
#ifdef __cplusplus
}
#endif
