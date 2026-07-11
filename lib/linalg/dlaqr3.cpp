#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static logical c_true = TRUE_;
static doublereal c_b17 = 0.;
static doublereal c_b18 = 1.;
static integer c__12 = 12;
int dlaqr3_(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot, integer *nw,
            doublereal *h__, integer *ldh, integer *iloz, integer *ihiz, doublereal *z__,
            integer *ldz, integer *ns, integer *nd, doublereal *sr, doublereal *si, doublereal *v,
            integer *ldv, integer *nh, doublereal *t, integer *ldt, integer *nv, doublereal *wv,
            integer *ldwv, doublereal *work, integer *lwork)
{
    integer h_dim1, h_offset, t_dim1, t_offset, v_dim1, v_offset, wv_dim1, wv_offset, z_dim1,
        z_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;
    double sqrt(doublereal);
    integer i__, j, k;
    doublereal s, aa, bb, cc, dd, cs, sn;
    integer jw;
    doublereal evi, evk, foo;
    integer kln;
    doublereal tau, ulp;
    integer lwk1, lwk2, lwk3;
    doublereal beta;
    integer kend, kcol, info, nmin, ifst, ilst, ltop, krow;
    extern int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
                      integer *, doublereal *, integer *, doublereal *, doublereal *, integer *,
                      ftnlen, ftnlen);
    logical bulge;
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer infqr, kwtop;
    extern int dlanv2_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                       doublereal *, doublereal *, doublereal *, doublereal *, doublereal *),
        dlaqr4_(logical *, logical *, integer *, integer *, integer *, doublereal *, integer *,
                doublereal *, doublereal *, integer *, integer *, doublereal *, integer *,
                doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern int dgehrd_(integer *, integer *, integer *, doublereal *, integer *, doublereal *,
                       doublereal *, integer *, integer *),
        dlarfg_(integer *, doublereal *, doublereal *, integer *, doublereal *),
        dlahqr_(logical *, logical *, integer *, integer *, integer *, doublereal *, integer *,
                doublereal *, doublereal *, integer *, integer *, doublereal *, integer *,
                integer *),
        dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *,
                ftnlen);
    doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    doublereal safmax;
    extern int dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *,
                       integer *, ftnlen),
        dtrexc_(char *, integer *, doublereal *, integer *, doublereal *, integer *, integer *,
                integer *, doublereal *, integer *, ftnlen),
        dormhr_(char *, char *, integer *, integer *, integer *, integer *, doublereal *, integer *,
                doublereal *, doublereal *, integer *, doublereal *, integer *, integer *, ftnlen,
                ftnlen);
    logical sorted;
    doublereal smlnum;
    extern int dlarf1f_(char *, integer *, integer *, doublereal *, integer *, doublereal *,
                        doublereal *, integer *, doublereal *, ftnlen);
    integer lwkopt;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --sr;
    --si;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    wv_dim1 = *ldwv;
    wv_offset = 1 + wv_dim1;
    wv -= wv_offset;
    --work;
    i__1 = *nw, i__2 = *kbot - *ktop + 1;
    jw = min(i__1, i__2);
    if (jw <= 2) {
        lwkopt = 1;
    } else {
        i__1 = jw - 1;
        dgehrd_(&jw, &c__1, &i__1, &t[t_offset], ldt, &work[1], &work[1], &c_n1, &info);
        lwk1 = (integer)work[1];
        i__1 = jw - 1;
        dormhr_((char *)"R", (char *)"N", &jw, &jw, &c__1, &i__1, &t[t_offset], ldt, &work[1], &v[v_offset], ldv,
                &work[1], &c_n1, &info, (ftnlen)1, (ftnlen)1);
        lwk2 = (integer)work[1];
        dlaqr4_(&c_true, &c_true, &jw, &c__1, &jw, &t[t_offset], ldt, &sr[1], &si[1], &c__1, &jw,
                &v[v_offset], ldv, &work[1], &c_n1, &infqr);
        lwk3 = (integer)work[1];
        i__1 = jw + max(lwk1, lwk2);
        lwkopt = max(i__1, lwk3);
    }
    if (*lwork == -1) {
        work[1] = (doublereal)lwkopt;
        return 0;
    }
    *ns = 0;
    *nd = 0;
    work[1] = 1.;
    if (*ktop > *kbot) {
        return 0;
    }
    if (*nw < 1) {
        return 0;
    }
    safmin = dlamch_((char *)"SAFE MINIMUM", (ftnlen)12);
    safmax = 1. / safmin;
    ulp = dlamch_((char *)"PRECISION", (ftnlen)9);
    smlnum = safmin * ((doublereal)(*n) / ulp);
    i__1 = *nw, i__2 = *kbot - *ktop + 1;
    jw = min(i__1, i__2);
    kwtop = *kbot - jw + 1;
    if (kwtop == *ktop) {
        s = 0.;
    } else {
        s = h__[kwtop + (kwtop - 1) * h_dim1];
    }
    if (*kbot == kwtop) {
        sr[kwtop] = h__[kwtop + kwtop * h_dim1];
        si[kwtop] = 0.;
        *ns = 1;
        *nd = 0;
        d__2 = smlnum, d__3 = ulp * (d__1 = h__[kwtop + kwtop * h_dim1], abs(d__1));
        if (abs(s) <= max(d__2, d__3)) {
            *ns = 0;
            *nd = 1;
            if (kwtop > *ktop) {
                h__[kwtop + (kwtop - 1) * h_dim1] = 0.;
            }
        }
        work[1] = 1.;
        return 0;
    }
    dlacpy_((char *)"U", &jw, &jw, &h__[kwtop + kwtop * h_dim1], ldh, &t[t_offset], ldt, (ftnlen)1);
    i__1 = jw - 1;
    i__2 = *ldh + 1;
    i__3 = *ldt + 1;
    dcopy_(&i__1, &h__[kwtop + 1 + kwtop * h_dim1], &i__2, &t[t_dim1 + 2], &i__3);
    dlaset_((char *)"A", &jw, &jw, &c_b17, &c_b18, &v[v_offset], ldv, (ftnlen)1);
    nmin = ilaenv_(&c__12, (char *)"DLAQR3", (char *)"SV", &jw, &c__1, &jw, lwork, (ftnlen)6, (ftnlen)2);
    if (jw > nmin) {
        dlaqr4_(&c_true, &c_true, &jw, &c__1, &jw, &t[t_offset], ldt, &sr[kwtop], &si[kwtop], &c__1,
                &jw, &v[v_offset], ldv, &work[1], lwork, &infqr);
    } else {
        dlahqr_(&c_true, &c_true, &jw, &c__1, &jw, &t[t_offset], ldt, &sr[kwtop], &si[kwtop], &c__1,
                &jw, &v[v_offset], ldv, &infqr);
    }
    i__1 = jw - 3;
    for (j = 1; j <= i__1; ++j) {
        t[j + 2 + j * t_dim1] = 0.;
        t[j + 3 + j * t_dim1] = 0.;
    }
    if (jw > 2) {
        t[jw + (jw - 2) * t_dim1] = 0.;
    }
    *ns = jw;
    ilst = infqr + 1;
L20:
    if (ilst <= *ns) {
        if (*ns == 1) {
            bulge = FALSE_;
        } else {
            bulge = t[*ns + (*ns - 1) * t_dim1] != 0.;
        }
        if (!bulge) {
            foo = (d__1 = t[*ns + *ns * t_dim1], abs(d__1));
            if (foo == 0.) {
                foo = abs(s);
            }
            d__2 = smlnum, d__3 = ulp * foo;
            if ((d__1 = s * v[*ns * v_dim1 + 1], abs(d__1)) <= max(d__2, d__3)) {
                --(*ns);
            } else {
                ifst = *ns;
                dtrexc_((char *)"V", &jw, &t[t_offset], ldt, &v[v_offset], ldv, &ifst, &ilst, &work[1],
                        &info, (ftnlen)1);
                ++ilst;
            }
        } else {
            foo = (d__3 = t[*ns + *ns * t_dim1], abs(d__3)) +
                  sqrt((d__1 = t[*ns + (*ns - 1) * t_dim1], abs(d__1))) *
                      sqrt((d__2 = t[*ns - 1 + *ns * t_dim1], abs(d__2)));
            if (foo == 0.) {
                foo = abs(s);
            }
            d__3 = (d__1 = s * v[*ns * v_dim1 + 1], abs(d__1)),
            d__4 = (d__2 = s * v[(*ns - 1) * v_dim1 + 1], abs(d__2));
            d__5 = smlnum, d__6 = ulp * foo;
            if (max(d__3, d__4) <= max(d__5, d__6)) {
                *ns += -2;
            } else {
                ifst = *ns;
                dtrexc_((char *)"V", &jw, &t[t_offset], ldt, &v[v_offset], ldv, &ifst, &ilst, &work[1],
                        &info, (ftnlen)1);
                ilst += 2;
            }
        }
        goto L20;
    }
    if (*ns == 0) {
        s = 0.;
    }
    if (*ns < jw) {
        sorted = FALSE_;
        i__ = *ns + 1;
    L30:
        if (sorted) {
            goto L50;
        }
        sorted = TRUE_;
        kend = i__ - 1;
        i__ = infqr + 1;
        if (i__ == *ns) {
            k = i__ + 1;
        } else if (t[i__ + 1 + i__ * t_dim1] == 0.) {
            k = i__ + 1;
        } else {
            k = i__ + 2;
        }
    L40:
        if (k <= kend) {
            if (k == i__ + 1) {
                evi = (d__1 = t[i__ + i__ * t_dim1], abs(d__1));
            } else {
                evi = (d__3 = t[i__ + i__ * t_dim1], abs(d__3)) +
                      sqrt((d__1 = t[i__ + 1 + i__ * t_dim1], abs(d__1))) *
                          sqrt((d__2 = t[i__ + (i__ + 1) * t_dim1], abs(d__2)));
            }
            if (k == kend) {
                evk = (d__1 = t[k + k * t_dim1], abs(d__1));
            } else if (t[k + 1 + k * t_dim1] == 0.) {
                evk = (d__1 = t[k + k * t_dim1], abs(d__1));
            } else {
                evk = (d__3 = t[k + k * t_dim1], abs(d__3)) +
                      sqrt((d__1 = t[k + 1 + k * t_dim1], abs(d__1))) *
                          sqrt((d__2 = t[k + (k + 1) * t_dim1], abs(d__2)));
            }
            if (evi >= evk) {
                i__ = k;
            } else {
                sorted = FALSE_;
                ifst = i__;
                ilst = k;
                dtrexc_((char *)"V", &jw, &t[t_offset], ldt, &v[v_offset], ldv, &ifst, &ilst, &work[1],
                        &info, (ftnlen)1);
                if (info == 0) {
                    i__ = ilst;
                } else {
                    i__ = k;
                }
            }
            if (i__ == kend) {
                k = i__ + 1;
            } else if (t[i__ + 1 + i__ * t_dim1] == 0.) {
                k = i__ + 1;
            } else {
                k = i__ + 2;
            }
            goto L40;
        }
        goto L30;
    L50:;
    }
    i__ = jw;
L60:
    if (i__ >= infqr + 1) {
        if (i__ == infqr + 1) {
            sr[kwtop + i__ - 1] = t[i__ + i__ * t_dim1];
            si[kwtop + i__ - 1] = 0.;
            --i__;
        } else if (t[i__ + (i__ - 1) * t_dim1] == 0.) {
            sr[kwtop + i__ - 1] = t[i__ + i__ * t_dim1];
            si[kwtop + i__ - 1] = 0.;
            --i__;
        } else {
            aa = t[i__ - 1 + (i__ - 1) * t_dim1];
            cc = t[i__ + (i__ - 1) * t_dim1];
            bb = t[i__ - 1 + i__ * t_dim1];
            dd = t[i__ + i__ * t_dim1];
            dlanv2_(&aa, &bb, &cc, &dd, &sr[kwtop + i__ - 2], &si[kwtop + i__ - 2],
                    &sr[kwtop + i__ - 1], &si[kwtop + i__ - 1], &cs, &sn);
            i__ += -2;
        }
        goto L60;
    }
    if (*ns < jw || s == 0.) {
        if (*ns > 1 && s != 0.) {
            dcopy_(ns, &v[v_offset], ldv, &work[1], &c__1);
            beta = work[1];
            dlarfg_(ns, &beta, &work[2], &c__1, &tau);
            i__1 = jw - 2;
            i__2 = jw - 2;
            dlaset_((char *)"L", &i__1, &i__2, &c_b17, &c_b17, &t[t_dim1 + 3], ldt, (ftnlen)1);
            dlarf1f_((char *)"L", ns, &jw, &work[1], &c__1, &tau, &t[t_offset], ldt, &work[jw + 1],
                     (ftnlen)1);
            dlarf1f_((char *)"R", ns, ns, &work[1], &c__1, &tau, &t[t_offset], ldt, &work[jw + 1],
                     (ftnlen)1);
            dlarf1f_((char *)"R", &jw, ns, &work[1], &c__1, &tau, &v[v_offset], ldv, &work[jw + 1],
                     (ftnlen)1);
            i__1 = *lwork - jw;
            dgehrd_(&jw, &c__1, ns, &t[t_offset], ldt, &work[1], &work[jw + 1], &i__1, &info);
        }
        if (kwtop > 1) {
            h__[kwtop + (kwtop - 1) * h_dim1] = s * v[v_dim1 + 1];
        }
        dlacpy_((char *)"U", &jw, &jw, &t[t_offset], ldt, &h__[kwtop + kwtop * h_dim1], ldh, (ftnlen)1);
        i__1 = jw - 1;
        i__2 = *ldt + 1;
        i__3 = *ldh + 1;
        dcopy_(&i__1, &t[t_dim1 + 2], &i__2, &h__[kwtop + 1 + kwtop * h_dim1], &i__3);
        if (*ns > 1 && s != 0.) {
            i__1 = *lwork - jw;
            dormhr_((char *)"R", (char *)"N", &jw, ns, &c__1, ns, &t[t_offset], ldt, &work[1], &v[v_offset], ldv,
                    &work[jw + 1], &i__1, &info, (ftnlen)1, (ftnlen)1);
        }
        if (*wantt) {
            ltop = 1;
        } else {
            ltop = *ktop;
        }
        i__1 = kwtop - 1;
        i__2 = *nv;
        for (krow = ltop; i__2 < 0 ? krow >= i__1 : krow <= i__1; krow += i__2) {
            i__3 = *nv, i__4 = kwtop - krow;
            kln = min(i__3, i__4);
            dgemm_((char *)"N", (char *)"N", &kln, &jw, &jw, &c_b18, &h__[krow + kwtop * h_dim1], ldh, &v[v_offset],
                   ldv, &c_b17, &wv[wv_offset], ldwv, (ftnlen)1, (ftnlen)1);
            dlacpy_((char *)"A", &kln, &jw, &wv[wv_offset], ldwv, &h__[krow + kwtop * h_dim1], ldh,
                    (ftnlen)1);
        }
        if (*wantt) {
            i__2 = *n;
            i__1 = *nh;
            for (kcol = *kbot + 1; i__1 < 0 ? kcol >= i__2 : kcol <= i__2; kcol += i__1) {
                i__3 = *nh, i__4 = *n - kcol + 1;
                kln = min(i__3, i__4);
                dgemm_((char *)"C", (char *)"N", &jw, &kln, &jw, &c_b18, &v[v_offset], ldv,
                       &h__[kwtop + kcol * h_dim1], ldh, &c_b17, &t[t_offset], ldt, (ftnlen)1,
                       (ftnlen)1);
                dlacpy_((char *)"A", &jw, &kln, &t[t_offset], ldt, &h__[kwtop + kcol * h_dim1], ldh,
                        (ftnlen)1);
            }
        }
        if (*wantz) {
            i__1 = *ihiz;
            i__2 = *nv;
            for (krow = *iloz; i__2 < 0 ? krow >= i__1 : krow <= i__1; krow += i__2) {
                i__3 = *nv, i__4 = *ihiz - krow + 1;
                kln = min(i__3, i__4);
                dgemm_((char *)"N", (char *)"N", &kln, &jw, &jw, &c_b18, &z__[krow + kwtop * z_dim1], ldz,
                       &v[v_offset], ldv, &c_b17, &wv[wv_offset], ldwv, (ftnlen)1, (ftnlen)1);
                dlacpy_((char *)"A", &kln, &jw, &wv[wv_offset], ldwv, &z__[krow + kwtop * z_dim1], ldz,
                        (ftnlen)1);
            }
        }
    }
    *nd = jw - *ns;
    *ns -= infqr;
    work[1] = (doublereal)lwkopt;
    return 0;
}
#ifdef __cplusplus
}
#endif
