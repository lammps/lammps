#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b7 = 0.;
static doublereal c_b8 = 1.;
static integer c__2 = 2;
static integer c__1 = 1;
static integer c__3 = 3;
int dlaqr5_(logical *wantt, logical *wantz, integer *kacc22, integer *n, integer *ktop,
            integer *kbot, integer *nshfts, doublereal *sr, doublereal *si, doublereal *h__,
            integer *ldh, integer *iloz, integer *ihiz, doublereal *z__, integer *ldz,
            doublereal *v, integer *ldv, doublereal *u, integer *ldu, integer *nv, doublereal *wv,
            integer *ldwv, integer *nh, doublereal *wh, integer *ldwh)
{
    integer h_dim1, h_offset, u_dim1, u_offset, v_dim1, v_offset, wh_dim1, wh_offset, wv_dim1,
        wv_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    doublereal d__1, d__2, d__3, d__4, d__5;
    integer i__, j, k, m, i2, k1, i4;
    doublereal t1, t2, t3, h11, h12, h21, h22;
    integer m22, ns, nu;
    doublereal vt[3], scl;
    integer kdu, kms;
    doublereal ulp, tst1, tst2, beta;
    logical bmp22;
    integer jcol, jlen, jbot, mbot;
    doublereal swap;
    integer jtop, jrow, mtop;
    doublereal alpha;
    logical accum;
    extern int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
                      integer *, doublereal *, integer *, doublereal *, doublereal *, integer *,
                      ftnlen, ftnlen);
    integer ndcol, incol, krcol, nbmps;
    extern int dlaqr1_(integer *, doublereal *, integer *, doublereal *, doublereal *, doublereal *,
                       doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern int dlarfg_(integer *, doublereal *, doublereal *, integer *, doublereal *),
        dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *,
                ftnlen);
    doublereal safmin;
    extern int dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *,
                       integer *, ftnlen);
    doublereal safmax, refsum, smlnum;
    --sr;
    --si;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    wv_dim1 = *ldwv;
    wv_offset = 1 + wv_dim1;
    wv -= wv_offset;
    wh_dim1 = *ldwh;
    wh_offset = 1 + wh_dim1;
    wh -= wh_offset;
    if (*nshfts < 2) {
        return 0;
    }
    if (*ktop >= *kbot) {
        return 0;
    }
    i__1 = *nshfts - 2;
    for (i__ = 1; i__ <= i__1; i__ += 2) {
        if (si[i__] != -si[i__ + 1]) {
            swap = sr[i__];
            sr[i__] = sr[i__ + 1];
            sr[i__ + 1] = sr[i__ + 2];
            sr[i__ + 2] = swap;
            swap = si[i__];
            si[i__] = si[i__ + 1];
            si[i__ + 1] = si[i__ + 2];
            si[i__ + 2] = swap;
        }
    }
    ns = *nshfts - *nshfts % 2;
    safmin = dlamch_((char *)"SAFE MINIMUM", (ftnlen)12);
    safmax = 1. / safmin;
    ulp = dlamch_((char *)"PRECISION", (ftnlen)9);
    smlnum = safmin * ((doublereal)(*n) / ulp);
    accum = *kacc22 == 1 || *kacc22 == 2;
    if (*ktop + 2 <= *kbot) {
        h__[*ktop + 2 + *ktop * h_dim1] = 0.;
    }
    nbmps = ns / 2;
    kdu = nbmps << 2;
    i__1 = *kbot - 2;
    i__2 = nbmps << 1;
    for (incol = *ktop - (nbmps << 1) + 1; i__2 < 0 ? incol >= i__1 : incol <= i__1;
         incol += i__2) {
        if (accum) {
            jtop = max(*ktop, incol);
        } else if (*wantt) {
            jtop = 1;
        } else {
            jtop = *ktop;
        }
        ndcol = incol + kdu;
        if (accum) {
            dlaset_((char *)"A", &kdu, &kdu, &c_b7, &c_b8, &u[u_offset], ldu, (ftnlen)1);
        }
        i__4 = incol + (nbmps << 1) - 1, i__5 = *kbot - 2;
        i__3 = min(i__4, i__5);
        for (krcol = incol; krcol <= i__3; ++krcol) {
            i__4 = 1, i__5 = (*ktop - krcol) / 2 + 1;
            mtop = max(i__4, i__5);
            i__4 = nbmps, i__5 = (*kbot - krcol - 1) / 2;
            mbot = min(i__4, i__5);
            m22 = mbot + 1;
            bmp22 = mbot < nbmps && krcol + (m22 - 1 << 1) == *kbot - 2;
            if (bmp22) {
                k = krcol + (m22 - 1 << 1);
                if (k == *ktop - 1) {
                    dlaqr1_(&c__2, &h__[k + 1 + (k + 1) * h_dim1], ldh, &sr[(m22 << 1) - 1],
                            &si[(m22 << 1) - 1], &sr[m22 * 2], &si[m22 * 2], &v[m22 * v_dim1 + 1]);
                    beta = v[m22 * v_dim1 + 1];
                    dlarfg_(&c__2, &beta, &v[m22 * v_dim1 + 2], &c__1, &v[m22 * v_dim1 + 1]);
                } else {
                    beta = h__[k + 1 + k * h_dim1];
                    v[m22 * v_dim1 + 2] = h__[k + 2 + k * h_dim1];
                    dlarfg_(&c__2, &beta, &v[m22 * v_dim1 + 2], &c__1, &v[m22 * v_dim1 + 1]);
                    h__[k + 1 + k * h_dim1] = beta;
                    h__[k + 2 + k * h_dim1] = 0.;
                }
                t1 = v[m22 * v_dim1 + 1];
                t2 = t1 * v[m22 * v_dim1 + 2];
                i__5 = *kbot, i__6 = k + 3;
                i__4 = min(i__5, i__6);
                for (j = jtop; j <= i__4; ++j) {
                    refsum =
                        h__[j + (k + 1) * h_dim1] + v[m22 * v_dim1 + 2] * h__[j + (k + 2) * h_dim1];
                    h__[j + (k + 1) * h_dim1] -= refsum * t1;
                    h__[j + (k + 2) * h_dim1] -= refsum * t2;
                }
                if (accum) {
                    jbot = min(ndcol, *kbot);
                } else if (*wantt) {
                    jbot = *n;
                } else {
                    jbot = *kbot;
                }
                t1 = v[m22 * v_dim1 + 1];
                t2 = t1 * v[m22 * v_dim1 + 2];
                i__4 = jbot;
                for (j = k + 1; j <= i__4; ++j) {
                    refsum =
                        h__[k + 1 + j * h_dim1] + v[m22 * v_dim1 + 2] * h__[k + 2 + j * h_dim1];
                    h__[k + 1 + j * h_dim1] -= refsum * t1;
                    h__[k + 2 + j * h_dim1] -= refsum * t2;
                }
                if (k >= *ktop) {
                    if (h__[k + 1 + k * h_dim1] != 0.) {
                        tst1 = (d__1 = h__[k + k * h_dim1], abs(d__1)) +
                               (d__2 = h__[k + 1 + (k + 1) * h_dim1], abs(d__2));
                        if (tst1 == 0.) {
                            if (k >= *ktop + 1) {
                                tst1 += (d__1 = h__[k + (k - 1) * h_dim1], abs(d__1));
                            }
                            if (k >= *ktop + 2) {
                                tst1 += (d__1 = h__[k + (k - 2) * h_dim1], abs(d__1));
                            }
                            if (k >= *ktop + 3) {
                                tst1 += (d__1 = h__[k + (k - 3) * h_dim1], abs(d__1));
                            }
                            if (k <= *kbot - 2) {
                                tst1 += (d__1 = h__[k + 2 + (k + 1) * h_dim1], abs(d__1));
                            }
                            if (k <= *kbot - 3) {
                                tst1 += (d__1 = h__[k + 3 + (k + 1) * h_dim1], abs(d__1));
                            }
                            if (k <= *kbot - 4) {
                                tst1 += (d__1 = h__[k + 4 + (k + 1) * h_dim1], abs(d__1));
                            }
                        }
                        d__2 = smlnum, d__3 = ulp * tst1;
                        if ((d__1 = h__[k + 1 + k * h_dim1], abs(d__1)) <= max(d__2, d__3)) {
                            d__3 = (d__1 = h__[k + 1 + k * h_dim1], abs(d__1)),
                            d__4 = (d__2 = h__[k + (k + 1) * h_dim1], abs(d__2));
                            h12 = max(d__3, d__4);
                            d__3 = (d__1 = h__[k + 1 + k * h_dim1], abs(d__1)),
                            d__4 = (d__2 = h__[k + (k + 1) * h_dim1], abs(d__2));
                            h21 = min(d__3, d__4);
                            d__3 = (d__1 = h__[k + 1 + (k + 1) * h_dim1], abs(d__1)),
                            d__4 = (d__2 = h__[k + k * h_dim1] - h__[k + 1 + (k + 1) * h_dim1],
                                    abs(d__2));
                            h11 = max(d__3, d__4);
                            d__3 = (d__1 = h__[k + 1 + (k + 1) * h_dim1], abs(d__1)),
                            d__4 = (d__2 = h__[k + k * h_dim1] - h__[k + 1 + (k + 1) * h_dim1],
                                    abs(d__2));
                            h22 = min(d__3, d__4);
                            scl = h11 + h12;
                            tst2 = h22 * (h11 / scl);
                            d__1 = smlnum, d__2 = ulp * tst2;
                            if (tst2 == 0. || h21 * (h12 / scl) <= max(d__1, d__2)) {
                                h__[k + 1 + k * h_dim1] = 0.;
                            }
                        }
                    }
                }
                if (accum) {
                    kms = k - incol;
                    t1 = v[m22 * v_dim1 + 1];
                    t2 = t1 * v[m22 * v_dim1 + 2];
                    i__4 = 1, i__5 = *ktop - incol;
                    i__6 = kdu;
                    for (j = max(i__4, i__5); j <= i__6; ++j) {
                        refsum = u[j + (kms + 1) * u_dim1] +
                                 v[m22 * v_dim1 + 2] * u[j + (kms + 2) * u_dim1];
                        u[j + (kms + 1) * u_dim1] -= refsum * t1;
                        u[j + (kms + 2) * u_dim1] -= refsum * t2;
                    }
                } else if (*wantz) {
                    t1 = v[m22 * v_dim1 + 1];
                    t2 = t1 * v[m22 * v_dim1 + 2];
                    i__6 = *ihiz;
                    for (j = *iloz; j <= i__6; ++j) {
                        refsum = z__[j + (k + 1) * z_dim1] +
                                 v[m22 * v_dim1 + 2] * z__[j + (k + 2) * z_dim1];
                        z__[j + (k + 1) * z_dim1] -= refsum * t1;
                        z__[j + (k + 2) * z_dim1] -= refsum * t2;
                    }
                }
            }
            i__6 = mtop;
            for (m = mbot; m >= i__6; --m) {
                k = krcol + (m - 1 << 1);
                if (k == *ktop - 1) {
                    dlaqr1_(&c__3, &h__[*ktop + *ktop * h_dim1], ldh, &sr[(m << 1) - 1],
                            &si[(m << 1) - 1], &sr[m * 2], &si[m * 2], &v[m * v_dim1 + 1]);
                    alpha = v[m * v_dim1 + 1];
                    dlarfg_(&c__3, &alpha, &v[m * v_dim1 + 2], &c__1, &v[m * v_dim1 + 1]);
                } else {
                    t1 = v[m * v_dim1 + 1];
                    t2 = t1 * v[m * v_dim1 + 2];
                    t3 = t1 * v[m * v_dim1 + 3];
                    refsum = v[m * v_dim1 + 3] * h__[k + 3 + (k + 2) * h_dim1];
                    h__[k + 3 + k * h_dim1] = -refsum * t1;
                    h__[k + 3 + (k + 1) * h_dim1] = -refsum * t2;
                    h__[k + 3 + (k + 2) * h_dim1] -= refsum * t3;
                    beta = h__[k + 1 + k * h_dim1];
                    v[m * v_dim1 + 2] = h__[k + 2 + k * h_dim1];
                    v[m * v_dim1 + 3] = h__[k + 3 + k * h_dim1];
                    dlarfg_(&c__3, &beta, &v[m * v_dim1 + 2], &c__1, &v[m * v_dim1 + 1]);
                    if (h__[k + 3 + k * h_dim1] != 0. || h__[k + 3 + (k + 1) * h_dim1] != 0. ||
                        h__[k + 3 + (k + 2) * h_dim1] == 0.) {
                        h__[k + 1 + k * h_dim1] = beta;
                        h__[k + 2 + k * h_dim1] = 0.;
                        h__[k + 3 + k * h_dim1] = 0.;
                    } else {
                        dlaqr1_(&c__3, &h__[k + 1 + (k + 1) * h_dim1], ldh, &sr[(m << 1) - 1],
                                &si[(m << 1) - 1], &sr[m * 2], &si[m * 2], vt);
                        alpha = vt[0];
                        dlarfg_(&c__3, &alpha, &vt[1], &c__1, vt);
                        t1 = vt[0];
                        t2 = t1 * vt[1];
                        t3 = t1 * vt[2];
                        refsum = h__[k + 1 + k * h_dim1] + vt[1] * h__[k + 2 + k * h_dim1];
                        if ((d__1 = h__[k + 2 + k * h_dim1] - refsum * t2, abs(d__1)) +
                                (d__2 = refsum * t3, abs(d__2)) >
                            ulp * ((d__3 = h__[k + k * h_dim1], abs(d__3)) +
                                   (d__4 = h__[k + 1 + (k + 1) * h_dim1], abs(d__4)) +
                                   (d__5 = h__[k + 2 + (k + 2) * h_dim1], abs(d__5)))) {
                            h__[k + 1 + k * h_dim1] = beta;
                            h__[k + 2 + k * h_dim1] = 0.;
                            h__[k + 3 + k * h_dim1] = 0.;
                        } else {
                            h__[k + 1 + k * h_dim1] -= refsum * t1;
                            h__[k + 2 + k * h_dim1] = 0.;
                            h__[k + 3 + k * h_dim1] = 0.;
                            v[m * v_dim1 + 1] = vt[0];
                            v[m * v_dim1 + 2] = vt[1];
                            v[m * v_dim1 + 3] = vt[2];
                        }
                    }
                }
                t1 = v[m * v_dim1 + 1];
                t2 = t1 * v[m * v_dim1 + 2];
                t3 = t1 * v[m * v_dim1 + 3];
                i__5 = *kbot, i__7 = k + 3;
                i__4 = min(i__5, i__7);
                for (j = jtop; j <= i__4; ++j) {
                    refsum = h__[j + (k + 1) * h_dim1] +
                             v[m * v_dim1 + 2] * h__[j + (k + 2) * h_dim1] +
                             v[m * v_dim1 + 3] * h__[j + (k + 3) * h_dim1];
                    h__[j + (k + 1) * h_dim1] -= refsum * t1;
                    h__[j + (k + 2) * h_dim1] -= refsum * t2;
                    h__[j + (k + 3) * h_dim1] -= refsum * t3;
                }
                refsum = h__[k + 1 + (k + 1) * h_dim1] +
                         v[m * v_dim1 + 2] * h__[k + 2 + (k + 1) * h_dim1] +
                         v[m * v_dim1 + 3] * h__[k + 3 + (k + 1) * h_dim1];
                h__[k + 1 + (k + 1) * h_dim1] -= refsum * t1;
                h__[k + 2 + (k + 1) * h_dim1] -= refsum * t2;
                h__[k + 3 + (k + 1) * h_dim1] -= refsum * t3;
                if (k < *ktop) {
                    goto L85;
                }
                if (h__[k + 1 + k * h_dim1] != 0.) {
                    tst1 = (d__1 = h__[k + k * h_dim1], abs(d__1)) +
                           (d__2 = h__[k + 1 + (k + 1) * h_dim1], abs(d__2));
                    if (tst1 == 0.) {
                        if (k >= *ktop + 1) {
                            tst1 += (d__1 = h__[k + (k - 1) * h_dim1], abs(d__1));
                        }
                        if (k >= *ktop + 2) {
                            tst1 += (d__1 = h__[k + (k - 2) * h_dim1], abs(d__1));
                        }
                        if (k >= *ktop + 3) {
                            tst1 += (d__1 = h__[k + (k - 3) * h_dim1], abs(d__1));
                        }
                        if (k <= *kbot - 2) {
                            tst1 += (d__1 = h__[k + 2 + (k + 1) * h_dim1], abs(d__1));
                        }
                        if (k <= *kbot - 3) {
                            tst1 += (d__1 = h__[k + 3 + (k + 1) * h_dim1], abs(d__1));
                        }
                        if (k <= *kbot - 4) {
                            tst1 += (d__1 = h__[k + 4 + (k + 1) * h_dim1], abs(d__1));
                        }
                    }
                    d__2 = smlnum, d__3 = ulp * tst1;
                    if ((d__1 = h__[k + 1 + k * h_dim1], abs(d__1)) <= max(d__2, d__3)) {
                        d__3 = (d__1 = h__[k + 1 + k * h_dim1], abs(d__1)),
                        d__4 = (d__2 = h__[k + (k + 1) * h_dim1], abs(d__2));
                        h12 = max(d__3, d__4);
                        d__3 = (d__1 = h__[k + 1 + k * h_dim1], abs(d__1)),
                        d__4 = (d__2 = h__[k + (k + 1) * h_dim1], abs(d__2));
                        h21 = min(d__3, d__4);
                        d__3 = (d__1 = h__[k + 1 + (k + 1) * h_dim1], abs(d__1)),
                        d__4 =
                            (d__2 = h__[k + k * h_dim1] - h__[k + 1 + (k + 1) * h_dim1], abs(d__2));
                        h11 = max(d__3, d__4);
                        d__3 = (d__1 = h__[k + 1 + (k + 1) * h_dim1], abs(d__1)),
                        d__4 =
                            (d__2 = h__[k + k * h_dim1] - h__[k + 1 + (k + 1) * h_dim1], abs(d__2));
                        h22 = min(d__3, d__4);
                        scl = h11 + h12;
                        tst2 = h22 * (h11 / scl);
                        d__1 = smlnum, d__2 = ulp * tst2;
                        if (tst2 == 0. || h21 * (h12 / scl) <= max(d__1, d__2)) {
                            h__[k + 1 + k * h_dim1] = 0.;
                        }
                    }
                }
            L85:;
            }
            if (accum) {
                jbot = min(ndcol, *kbot);
            } else if (*wantt) {
                jbot = *n;
            } else {
                jbot = *kbot;
            }
            i__6 = mtop;
            for (m = mbot; m >= i__6; --m) {
                k = krcol + (m - 1 << 1);
                t1 = v[m * v_dim1 + 1];
                t2 = t1 * v[m * v_dim1 + 2];
                t3 = t1 * v[m * v_dim1 + 3];
                i__4 = *ktop, i__5 = krcol + (m << 1);
                i__7 = jbot;
                for (j = max(i__4, i__5); j <= i__7; ++j) {
                    refsum = h__[k + 1 + j * h_dim1] + v[m * v_dim1 + 2] * h__[k + 2 + j * h_dim1] +
                             v[m * v_dim1 + 3] * h__[k + 3 + j * h_dim1];
                    h__[k + 1 + j * h_dim1] -= refsum * t1;
                    h__[k + 2 + j * h_dim1] -= refsum * t2;
                    h__[k + 3 + j * h_dim1] -= refsum * t3;
                }
            }
            if (accum) {
                i__6 = mtop;
                for (m = mbot; m >= i__6; --m) {
                    k = krcol + (m - 1 << 1);
                    kms = k - incol;
                    i__7 = 1, i__4 = *ktop - incol;
                    i2 = max(i__7, i__4);
                    i__7 = i2, i__4 = kms - (krcol - incol) + 1;
                    i2 = max(i__7, i__4);
                    i__7 = kdu, i__4 = krcol + (mbot - 1 << 1) - incol + 5;
                    i4 = min(i__7, i__4);
                    t1 = v[m * v_dim1 + 1];
                    t2 = t1 * v[m * v_dim1 + 2];
                    t3 = t1 * v[m * v_dim1 + 3];
                    i__7 = i4;
                    for (j = i2; j <= i__7; ++j) {
                        refsum = u[j + (kms + 1) * u_dim1] +
                                 v[m * v_dim1 + 2] * u[j + (kms + 2) * u_dim1] +
                                 v[m * v_dim1 + 3] * u[j + (kms + 3) * u_dim1];
                        u[j + (kms + 1) * u_dim1] -= refsum * t1;
                        u[j + (kms + 2) * u_dim1] -= refsum * t2;
                        u[j + (kms + 3) * u_dim1] -= refsum * t3;
                    }
                }
            } else if (*wantz) {
                i__6 = mtop;
                for (m = mbot; m >= i__6; --m) {
                    k = krcol + (m - 1 << 1);
                    t1 = v[m * v_dim1 + 1];
                    t2 = t1 * v[m * v_dim1 + 2];
                    t3 = t1 * v[m * v_dim1 + 3];
                    i__7 = *ihiz;
                    for (j = *iloz; j <= i__7; ++j) {
                        refsum = z__[j + (k + 1) * z_dim1] +
                                 v[m * v_dim1 + 2] * z__[j + (k + 2) * z_dim1] +
                                 v[m * v_dim1 + 3] * z__[j + (k + 3) * z_dim1];
                        z__[j + (k + 1) * z_dim1] -= refsum * t1;
                        z__[j + (k + 2) * z_dim1] -= refsum * t2;
                        z__[j + (k + 3) * z_dim1] -= refsum * t3;
                    }
                }
            }
        }
        if (accum) {
            if (*wantt) {
                jtop = 1;
                jbot = *n;
            } else {
                jtop = *ktop;
                jbot = *kbot;
            }
            i__3 = 1, i__6 = *ktop - incol;
            k1 = max(i__3, i__6);
            i__3 = 0, i__6 = ndcol - *kbot;
            nu = kdu - max(i__3, i__6) - k1 + 1;
            i__3 = jbot;
            i__6 = *nh;
            for (jcol = min(ndcol, *kbot) + 1; i__6 < 0 ? jcol >= i__3 : jcol <= i__3;
                 jcol += i__6) {
                i__7 = *nh, i__4 = jbot - jcol + 1;
                jlen = min(i__7, i__4);
                dgemm_((char *)"C", (char *)"N", &nu, &jlen, &nu, &c_b8, &u[k1 + k1 * u_dim1], ldu,
                       &h__[incol + k1 + jcol * h_dim1], ldh, &c_b7, &wh[wh_offset], ldwh,
                       (ftnlen)1, (ftnlen)1);
                dlacpy_((char *)"A", &nu, &jlen, &wh[wh_offset], ldwh, &h__[incol + k1 + jcol * h_dim1],
                        ldh, (ftnlen)1);
            }
            i__6 = max(*ktop, incol) - 1;
            i__3 = *nv;
            for (jrow = jtop; i__3 < 0 ? jrow >= i__6 : jrow <= i__6; jrow += i__3) {
                i__7 = *nv, i__4 = max(*ktop, incol) - jrow;
                jlen = min(i__7, i__4);
                dgemm_((char *)"N", (char *)"N", &jlen, &nu, &nu, &c_b8, &h__[jrow + (incol + k1) * h_dim1], ldh,
                       &u[k1 + k1 * u_dim1], ldu, &c_b7, &wv[wv_offset], ldwv, (ftnlen)1,
                       (ftnlen)1);
                dlacpy_((char *)"A", &jlen, &nu, &wv[wv_offset], ldwv, &h__[jrow + (incol + k1) * h_dim1],
                        ldh, (ftnlen)1);
            }
            if (*wantz) {
                i__3 = *ihiz;
                i__6 = *nv;
                for (jrow = *iloz; i__6 < 0 ? jrow >= i__3 : jrow <= i__3; jrow += i__6) {
                    i__7 = *nv, i__4 = *ihiz - jrow + 1;
                    jlen = min(i__7, i__4);
                    dgemm_((char *)"N", (char *)"N", &jlen, &nu, &nu, &c_b8, &z__[jrow + (incol + k1) * z_dim1],
                           ldz, &u[k1 + k1 * u_dim1], ldu, &c_b7, &wv[wv_offset], ldwv, (ftnlen)1,
                           (ftnlen)1);
                    dlacpy_((char *)"A", &jlen, &nu, &wv[wv_offset], ldwv,
                            &z__[jrow + (incol + k1) * z_dim1], ldz, (ftnlen)1);
                }
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
