#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
int dlahqr_(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, doublereal *h__,
            integer *ldh, doublereal *wr, doublereal *wi, integer *iloz, integer *ihiz,
            doublereal *z__, integer *ldz, integer *info)
{
    integer h_dim1, h_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;
    double sqrt(doublereal);
    integer i__, j, k, l, m;
    doublereal s, v[3];
    integer i1, i2;
    doublereal t1, t2, t3, v2, v3, aa, ab, ba, bb, h11, h12, h21, h22, cs;
    integer nh;
    doublereal sn;
    integer nr;
    doublereal tr;
    integer nz;
    doublereal det, h21s;
    integer its;
    doublereal ulp, sum, tst, rt1i, rt2i, rt1r, rt2r;
    extern int drot_(integer *, doublereal *, integer *, doublereal *, integer *, doublereal *,
                     doublereal *);
    integer kdefl;
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer itmax;
    extern int dlanv2_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                       doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern int dlarfg_(integer *, doublereal *, doublereal *, integer *, doublereal *);
    doublereal safmin, safmax, rtdisc, smlnum;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --wr;
    --wi;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    *info = 0;
    if (*n == 0) {
        return 0;
    }
    if (*ilo == *ihi) {
        wr[*ilo] = h__[*ilo + *ilo * h_dim1];
        wi[*ilo] = 0.;
        return 0;
    }
    i__1 = *ihi - 3;
    for (j = *ilo; j <= i__1; ++j) {
        h__[j + 2 + j * h_dim1] = 0.;
        h__[j + 3 + j * h_dim1] = 0.;
    }
    if (*ilo <= *ihi - 2) {
        h__[*ihi + (*ihi - 2) * h_dim1] = 0.;
    }
    nh = *ihi - *ilo + 1;
    nz = *ihiz - *iloz + 1;
    safmin = dlamch_((char *)"SAFE MINIMUM", (ftnlen)12);
    safmax = 1. / safmin;
    ulp = dlamch_((char *)"PRECISION", (ftnlen)9);
    smlnum = safmin * ((doublereal)nh / ulp);
    if (*wantt) {
        i1 = 1;
        i2 = *n;
    }
    itmax = max(10, nh) * 30;
    kdefl = 0;
    i__ = *ihi;
L20:
    l = *ilo;
    if (i__ < *ilo) {
        goto L160;
    }
    i__1 = itmax;
    for (its = 0; its <= i__1; ++its) {
        i__2 = l + 1;
        for (k = i__; k >= i__2; --k) {
            if ((d__1 = h__[k + (k - 1) * h_dim1], abs(d__1)) <= smlnum) {
                goto L40;
            }
            tst = (d__1 = h__[k - 1 + (k - 1) * h_dim1], abs(d__1)) +
                  (d__2 = h__[k + k * h_dim1], abs(d__2));
            if (tst == 0.) {
                if (k - 2 >= *ilo) {
                    tst += (d__1 = h__[k - 1 + (k - 2) * h_dim1], abs(d__1));
                }
                if (k + 1 <= *ihi) {
                    tst += (d__1 = h__[k + 1 + k * h_dim1], abs(d__1));
                }
            }
            if ((d__1 = h__[k + (k - 1) * h_dim1], abs(d__1)) <= ulp * tst) {
                d__3 = (d__1 = h__[k + (k - 1) * h_dim1], abs(d__1)),
                d__4 = (d__2 = h__[k - 1 + k * h_dim1], abs(d__2));
                ab = max(d__3, d__4);
                d__3 = (d__1 = h__[k + (k - 1) * h_dim1], abs(d__1)),
                d__4 = (d__2 = h__[k - 1 + k * h_dim1], abs(d__2));
                ba = min(d__3, d__4);
                d__3 = (d__1 = h__[k + k * h_dim1], abs(d__1)),
                d__4 = (d__2 = h__[k - 1 + (k - 1) * h_dim1] - h__[k + k * h_dim1], abs(d__2));
                aa = max(d__3, d__4);
                d__3 = (d__1 = h__[k + k * h_dim1], abs(d__1)),
                d__4 = (d__2 = h__[k - 1 + (k - 1) * h_dim1] - h__[k + k * h_dim1], abs(d__2));
                bb = min(d__3, d__4);
                s = aa + ab;
                d__1 = smlnum, d__2 = ulp * (bb * (aa / s));
                if (ba * (ab / s) <= max(d__1, d__2)) {
                    goto L40;
                }
            }
        }
    L40:
        l = k;
        if (l > *ilo) {
            h__[l + (l - 1) * h_dim1] = 0.;
        }
        if (l >= i__ - 1) {
            goto L150;
        }
        ++kdefl;
        if (!(*wantt)) {
            i1 = l;
            i2 = i__;
        }
        if (kdefl % 20 == 0) {
            s = (d__1 = h__[i__ + (i__ - 1) * h_dim1], abs(d__1)) +
                (d__2 = h__[i__ - 1 + (i__ - 2) * h_dim1], abs(d__2));
            h11 = s * .75 + h__[i__ + i__ * h_dim1];
            h12 = s * -.4375;
            h21 = s;
            h22 = h11;
        } else if (kdefl % 10 == 0) {
            s = (d__1 = h__[l + 1 + l * h_dim1], abs(d__1)) +
                (d__2 = h__[l + 2 + (l + 1) * h_dim1], abs(d__2));
            h11 = s * .75 + h__[l + l * h_dim1];
            h12 = s * -.4375;
            h21 = s;
            h22 = h11;
        } else {
            h11 = h__[i__ - 1 + (i__ - 1) * h_dim1];
            h21 = h__[i__ + (i__ - 1) * h_dim1];
            h12 = h__[i__ - 1 + i__ * h_dim1];
            h22 = h__[i__ + i__ * h_dim1];
        }
        s = abs(h11) + abs(h12) + abs(h21) + abs(h22);
        if (s == 0.) {
            rt1r = 0.;
            rt1i = 0.;
            rt2r = 0.;
            rt2i = 0.;
        } else {
            h11 /= s;
            h21 /= s;
            h12 /= s;
            h22 /= s;
            tr = (h11 + h22) / 2.;
            det = (h11 - tr) * (h22 - tr) - h12 * h21;
            rtdisc = sqrt((abs(det)));
            if (det >= 0.) {
                rt1r = tr * s;
                rt2r = rt1r;
                rt1i = rtdisc * s;
                rt2i = -rt1i;
            } else {
                rt1r = tr + rtdisc;
                rt2r = tr - rtdisc;
                if ((d__1 = rt1r - h22, abs(d__1)) <= (d__2 = rt2r - h22, abs(d__2))) {
                    rt1r *= s;
                    rt2r = rt1r;
                } else {
                    rt2r *= s;
                    rt1r = rt2r;
                }
                rt1i = 0.;
                rt2i = 0.;
            }
        }
        i__2 = l;
        for (m = i__ - 2; m >= i__2; --m) {
            h21s = h__[m + 1 + m * h_dim1];
            s = (d__1 = h__[m + m * h_dim1] - rt2r, abs(d__1)) + abs(rt2i) + abs(h21s);
            h21s = h__[m + 1 + m * h_dim1] / s;
            v[0] = h21s * h__[m + (m + 1) * h_dim1] +
                   (h__[m + m * h_dim1] - rt1r) * ((h__[m + m * h_dim1] - rt2r) / s) -
                   rt1i * (rt2i / s);
            v[1] = h21s * (h__[m + m * h_dim1] + h__[m + 1 + (m + 1) * h_dim1] - rt1r - rt2r);
            v[2] = h21s * h__[m + 2 + (m + 1) * h_dim1];
            s = abs(v[0]) + abs(v[1]) + abs(v[2]);
            v[0] /= s;
            v[1] /= s;
            v[2] /= s;
            if (m == l) {
                goto L60;
            }
            if ((d__1 = h__[m + (m - 1) * h_dim1], abs(d__1)) * (abs(v[1]) + abs(v[2])) <=
                ulp * abs(v[0]) *
                    ((d__2 = h__[m - 1 + (m - 1) * h_dim1], abs(d__2)) +
                     (d__3 = h__[m + m * h_dim1], abs(d__3)) +
                     (d__4 = h__[m + 1 + (m + 1) * h_dim1], abs(d__4)))) {
                goto L60;
            }
        }
    L60:
        i__2 = i__ - 1;
        for (k = m; k <= i__2; ++k) {
            i__3 = 3, i__4 = i__ - k + 1;
            nr = min(i__3, i__4);
            if (k > m) {
                dcopy_(&nr, &h__[k + (k - 1) * h_dim1], &c__1, v, &c__1);
            }
            dlarfg_(&nr, v, &v[1], &c__1, &t1);
            if (k > m) {
                h__[k + (k - 1) * h_dim1] = v[0];
                h__[k + 1 + (k - 1) * h_dim1] = 0.;
                if (k < i__ - 1) {
                    h__[k + 2 + (k - 1) * h_dim1] = 0.;
                }
            } else if (m > l) {
                h__[k + (k - 1) * h_dim1] *= 1. - t1;
            }
            v2 = v[1];
            t2 = t1 * v2;
            if (nr == 3) {
                v3 = v[2];
                t3 = t1 * v3;
                i__3 = i2;
                for (j = k; j <= i__3; ++j) {
                    sum = h__[k + j * h_dim1] + v2 * h__[k + 1 + j * h_dim1] +
                          v3 * h__[k + 2 + j * h_dim1];
                    h__[k + j * h_dim1] -= sum * t1;
                    h__[k + 1 + j * h_dim1] -= sum * t2;
                    h__[k + 2 + j * h_dim1] -= sum * t3;
                }
                i__4 = k + 3;
                i__3 = min(i__4, i__);
                for (j = i1; j <= i__3; ++j) {
                    sum = h__[j + k * h_dim1] + v2 * h__[j + (k + 1) * h_dim1] +
                          v3 * h__[j + (k + 2) * h_dim1];
                    h__[j + k * h_dim1] -= sum * t1;
                    h__[j + (k + 1) * h_dim1] -= sum * t2;
                    h__[j + (k + 2) * h_dim1] -= sum * t3;
                }
                if (*wantz) {
                    i__3 = *ihiz;
                    for (j = *iloz; j <= i__3; ++j) {
                        sum = z__[j + k * z_dim1] + v2 * z__[j + (k + 1) * z_dim1] +
                              v3 * z__[j + (k + 2) * z_dim1];
                        z__[j + k * z_dim1] -= sum * t1;
                        z__[j + (k + 1) * z_dim1] -= sum * t2;
                        z__[j + (k + 2) * z_dim1] -= sum * t3;
                    }
                }
            } else if (nr == 2) {
                i__3 = i2;
                for (j = k; j <= i__3; ++j) {
                    sum = h__[k + j * h_dim1] + v2 * h__[k + 1 + j * h_dim1];
                    h__[k + j * h_dim1] -= sum * t1;
                    h__[k + 1 + j * h_dim1] -= sum * t2;
                }
                i__3 = i__;
                for (j = i1; j <= i__3; ++j) {
                    sum = h__[j + k * h_dim1] + v2 * h__[j + (k + 1) * h_dim1];
                    h__[j + k * h_dim1] -= sum * t1;
                    h__[j + (k + 1) * h_dim1] -= sum * t2;
                }
                if (*wantz) {
                    i__3 = *ihiz;
                    for (j = *iloz; j <= i__3; ++j) {
                        sum = z__[j + k * z_dim1] + v2 * z__[j + (k + 1) * z_dim1];
                        z__[j + k * z_dim1] -= sum * t1;
                        z__[j + (k + 1) * z_dim1] -= sum * t2;
                    }
                }
            }
        }
    }
    *info = i__;
    return 0;
L150:
    if (l == i__) {
        wr[i__] = h__[i__ + i__ * h_dim1];
        wi[i__] = 0.;
    } else if (l == i__ - 1) {
        dlanv2_(&h__[i__ - 1 + (i__ - 1) * h_dim1], &h__[i__ - 1 + i__ * h_dim1],
                &h__[i__ + (i__ - 1) * h_dim1], &h__[i__ + i__ * h_dim1], &wr[i__ - 1],
                &wi[i__ - 1], &wr[i__], &wi[i__], &cs, &sn);
        if (*wantt) {
            if (i2 > i__) {
                i__1 = i2 - i__;
                drot_(&i__1, &h__[i__ - 1 + (i__ + 1) * h_dim1], ldh,
                      &h__[i__ + (i__ + 1) * h_dim1], ldh, &cs, &sn);
            }
            i__1 = i__ - i1 - 1;
            drot_(&i__1, &h__[i1 + (i__ - 1) * h_dim1], &c__1, &h__[i1 + i__ * h_dim1], &c__1, &cs,
                  &sn);
        }
        if (*wantz) {
            drot_(&nz, &z__[*iloz + (i__ - 1) * z_dim1], &c__1, &z__[*iloz + i__ * z_dim1], &c__1,
                  &cs, &sn);
        }
    }
    kdefl = 0;
    i__ = l - 1;
    goto L20;
L160:
    return 0;
}
#ifdef __cplusplus
}
#endif
