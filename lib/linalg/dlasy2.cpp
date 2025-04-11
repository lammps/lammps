#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__4 = 4;
static integer c__1 = 1;
static integer c__16 = 16;
static integer c__0 = 0;
int dlasy2_(logical *ltranl, logical *ltranr, integer *isgn, integer *n1, integer *n2,
            doublereal *tl, integer *ldtl, doublereal *tr, integer *ldtr, doublereal *b,
            integer *ldb, doublereal *scale, doublereal *x, integer *ldx, doublereal *xnorm,
            integer *info)
{
    static integer locu12[4] = {3, 4, 1, 2};
    static integer locl21[4] = {2, 1, 4, 3};
    static integer locu22[4] = {4, 3, 2, 1};
    static logical xswpiv[4] = {FALSE_, FALSE_, TRUE_, TRUE_};
    static logical bswpiv[4] = {FALSE_, TRUE_, FALSE_, TRUE_};
    integer b_dim1, b_offset, tl_dim1, tl_offset, tr_dim1, tr_offset, x_dim1, x_offset;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;
    integer i__, j, k;
    doublereal x2[2], l21, u11, u12;
    integer ip, jp;
    doublereal u22, t16[16], gam, bet, eps, sgn, tmp[4], tau1, btmp[4], smin;
    integer ipiv;
    doublereal temp;
    integer jpiv[4];
    doublereal xmax;
    integer ipsv, jpsv;
    logical bswap;
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *),
        dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
    logical xswap;
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    doublereal smlnum;
    tl_dim1 = *ldtl;
    tl_offset = 1 + tl_dim1;
    tl -= tl_offset;
    tr_dim1 = *ldtr;
    tr_offset = 1 + tr_dim1;
    tr -= tr_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    *info = 0;
    if (*n1 == 0 || *n2 == 0) {
        return 0;
    }
    eps = dlamch_((char *)"P", (ftnlen)1);
    smlnum = dlamch_((char *)"S", (ftnlen)1) / eps;
    sgn = (doublereal)(*isgn);
    k = *n1 + *n1 + *n2 - 2;
    switch (k) {
        case 1:
            goto L10;
        case 2:
            goto L20;
        case 3:
            goto L30;
        case 4:
            goto L50;
    }
L10:
    tau1 = tl[tl_dim1 + 1] + sgn * tr[tr_dim1 + 1];
    bet = abs(tau1);
    if (bet <= smlnum) {
        tau1 = smlnum;
        bet = smlnum;
        *info = 1;
    }
    *scale = 1.;
    gam = (d__1 = b[b_dim1 + 1], abs(d__1));
    if (smlnum * gam > bet) {
        *scale = 1. / gam;
    }
    x[x_dim1 + 1] = b[b_dim1 + 1] * *scale / tau1;
    *xnorm = (d__1 = x[x_dim1 + 1], abs(d__1));
    return 0;
L20:
    d__7 = (d__1 = tl[tl_dim1 + 1], abs(d__1)), d__8 = (d__2 = tr[tr_dim1 + 1], abs(d__2)),
    d__7 = max(d__7, d__8), d__8 = (d__3 = tr[(tr_dim1 << 1) + 1], abs(d__3)),
    d__7 = max(d__7, d__8), d__8 = (d__4 = tr[tr_dim1 + 2], abs(d__4)), d__7 = max(d__7, d__8),
    d__8 = (d__5 = tr[(tr_dim1 << 1) + 2], abs(d__5));
    d__6 = eps * max(d__7, d__8);
    smin = max(d__6, smlnum);
    tmp[0] = tl[tl_dim1 + 1] + sgn * tr[tr_dim1 + 1];
    tmp[3] = tl[tl_dim1 + 1] + sgn * tr[(tr_dim1 << 1) + 2];
    if (*ltranr) {
        tmp[1] = sgn * tr[tr_dim1 + 2];
        tmp[2] = sgn * tr[(tr_dim1 << 1) + 1];
    } else {
        tmp[1] = sgn * tr[(tr_dim1 << 1) + 1];
        tmp[2] = sgn * tr[tr_dim1 + 2];
    }
    btmp[0] = b[b_dim1 + 1];
    btmp[1] = b[(b_dim1 << 1) + 1];
    goto L40;
L30:
    d__7 = (d__1 = tr[tr_dim1 + 1], abs(d__1)), d__8 = (d__2 = tl[tl_dim1 + 1], abs(d__2)),
    d__7 = max(d__7, d__8), d__8 = (d__3 = tl[(tl_dim1 << 1) + 1], abs(d__3)),
    d__7 = max(d__7, d__8), d__8 = (d__4 = tl[tl_dim1 + 2], abs(d__4)), d__7 = max(d__7, d__8),
    d__8 = (d__5 = tl[(tl_dim1 << 1) + 2], abs(d__5));
    d__6 = eps * max(d__7, d__8);
    smin = max(d__6, smlnum);
    tmp[0] = tl[tl_dim1 + 1] + sgn * tr[tr_dim1 + 1];
    tmp[3] = tl[(tl_dim1 << 1) + 2] + sgn * tr[tr_dim1 + 1];
    if (*ltranl) {
        tmp[1] = tl[(tl_dim1 << 1) + 1];
        tmp[2] = tl[tl_dim1 + 2];
    } else {
        tmp[1] = tl[tl_dim1 + 2];
        tmp[2] = tl[(tl_dim1 << 1) + 1];
    }
    btmp[0] = b[b_dim1 + 1];
    btmp[1] = b[b_dim1 + 2];
L40:
    ipiv = idamax_(&c__4, tmp, &c__1);
    u11 = tmp[ipiv - 1];
    if (abs(u11) <= smin) {
        *info = 1;
        u11 = smin;
    }
    u12 = tmp[locu12[ipiv - 1] - 1];
    l21 = tmp[locl21[ipiv - 1] - 1] / u11;
    u22 = tmp[locu22[ipiv - 1] - 1] - u12 * l21;
    xswap = xswpiv[ipiv - 1];
    bswap = bswpiv[ipiv - 1];
    if (abs(u22) <= smin) {
        *info = 1;
        u22 = smin;
    }
    if (bswap) {
        temp = btmp[1];
        btmp[1] = btmp[0] - l21 * temp;
        btmp[0] = temp;
    } else {
        btmp[1] -= l21 * btmp[0];
    }
    *scale = 1.;
    if (smlnum * 2. * abs(btmp[1]) > abs(u22) || smlnum * 2. * abs(btmp[0]) > abs(u11)) {
        d__1 = abs(btmp[0]), d__2 = abs(btmp[1]);
        *scale = .5 / max(d__1, d__2);
        btmp[0] *= *scale;
        btmp[1] *= *scale;
    }
    x2[1] = btmp[1] / u22;
    x2[0] = btmp[0] / u11 - u12 / u11 * x2[1];
    if (xswap) {
        temp = x2[1];
        x2[1] = x2[0];
        x2[0] = temp;
    }
    x[x_dim1 + 1] = x2[0];
    if (*n1 == 1) {
        x[(x_dim1 << 1) + 1] = x2[1];
        *xnorm = (d__1 = x[x_dim1 + 1], abs(d__1)) + (d__2 = x[(x_dim1 << 1) + 1], abs(d__2));
    } else {
        x[x_dim1 + 2] = x2[1];
        d__3 = (d__1 = x[x_dim1 + 1], abs(d__1)), d__4 = (d__2 = x[x_dim1 + 2], abs(d__2));
        *xnorm = max(d__3, d__4);
    }
    return 0;
L50:
    d__5 = (d__1 = tr[tr_dim1 + 1], abs(d__1)), d__6 = (d__2 = tr[(tr_dim1 << 1) + 1], abs(d__2)),
    d__5 = max(d__5, d__6), d__6 = (d__3 = tr[tr_dim1 + 2], abs(d__3)), d__5 = max(d__5, d__6),
    d__6 = (d__4 = tr[(tr_dim1 << 1) + 2], abs(d__4));
    smin = max(d__5, d__6);
    d__5 = smin, d__6 = (d__1 = tl[tl_dim1 + 1], abs(d__1)), d__5 = max(d__5, d__6),
    d__6 = (d__2 = tl[(tl_dim1 << 1) + 1], abs(d__2)), d__5 = max(d__5, d__6),
    d__6 = (d__3 = tl[tl_dim1 + 2], abs(d__3)), d__5 = max(d__5, d__6),
    d__6 = (d__4 = tl[(tl_dim1 << 1) + 2], abs(d__4));
    smin = max(d__5, d__6);
    d__1 = eps * smin;
    smin = max(d__1, smlnum);
    btmp[0] = 0.;
    dcopy_(&c__16, btmp, &c__0, t16, &c__1);
    t16[0] = tl[tl_dim1 + 1] + sgn * tr[tr_dim1 + 1];
    t16[5] = tl[(tl_dim1 << 1) + 2] + sgn * tr[tr_dim1 + 1];
    t16[10] = tl[tl_dim1 + 1] + sgn * tr[(tr_dim1 << 1) + 2];
    t16[15] = tl[(tl_dim1 << 1) + 2] + sgn * tr[(tr_dim1 << 1) + 2];
    if (*ltranl) {
        t16[4] = tl[tl_dim1 + 2];
        t16[1] = tl[(tl_dim1 << 1) + 1];
        t16[14] = tl[tl_dim1 + 2];
        t16[11] = tl[(tl_dim1 << 1) + 1];
    } else {
        t16[4] = tl[(tl_dim1 << 1) + 1];
        t16[1] = tl[tl_dim1 + 2];
        t16[14] = tl[(tl_dim1 << 1) + 1];
        t16[11] = tl[tl_dim1 + 2];
    }
    if (*ltranr) {
        t16[8] = sgn * tr[(tr_dim1 << 1) + 1];
        t16[13] = sgn * tr[(tr_dim1 << 1) + 1];
        t16[2] = sgn * tr[tr_dim1 + 2];
        t16[7] = sgn * tr[tr_dim1 + 2];
    } else {
        t16[8] = sgn * tr[tr_dim1 + 2];
        t16[13] = sgn * tr[tr_dim1 + 2];
        t16[2] = sgn * tr[(tr_dim1 << 1) + 1];
        t16[7] = sgn * tr[(tr_dim1 << 1) + 1];
    }
    btmp[0] = b[b_dim1 + 1];
    btmp[1] = b[b_dim1 + 2];
    btmp[2] = b[(b_dim1 << 1) + 1];
    btmp[3] = b[(b_dim1 << 1) + 2];
    for (i__ = 1; i__ <= 3; ++i__) {
        xmax = 0.;
        for (ip = i__; ip <= 4; ++ip) {
            for (jp = i__; jp <= 4; ++jp) {
                if ((d__1 = t16[ip + (jp << 2) - 5], abs(d__1)) >= xmax) {
                    xmax = (d__1 = t16[ip + (jp << 2) - 5], abs(d__1));
                    ipsv = ip;
                    jpsv = jp;
                }
            }
        }
        if (ipsv != i__) {
            dswap_(&c__4, &t16[ipsv - 1], &c__4, &t16[i__ - 1], &c__4);
            temp = btmp[i__ - 1];
            btmp[i__ - 1] = btmp[ipsv - 1];
            btmp[ipsv - 1] = temp;
        }
        if (jpsv != i__) {
            dswap_(&c__4, &t16[(jpsv << 2) - 4], &c__1, &t16[(i__ << 2) - 4], &c__1);
        }
        jpiv[i__ - 1] = jpsv;
        if ((d__1 = t16[i__ + (i__ << 2) - 5], abs(d__1)) < smin) {
            *info = 1;
            t16[i__ + (i__ << 2) - 5] = smin;
        }
        for (j = i__ + 1; j <= 4; ++j) {
            t16[j + (i__ << 2) - 5] /= t16[i__ + (i__ << 2) - 5];
            btmp[j - 1] -= t16[j + (i__ << 2) - 5] * btmp[i__ - 1];
            for (k = i__ + 1; k <= 4; ++k) {
                t16[j + (k << 2) - 5] -= t16[j + (i__ << 2) - 5] * t16[i__ + (k << 2) - 5];
            }
        }
    }
    if (abs(t16[15]) < smin) {
        *info = 1;
        t16[15] = smin;
    }
    *scale = 1.;
    if (smlnum * 8. * abs(btmp[0]) > abs(t16[0]) || smlnum * 8. * abs(btmp[1]) > abs(t16[5]) ||
        smlnum * 8. * abs(btmp[2]) > abs(t16[10]) || smlnum * 8. * abs(btmp[3]) > abs(t16[15])) {
        d__1 = abs(btmp[0]), d__2 = abs(btmp[1]), d__1 = max(d__1, d__2), d__2 = abs(btmp[2]),
        d__1 = max(d__1, d__2), d__2 = abs(btmp[3]);
        *scale = .125 / max(d__1, d__2);
        btmp[0] *= *scale;
        btmp[1] *= *scale;
        btmp[2] *= *scale;
        btmp[3] *= *scale;
    }
    for (i__ = 1; i__ <= 4; ++i__) {
        k = 5 - i__;
        temp = 1. / t16[k + (k << 2) - 5];
        tmp[k - 1] = btmp[k - 1] * temp;
        for (j = k + 1; j <= 4; ++j) {
            tmp[k - 1] -= temp * t16[k + (j << 2) - 5] * tmp[j - 1];
        }
    }
    for (i__ = 1; i__ <= 3; ++i__) {
        if (jpiv[4 - i__ - 1] != 4 - i__) {
            temp = tmp[4 - i__ - 1];
            tmp[4 - i__ - 1] = tmp[jpiv[4 - i__ - 1] - 1];
            tmp[jpiv[4 - i__ - 1] - 1] = temp;
        }
    }
    x[x_dim1 + 1] = tmp[0];
    x[x_dim1 + 2] = tmp[1];
    x[(x_dim1 << 1) + 1] = tmp[2];
    x[(x_dim1 << 1) + 2] = tmp[3];
    d__1 = abs(tmp[0]) + abs(tmp[2]), d__2 = abs(tmp[1]) + abs(tmp[3]);
    *xnorm = max(d__1, d__2);
    return 0;
}
#ifdef __cplusplus
}
#endif
