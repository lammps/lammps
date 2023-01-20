#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlartg_(doublereal *f, doublereal *g, doublereal *cs, doublereal *sn, doublereal *r__)
{
    integer i__1;
    doublereal d__1, d__2;
    double log(doublereal), pow_lmp_di(doublereal *, integer *), sqrt(doublereal);
    integer i__;
    doublereal f1, g1, eps, scale;
    integer count;
    doublereal safmn2, safmx2;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal safmin;
    safmin = dlamch_((char *)"S", (ftnlen)1);
    eps = dlamch_((char *)"E", (ftnlen)1);
    d__1 = dlamch_((char *)"B", (ftnlen)1);
    i__1 = (integer)(log(safmin / eps) / log(dlamch_((char *)"B", (ftnlen)1)) / 2.);
    safmn2 = pow_lmp_di(&d__1, &i__1);
    safmx2 = 1. / safmn2;
    if (*g == 0.) {
        *cs = 1.;
        *sn = 0.;
        *r__ = *f;
    } else if (*f == 0.) {
        *cs = 0.;
        *sn = 1.;
        *r__ = *g;
    } else {
        f1 = *f;
        g1 = *g;
        d__1 = abs(f1), d__2 = abs(g1);
        scale = max(d__1, d__2);
        if (scale >= safmx2) {
            count = 0;
        L10:
            ++count;
            f1 *= safmn2;
            g1 *= safmn2;
            d__1 = abs(f1), d__2 = abs(g1);
            scale = max(d__1, d__2);
            if (scale >= safmx2) {
                goto L10;
            }
            d__1 = f1;
            d__2 = g1;
            *r__ = sqrt(d__1 * d__1 + d__2 * d__2);
            *cs = f1 / *r__;
            *sn = g1 / *r__;
            i__1 = count;
            for (i__ = 1; i__ <= i__1; ++i__) {
                *r__ *= safmx2;
            }
        } else if (scale <= safmn2) {
            count = 0;
        L30:
            ++count;
            f1 *= safmx2;
            g1 *= safmx2;
            d__1 = abs(f1), d__2 = abs(g1);
            scale = max(d__1, d__2);
            if (scale <= safmn2) {
                goto L30;
            }
            d__1 = f1;
            d__2 = g1;
            *r__ = sqrt(d__1 * d__1 + d__2 * d__2);
            *cs = f1 / *r__;
            *sn = g1 / *r__;
            i__1 = count;
            for (i__ = 1; i__ <= i__1; ++i__) {
                *r__ *= safmn2;
            }
        } else {
            d__1 = f1;
            d__2 = g1;
            *r__ = sqrt(d__1 * d__1 + d__2 * d__2);
            *cs = f1 / *r__;
            *sn = g1 / *r__;
        }
        if (abs(*f) > abs(*g) && *cs < 0.) {
            *cs = -(*cs);
            *sn = -(*sn);
            *r__ = -(*r__);
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
