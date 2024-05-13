#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dladiv_(doublereal *a, doublereal *b, doublereal *c__, doublereal *d__, doublereal *p,
            doublereal *q)
{
    doublereal d__1, d__2;
    doublereal s, aa, ab, bb, cc, cd, dd, be, un, ov, eps;
    extern doublereal dlamch_(char *, ftnlen);
    extern int dladiv1_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                        doublereal *);
    aa = *a;
    bb = *b;
    cc = *c__;
    dd = *d__;
    d__1 = abs(*a), d__2 = abs(*b);
    ab = max(d__1, d__2);
    d__1 = abs(*c__), d__2 = abs(*d__);
    cd = max(d__1, d__2);
    s = 1.;
    ov = dlamch_((char *)"Overflow threshold", (ftnlen)18);
    un = dlamch_((char *)"Safe minimum", (ftnlen)12);
    eps = dlamch_((char *)"Epsilon", (ftnlen)7);
    be = 2. / (eps * eps);
    if (ab >= ov * .5) {
        aa *= .5;
        bb *= .5;
        s *= 2.;
    }
    if (cd >= ov * .5) {
        cc *= .5;
        dd *= .5;
        s *= .5;
    }
    if (ab <= un * 2. / eps) {
        aa *= be;
        bb *= be;
        s /= be;
    }
    if (cd <= un * 2. / eps) {
        cc *= be;
        dd *= be;
        s *= be;
    }
    if (abs(*d__) <= abs(*c__)) {
        dladiv1_(&aa, &bb, &cc, &dd, p, q);
    } else {
        dladiv1_(&bb, &aa, &dd, &cc, p, q);
        *q = -(*q);
    }
    *p *= s;
    *q *= s;
    return 0;
}
int dladiv1_(doublereal *a, doublereal *b, doublereal *c__, doublereal *d__, doublereal *p,
             doublereal *q)
{
    doublereal r__, t;
    extern doublereal dladiv2_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                               doublereal *);
    r__ = *d__ / *c__;
    t = 1. / (*c__ + *d__ * r__);
    *p = dladiv2_(a, b, c__, d__, &r__, &t);
    *a = -(*a);
    *q = dladiv2_(b, a, c__, d__, &r__, &t);
    return 0;
}
doublereal dladiv2_(doublereal *a, doublereal *b, doublereal *c__, doublereal *d__, doublereal *r__,
                    doublereal *t)
{
    doublereal ret_val;
    doublereal br;
    if (*r__ != 0.) {
        br = *b * *r__;
        if (br != 0.) {
            ret_val = (*a + br) * *t;
        } else {
            ret_val = *a * *t + *b * *t * *r__;
        }
    } else {
        ret_val = (*a + *d__ * (*b / *c__)) * *t;
    }
    return ret_val;
}
#ifdef __cplusplus
}
#endif
