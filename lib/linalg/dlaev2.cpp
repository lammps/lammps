#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlaev2_(doublereal *a, doublereal *b, doublereal *c__, doublereal *rt1, doublereal *rt2,
            doublereal *cs1, doublereal *sn1)
{
    doublereal d__1;
    double sqrt(doublereal);
    doublereal ab, df, cs, ct, tb, sm, tn, rt, adf, acs;
    integer sgn1, sgn2;
    doublereal acmn, acmx;
    sm = *a + *c__;
    df = *a - *c__;
    adf = abs(df);
    tb = *b + *b;
    ab = abs(tb);
    if (abs(*a) > abs(*c__)) {
        acmx = *a;
        acmn = *c__;
    } else {
        acmx = *c__;
        acmn = *a;
    }
    if (adf > ab) {
        d__1 = ab / adf;
        rt = adf * sqrt(d__1 * d__1 + 1.);
    } else if (adf < ab) {
        d__1 = adf / ab;
        rt = ab * sqrt(d__1 * d__1 + 1.);
    } else {
        rt = ab * sqrt(2.);
    }
    if (sm < 0.) {
        *rt1 = (sm - rt) * .5;
        sgn1 = -1;
        *rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else if (sm > 0.) {
        *rt1 = (sm + rt) * .5;
        sgn1 = 1;
        *rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else {
        *rt1 = rt * .5;
        *rt2 = rt * -.5;
        sgn1 = 1;
    }
    if (df >= 0.) {
        cs = df + rt;
        sgn2 = 1;
    } else {
        cs = df - rt;
        sgn2 = -1;
    }
    acs = abs(cs);
    if (acs > ab) {
        ct = -tb / cs;
        *sn1 = 1. / sqrt(ct * ct + 1.);
        *cs1 = ct * *sn1;
    } else {
        if (ab == 0.) {
            *cs1 = 1.;
            *sn1 = 0.;
        } else {
            tn = -cs / tb;
            *cs1 = 1. / sqrt(tn * tn + 1.);
            *sn1 = tn * *cs1;
        }
    }
    if (sgn1 == sgn2) {
        tn = *cs1;
        *cs1 = -(*sn1);
        *sn1 = tn;
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
