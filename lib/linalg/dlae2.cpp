#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlae2_(doublereal *a, doublereal *b, doublereal *c__, doublereal *rt1, doublereal *rt2)
{
    doublereal d__1;
    double sqrt(doublereal);
    doublereal ab, df, tb, sm, rt, adf, acmn, acmx;
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
        *rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else if (sm > 0.) {
        *rt1 = (sm + rt) * .5;
        *rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else {
        *rt1 = rt * .5;
        *rt2 = rt * -.5;
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
