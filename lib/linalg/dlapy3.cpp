#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
doublereal dlapy3_(doublereal *x, doublereal *y, doublereal *z__)
{
    doublereal ret_val, d__1, d__2, d__3;
    double sqrt(doublereal);
    doublereal w, xabs, yabs, zabs;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal hugeval;
    hugeval = dlamch_((char *)"Overflow", (ftnlen)8);
    xabs = abs(*x);
    yabs = abs(*y);
    zabs = abs(*z__);
    d__1 = max(xabs, yabs);
    w = max(d__1, zabs);
    if (w == 0. || w > hugeval) {
        ret_val = xabs + yabs + zabs;
    } else {
        d__1 = xabs / w;
        d__2 = yabs / w;
        d__3 = zabs / w;
        ret_val = w * sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    }
    return ret_val;
}
#ifdef __cplusplus
}
#endif
