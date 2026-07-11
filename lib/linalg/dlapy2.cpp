#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
doublereal dlapy2_(doublereal *x, doublereal *y)
{
    doublereal ret_val, d__1;
    double sqrt(doublereal);
    logical x_is_nan__, y_is_nan__;
    doublereal w, z__, xabs, yabs;
    extern doublereal dlamch_(char *, ftnlen);
    extern logical disnan_(doublereal *);
    doublereal hugeval;
    x_is_nan__ = disnan_(x);
    y_is_nan__ = disnan_(y);
    if (x_is_nan__) {
        ret_val = *x;
    }
    if (y_is_nan__) {
        ret_val = *y;
    }
    hugeval = dlamch_((char *)"Overflow", (ftnlen)8);
    if (!(x_is_nan__ || y_is_nan__)) {
        xabs = abs(*x);
        yabs = abs(*y);
        w = max(xabs, yabs);
        z__ = min(xabs, yabs);
        if (z__ == 0. || w > hugeval) {
            ret_val = w;
        } else {
            d__1 = z__ / w;
            ret_val = w * sqrt(d__1 * d__1 + 1.);
        }
    }
    return ret_val;
}
#ifdef __cplusplus
}
#endif
