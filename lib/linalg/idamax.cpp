#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
integer idamax_(integer *n, doublereal *dx, integer *incx)
{
    integer ret_val, i__1;
    doublereal d__1;
    integer i__, ix;
    doublereal dmax__;
    --dx;
    ret_val = 0;
    if (*n < 1 || *incx <= 0) {
        return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
        return ret_val;
    }
    if (*incx == 1) {
        dmax__ = abs(dx[1]);
        i__1 = *n;
        for (i__ = 2; i__ <= i__1; ++i__) {
            if ((d__1 = dx[i__], abs(d__1)) > dmax__) {
                ret_val = i__;
                dmax__ = (d__1 = dx[i__], abs(d__1));
            }
        }
    } else {
        ix = 1;
        dmax__ = abs(dx[1]);
        ix += *incx;
        i__1 = *n;
        for (i__ = 2; i__ <= i__1; ++i__) {
            if ((d__1 = dx[ix], abs(d__1)) > dmax__) {
                ret_val = i__;
                dmax__ = (d__1 = dx[ix], abs(d__1));
            }
            ix += *incx;
        }
    }
    return ret_val;
}
#ifdef __cplusplus
}
#endif
