#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
doublereal ddot_(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy)
{
    integer i__1;
    doublereal ret_val;
    integer i__, m, ix, iy, mp1;
    doublereal dtemp;
    --dy;
    --dx;
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
        return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
        m = *n % 5;
        if (m != 0) {
            i__1 = m;
            for (i__ = 1; i__ <= i__1; ++i__) {
                dtemp += dx[i__] * dy[i__];
            }
            if (*n < 5) {
                ret_val = dtemp;
                return ret_val;
            }
        }
        mp1 = m + 1;
        i__1 = *n;
        for (i__ = mp1; i__ <= i__1; i__ += 5) {
            dtemp = dtemp + dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1] +
                    dx[i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3] +
                    dx[i__ + 4] * dy[i__ + 4];
        }
    } else {
        ix = 1;
        iy = 1;
        if (*incx < 0) {
            ix = (-(*n) + 1) * *incx + 1;
        }
        if (*incy < 0) {
            iy = (-(*n) + 1) * *incy + 1;
        }
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            dtemp += dx[ix] * dy[iy];
            ix += *incx;
            iy += *incy;
        }
    }
    ret_val = dtemp;
    return ret_val;
}
#ifdef __cplusplus
}
#endif
