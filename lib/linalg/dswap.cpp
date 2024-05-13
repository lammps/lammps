#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dswap_(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy)
{
    integer i__1;
    integer i__, m, ix, iy, mp1;
    doublereal dtemp;
    --dy;
    --dx;
    if (*n <= 0) {
        return 0;
    }
    if (*incx == 1 && *incy == 1) {
        m = *n % 3;
        if (m != 0) {
            i__1 = m;
            for (i__ = 1; i__ <= i__1; ++i__) {
                dtemp = dx[i__];
                dx[i__] = dy[i__];
                dy[i__] = dtemp;
            }
            if (*n < 3) {
                return 0;
            }
        }
        mp1 = m + 1;
        i__1 = *n;
        for (i__ = mp1; i__ <= i__1; i__ += 3) {
            dtemp = dx[i__];
            dx[i__] = dy[i__];
            dy[i__] = dtemp;
            dtemp = dx[i__ + 1];
            dx[i__ + 1] = dy[i__ + 1];
            dy[i__ + 1] = dtemp;
            dtemp = dx[i__ + 2];
            dx[i__ + 2] = dy[i__ + 2];
            dy[i__ + 2] = dtemp;
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
            dtemp = dx[ix];
            dx[ix] = dy[iy];
            dy[iy] = dtemp;
            ix += *incx;
            iy += *incy;
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
