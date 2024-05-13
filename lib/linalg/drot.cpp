#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int drot_(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy, doublereal *c__,
          doublereal *s)
{
    integer i__1;
    integer i__, ix, iy;
    doublereal dtemp;
    --dy;
    --dx;
    if (*n <= 0) {
        return 0;
    }
    if (*incx == 1 && *incy == 1) {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            dtemp = *c__ * dx[i__] + *s * dy[i__];
            dy[i__] = *c__ * dy[i__] - *s * dx[i__];
            dx[i__] = dtemp;
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
            dtemp = *c__ * dx[ix] + *s * dy[iy];
            dy[iy] = *c__ * dy[iy] - *s * dx[ix];
            dx[ix] = dtemp;
            ix += *incx;
            iy += *incy;
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
