#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int zswap_(integer *n, doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy)
{
    integer i__1, i__2, i__3;
    integer i__, ix, iy;
    doublecomplex ztemp;
    --zy;
    --zx;
    if (*n <= 0) {
        return 0;
    }
    if (*incx == 1 && *incy == 1) {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            i__2 = i__;
            ztemp.r = zx[i__2].r, ztemp.i = zx[i__2].i;
            i__2 = i__;
            i__3 = i__;
            zx[i__2].r = zy[i__3].r, zx[i__2].i = zy[i__3].i;
            i__2 = i__;
            zy[i__2].r = ztemp.r, zy[i__2].i = ztemp.i;
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
            i__2 = ix;
            ztemp.r = zx[i__2].r, ztemp.i = zx[i__2].i;
            i__2 = ix;
            i__3 = iy;
            zx[i__2].r = zy[i__3].r, zx[i__2].i = zy[i__3].i;
            i__2 = iy;
            zy[i__2].r = ztemp.r, zy[i__2].i = ztemp.i;
            ix += *incx;
            iy += *incy;
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
