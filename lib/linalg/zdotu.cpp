#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
VOID zdotu_(doublecomplex *ret_val, integer *n, doublecomplex *zx, integer *incx, doublecomplex *zy,
            integer *incy)
{
    integer i__1, i__2, i__3;
    doublecomplex z__1, z__2;
    integer i__, ix, iy;
    doublecomplex ztemp;
    --zy;
    --zx;
    ztemp.r = 0., ztemp.i = 0.;
    ret_val->r = 0., ret_val->i = 0.;
    if (*n <= 0) {
        return;
    }
    if (*incx == 1 && *incy == 1) {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            i__2 = i__;
            i__3 = i__;
            z__2.r = zx[i__2].r * zy[i__3].r - zx[i__2].i * zy[i__3].i,
            z__2.i = zx[i__2].r * zy[i__3].i + zx[i__2].i * zy[i__3].r;
            z__1.r = ztemp.r + z__2.r, z__1.i = ztemp.i + z__2.i;
            ztemp.r = z__1.r, ztemp.i = z__1.i;
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
            i__3 = iy;
            z__2.r = zx[i__2].r * zy[i__3].r - zx[i__2].i * zy[i__3].i,
            z__2.i = zx[i__2].r * zy[i__3].i + zx[i__2].i * zy[i__3].r;
            z__1.r = ztemp.r + z__2.r, z__1.i = ztemp.i + z__2.i;
            ztemp.r = z__1.r, ztemp.i = z__1.i;
            ix += *incx;
            iy += *incy;
        }
    }
    ret_val->r = ztemp.r, ret_val->i = ztemp.i;
    return;
}
#ifdef __cplusplus
}
#endif
