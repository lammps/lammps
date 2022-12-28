#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int zscal_(integer *n, doublecomplex *za, doublecomplex *zx, integer *incx)
{
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1;
    integer i__, nincx;
    --zx;
    if (*n <= 0 || *incx <= 0 || za->r == 1. && za->i == 0.) {
        return 0;
    }
    if (*incx == 1) {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            i__2 = i__;
            i__3 = i__;
            z__1.r = za->r * zx[i__3].r - za->i * zx[i__3].i,
            z__1.i = za->r * zx[i__3].i + za->i * zx[i__3].r;
            zx[i__2].r = z__1.r, zx[i__2].i = z__1.i;
        }
    } else {
        nincx = *n * *incx;
        i__1 = nincx;
        i__2 = *incx;
        for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
            i__3 = i__;
            i__4 = i__;
            z__1.r = za->r * zx[i__4].r - za->i * zx[i__4].i,
            z__1.i = za->r * zx[i__4].i + za->i * zx[i__4].r;
            zx[i__3].r = z__1.r, zx[i__3].i = z__1.i;
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
