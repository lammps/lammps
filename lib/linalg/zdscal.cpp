#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int zdscal_(integer *n, doublereal *da, doublecomplex *zx, integer *incx)
{
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    doublecomplex z__1;
    double d_lmp_imag(doublecomplex *);
    integer i__, nincx;
    --zx;
    if (*n <= 0 || *incx <= 0 || *da == 1.) {
        return 0;
    }
    if (*incx == 1) {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            i__2 = i__;
            i__3 = i__;
            d__1 = *da * zx[i__3].r;
            d__2 = *da * d_lmp_imag(&zx[i__]);
            z__1.r = d__1, z__1.i = d__2;
            zx[i__2].r = z__1.r, zx[i__2].i = z__1.i;
        }
    } else {
        nincx = *n * *incx;
        i__1 = nincx;
        i__2 = *incx;
        for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
            i__3 = i__;
            i__4 = i__;
            d__1 = *da * zx[i__4].r;
            d__2 = *da * d_lmp_imag(&zx[i__]);
            z__1.r = d__1, z__1.i = d__2;
            zx[i__3].r = z__1.r, zx[i__3].i = z__1.i;
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
