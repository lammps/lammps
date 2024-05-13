#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int zlacgv_(integer *n, doublecomplex *x, integer *incx)
{
    integer i__1, i__2;
    doublecomplex z__1;
    void d_lmp_cnjg(doublecomplex *, doublecomplex *);
    integer i__, ioff;
    --x;
    if (*incx == 1) {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            i__2 = i__;
            d_lmp_cnjg(&z__1, &x[i__]);
            x[i__2].r = z__1.r, x[i__2].i = z__1.i;
        }
    } else {
        ioff = 1;
        if (*incx < 0) {
            ioff = 1 - (*n - 1) * *incx;
        }
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            i__2 = ioff;
            d_lmp_cnjg(&z__1, &x[ioff]);
            x[i__2].r = z__1.r, x[i__2].i = z__1.i;
            ioff += *incx;
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
