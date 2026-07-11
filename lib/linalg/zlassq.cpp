#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int zlassq_(integer *n, doublecomplex *x, integer *incx, doublereal *scale, doublereal *sumsq)
{
    integer i__1, i__2, i__3;
    doublereal d__1;
    double d_lmp_imag(doublecomplex *);
    integer ix;
    doublereal temp1;
    extern logical disnan_(doublereal *);
    --x;
    if (*n > 0) {
        i__1 = (*n - 1) * *incx + 1;
        i__2 = *incx;
        for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
            i__3 = ix;
            temp1 = (d__1 = x[i__3].r, abs(d__1));
            if (temp1 > 0. || disnan_(&temp1)) {
                if (*scale < temp1) {
                    d__1 = *scale / temp1;
                    *sumsq = *sumsq * (d__1 * d__1) + 1;
                    *scale = temp1;
                } else {
                    d__1 = temp1 / *scale;
                    *sumsq += d__1 * d__1;
                }
            }
            temp1 = (d__1 = d_lmp_imag(&x[ix]), abs(d__1));
            if (temp1 > 0. || disnan_(&temp1)) {
                if (*scale < temp1) {
                    d__1 = *scale / temp1;
                    *sumsq = *sumsq * (d__1 * d__1) + 1;
                    *scale = temp1;
                } else {
                    d__1 = temp1 / *scale;
                    *sumsq += d__1 * d__1;
                }
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
