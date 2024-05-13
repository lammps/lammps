#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlassq_(integer *n, doublereal *x, integer *incx, doublereal *scale, doublereal *sumsq)
{
    integer i__1, i__2;
    doublereal d__1;
    integer ix;
    doublereal absxi;
    extern logical disnan_(doublereal *);
    --x;
    if (*n > 0) {
        i__1 = (*n - 1) * *incx + 1;
        i__2 = *incx;
        for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
            absxi = (d__1 = x[ix], abs(d__1));
            if (absxi > 0. || disnan_(&absxi)) {
                if (*scale < absxi) {
                    d__1 = *scale / absxi;
                    *sumsq = *sumsq * (d__1 * d__1) + 1;
                    *scale = absxi;
                } else {
                    d__1 = absxi / *scale;
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
