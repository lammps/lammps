#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dscal_(integer *n, doublereal *da, doublereal *dx, integer *incx)
{
    integer i__1, i__2;
    integer i__, m, mp1, nincx;
    --dx;
    if (*n <= 0 || *incx <= 0 || *da == 1.) {
        return 0;
    }
    if (*incx == 1) {
        m = *n % 5;
        if (m != 0) {
            i__1 = m;
            for (i__ = 1; i__ <= i__1; ++i__) {
                dx[i__] = *da * dx[i__];
            }
            if (*n < 5) {
                return 0;
            }
        }
        mp1 = m + 1;
        i__1 = *n;
        for (i__ = mp1; i__ <= i__1; i__ += 5) {
            dx[i__] = *da * dx[i__];
            dx[i__ + 1] = *da * dx[i__ + 1];
            dx[i__ + 2] = *da * dx[i__ + 2];
            dx[i__ + 3] = *da * dx[i__ + 3];
            dx[i__ + 4] = *da * dx[i__ + 4];
        }
    } else {
        nincx = *n * *incx;
        i__1 = nincx;
        i__2 = *incx;
        for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
            dx[i__] = *da * dx[i__];
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
