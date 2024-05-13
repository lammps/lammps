#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
doublereal dasum_(integer *n, doublereal *dx, integer *incx)
{
    integer i__1, i__2;
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5, d__6;
    integer i__, m, mp1;
    doublereal dtemp;
    integer nincx;
    --dx;
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0 || *incx <= 0) {
        return ret_val;
    }
    if (*incx == 1) {
        m = *n % 6;
        if (m != 0) {
            i__1 = m;
            for (i__ = 1; i__ <= i__1; ++i__) {
                dtemp += (d__1 = dx[i__], abs(d__1));
            }
            if (*n < 6) {
                ret_val = dtemp;
                return ret_val;
            }
        }
        mp1 = m + 1;
        i__1 = *n;
        for (i__ = mp1; i__ <= i__1; i__ += 6) {
            dtemp = dtemp + (d__1 = dx[i__], abs(d__1)) + (d__2 = dx[i__ + 1], abs(d__2)) +
                    (d__3 = dx[i__ + 2], abs(d__3)) + (d__4 = dx[i__ + 3], abs(d__4)) +
                    (d__5 = dx[i__ + 4], abs(d__5)) + (d__6 = dx[i__ + 5], abs(d__6));
        }
    } else {
        nincx = *n * *incx;
        i__1 = nincx;
        i__2 = *incx;
        for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
            dtemp += (d__1 = dx[i__], abs(d__1));
        }
    }
    ret_val = dtemp;
    return ret_val;
}
#ifdef __cplusplus
}
#endif
