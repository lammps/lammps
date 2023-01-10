#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlaswp_(integer *n, doublereal *a, integer *lda, integer *k1, integer *k2, integer *ipiv,
            integer *incx)
{
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    integer i__, j, k, i1, i2, n32, ip, ix, ix0, inc;
    doublereal temp;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    if (*incx > 0) {
        ix0 = *k1;
        i1 = *k1;
        i2 = *k2;
        inc = 1;
    } else if (*incx < 0) {
        ix0 = *k1 + (*k1 - *k2) * *incx;
        i1 = *k2;
        i2 = *k1;
        inc = -1;
    } else {
        return 0;
    }
    n32 = *n / 32 << 5;
    if (n32 != 0) {
        i__1 = n32;
        for (j = 1; j <= i__1; j += 32) {
            ix = ix0;
            i__2 = i2;
            i__3 = inc;
            for (i__ = i1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) {
                ip = ipiv[ix];
                if (ip != i__) {
                    i__4 = j + 31;
                    for (k = j; k <= i__4; ++k) {
                        temp = a[i__ + k * a_dim1];
                        a[i__ + k * a_dim1] = a[ip + k * a_dim1];
                        a[ip + k * a_dim1] = temp;
                    }
                }
                ix += *incx;
            }
        }
    }
    if (n32 != *n) {
        ++n32;
        ix = ix0;
        i__1 = i2;
        i__3 = inc;
        for (i__ = i1; i__3 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__3) {
            ip = ipiv[ix];
            if (ip != i__) {
                i__2 = *n;
                for (k = n32; k <= i__2; ++k) {
                    temp = a[i__ + k * a_dim1];
                    a[i__ + k * a_dim1] = a[ip + k * a_dim1];
                    a[ip + k * a_dim1] = temp;
                }
            }
            ix += *incx;
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
