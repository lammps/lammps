#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dger_(integer *m, integer *n, doublereal *alpha, doublereal *x, integer *incx, doublereal *y,
          integer *incy, doublereal *a, integer *lda)
{
    integer a_dim1, a_offset, i__1, i__2;
    integer i__, j, ix, jy, kx, info;
    doublereal temp;
    extern int xerbla_(char *, integer *, ftnlen);
    --x;
    --y;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    info = 0;
    if (*m < 0) {
        info = 1;
    } else if (*n < 0) {
        info = 2;
    } else if (*incx == 0) {
        info = 5;
    } else if (*incy == 0) {
        info = 7;
    } else if (*lda < max(1, *m)) {
        info = 9;
    }
    if (info != 0) {
        xerbla_((char *)"DGER  ", &info, (ftnlen)6);
        return 0;
    }
    if (*m == 0 || *n == 0 || *alpha == 0.) {
        return 0;
    }
    if (*incy > 0) {
        jy = 1;
    } else {
        jy = 1 - (*n - 1) * *incy;
    }
    if (*incx == 1) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            if (y[jy] != 0.) {
                temp = *alpha * y[jy];
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    a[i__ + j * a_dim1] += x[i__] * temp;
                }
            }
            jy += *incy;
        }
    } else {
        if (*incx > 0) {
            kx = 1;
        } else {
            kx = 1 - (*m - 1) * *incx;
        }
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            if (y[jy] != 0.) {
                temp = *alpha * y[jy];
                ix = kx;
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    a[i__ + j * a_dim1] += x[ix] * temp;
                    ix += *incx;
                }
            }
            jy += *incy;
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
