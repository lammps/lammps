#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int zgerc_(integer *m, integer *n, doublecomplex *alpha, doublecomplex *x, integer *incx,
           doublecomplex *y, integer *incy, doublecomplex *a, integer *lda)
{
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2;
    void d_lmp_cnjg(doublecomplex *, doublecomplex *);
    integer i__, j, ix, jy, kx, info;
    doublecomplex temp;
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
        xerbla_((char *)"ZGERC ", &info, (ftnlen)6);
        return 0;
    }
    if (*m == 0 || *n == 0 || alpha->r == 0. && alpha->i == 0.) {
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
            i__2 = jy;
            if (y[i__2].r != 0. || y[i__2].i != 0.) {
                d_lmp_cnjg(&z__2, &y[jy]);
                z__1.r = alpha->r * z__2.r - alpha->i * z__2.i,
                z__1.i = alpha->r * z__2.i + alpha->i * z__2.r;
                temp.r = z__1.r, temp.i = z__1.i;
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    i__3 = i__ + j * a_dim1;
                    i__4 = i__ + j * a_dim1;
                    i__5 = i__;
                    z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i,
                    z__2.i = x[i__5].r * temp.i + x[i__5].i * temp.r;
                    z__1.r = a[i__4].r + z__2.r, z__1.i = a[i__4].i + z__2.i;
                    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
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
            i__2 = jy;
            if (y[i__2].r != 0. || y[i__2].i != 0.) {
                d_lmp_cnjg(&z__2, &y[jy]);
                z__1.r = alpha->r * z__2.r - alpha->i * z__2.i,
                z__1.i = alpha->r * z__2.i + alpha->i * z__2.r;
                temp.r = z__1.r, temp.i = z__1.i;
                ix = kx;
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    i__3 = i__ + j * a_dim1;
                    i__4 = i__ + j * a_dim1;
                    i__5 = ix;
                    z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i,
                    z__2.i = x[i__5].r * temp.i + x[i__5].i * temp.r;
                    z__1.r = a[i__4].r + z__2.r, z__1.i = a[i__4].i + z__2.i;
                    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
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
