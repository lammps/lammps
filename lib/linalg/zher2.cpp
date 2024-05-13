#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int zher2_(char *uplo, integer *n, doublecomplex *alpha, doublecomplex *x, integer *incx,
           doublecomplex *y, integer *incy, doublecomplex *a, integer *lda, ftnlen uplo_len)
{
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4;
    void d_lmp_cnjg(doublecomplex *, doublecomplex *);
    integer i__, j, ix, iy, jx, jy, kx, ky, info;
    doublecomplex temp1, temp2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int xerbla_(char *, integer *, ftnlen);
    --x;
    --y;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    info = 0;
    if (!lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1) && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        info = 1;
    } else if (*n < 0) {
        info = 2;
    } else if (*incx == 0) {
        info = 5;
    } else if (*incy == 0) {
        info = 7;
    } else if (*lda < max(1, *n)) {
        info = 9;
    }
    if (info != 0) {
        xerbla_((char *)"ZHER2 ", &info, (ftnlen)6);
        return 0;
    }
    if (*n == 0 || alpha->r == 0. && alpha->i == 0.) {
        return 0;
    }
    if (*incx != 1 || *incy != 1) {
        if (*incx > 0) {
            kx = 1;
        } else {
            kx = 1 - (*n - 1) * *incx;
        }
        if (*incy > 0) {
            ky = 1;
        } else {
            ky = 1 - (*n - 1) * *incy;
        }
        jx = kx;
        jy = ky;
    }
    if (lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1)) {
        if (*incx == 1 && *incy == 1) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = j;
                i__3 = j;
                if (x[i__2].r != 0. || x[i__2].i != 0. || (y[i__3].r != 0. || y[i__3].i != 0.)) {
                    d_lmp_cnjg(&z__2, &y[j]);
                    z__1.r = alpha->r * z__2.r - alpha->i * z__2.i,
                    z__1.i = alpha->r * z__2.i + alpha->i * z__2.r;
                    temp1.r = z__1.r, temp1.i = z__1.i;
                    i__2 = j;
                    z__2.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i,
                    z__2.i = alpha->r * x[i__2].i + alpha->i * x[i__2].r;
                    d_lmp_cnjg(&z__1, &z__2);
                    temp2.r = z__1.r, temp2.i = z__1.i;
                    i__2 = j - 1;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * a_dim1;
                        i__4 = i__ + j * a_dim1;
                        i__5 = i__;
                        z__3.r = x[i__5].r * temp1.r - x[i__5].i * temp1.i,
                        z__3.i = x[i__5].r * temp1.i + x[i__5].i * temp1.r;
                        z__2.r = a[i__4].r + z__3.r, z__2.i = a[i__4].i + z__3.i;
                        i__6 = i__;
                        z__4.r = y[i__6].r * temp2.r - y[i__6].i * temp2.i,
                        z__4.i = y[i__6].r * temp2.i + y[i__6].i * temp2.r;
                        z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
                        a[i__3].r = z__1.r, a[i__3].i = z__1.i;
                    }
                    i__2 = j + j * a_dim1;
                    i__3 = j + j * a_dim1;
                    i__4 = j;
                    z__2.r = x[i__4].r * temp1.r - x[i__4].i * temp1.i,
                    z__2.i = x[i__4].r * temp1.i + x[i__4].i * temp1.r;
                    i__5 = j;
                    z__3.r = y[i__5].r * temp2.r - y[i__5].i * temp2.i,
                    z__3.i = y[i__5].r * temp2.i + y[i__5].i * temp2.r;
                    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                    d__1 = a[i__3].r + z__1.r;
                    a[i__2].r = d__1, a[i__2].i = 0.;
                } else {
                    i__2 = j + j * a_dim1;
                    i__3 = j + j * a_dim1;
                    d__1 = a[i__3].r;
                    a[i__2].r = d__1, a[i__2].i = 0.;
                }
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = jx;
                i__3 = jy;
                if (x[i__2].r != 0. || x[i__2].i != 0. || (y[i__3].r != 0. || y[i__3].i != 0.)) {
                    d_lmp_cnjg(&z__2, &y[jy]);
                    z__1.r = alpha->r * z__2.r - alpha->i * z__2.i,
                    z__1.i = alpha->r * z__2.i + alpha->i * z__2.r;
                    temp1.r = z__1.r, temp1.i = z__1.i;
                    i__2 = jx;
                    z__2.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i,
                    z__2.i = alpha->r * x[i__2].i + alpha->i * x[i__2].r;
                    d_lmp_cnjg(&z__1, &z__2);
                    temp2.r = z__1.r, temp2.i = z__1.i;
                    ix = kx;
                    iy = ky;
                    i__2 = j - 1;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * a_dim1;
                        i__4 = i__ + j * a_dim1;
                        i__5 = ix;
                        z__3.r = x[i__5].r * temp1.r - x[i__5].i * temp1.i,
                        z__3.i = x[i__5].r * temp1.i + x[i__5].i * temp1.r;
                        z__2.r = a[i__4].r + z__3.r, z__2.i = a[i__4].i + z__3.i;
                        i__6 = iy;
                        z__4.r = y[i__6].r * temp2.r - y[i__6].i * temp2.i,
                        z__4.i = y[i__6].r * temp2.i + y[i__6].i * temp2.r;
                        z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
                        a[i__3].r = z__1.r, a[i__3].i = z__1.i;
                        ix += *incx;
                        iy += *incy;
                    }
                    i__2 = j + j * a_dim1;
                    i__3 = j + j * a_dim1;
                    i__4 = jx;
                    z__2.r = x[i__4].r * temp1.r - x[i__4].i * temp1.i,
                    z__2.i = x[i__4].r * temp1.i + x[i__4].i * temp1.r;
                    i__5 = jy;
                    z__3.r = y[i__5].r * temp2.r - y[i__5].i * temp2.i,
                    z__3.i = y[i__5].r * temp2.i + y[i__5].i * temp2.r;
                    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                    d__1 = a[i__3].r + z__1.r;
                    a[i__2].r = d__1, a[i__2].i = 0.;
                } else {
                    i__2 = j + j * a_dim1;
                    i__3 = j + j * a_dim1;
                    d__1 = a[i__3].r;
                    a[i__2].r = d__1, a[i__2].i = 0.;
                }
                jx += *incx;
                jy += *incy;
            }
        }
    } else {
        if (*incx == 1 && *incy == 1) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = j;
                i__3 = j;
                if (x[i__2].r != 0. || x[i__2].i != 0. || (y[i__3].r != 0. || y[i__3].i != 0.)) {
                    d_lmp_cnjg(&z__2, &y[j]);
                    z__1.r = alpha->r * z__2.r - alpha->i * z__2.i,
                    z__1.i = alpha->r * z__2.i + alpha->i * z__2.r;
                    temp1.r = z__1.r, temp1.i = z__1.i;
                    i__2 = j;
                    z__2.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i,
                    z__2.i = alpha->r * x[i__2].i + alpha->i * x[i__2].r;
                    d_lmp_cnjg(&z__1, &z__2);
                    temp2.r = z__1.r, temp2.i = z__1.i;
                    i__2 = j + j * a_dim1;
                    i__3 = j + j * a_dim1;
                    i__4 = j;
                    z__2.r = x[i__4].r * temp1.r - x[i__4].i * temp1.i,
                    z__2.i = x[i__4].r * temp1.i + x[i__4].i * temp1.r;
                    i__5 = j;
                    z__3.r = y[i__5].r * temp2.r - y[i__5].i * temp2.i,
                    z__3.i = y[i__5].r * temp2.i + y[i__5].i * temp2.r;
                    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                    d__1 = a[i__3].r + z__1.r;
                    a[i__2].r = d__1, a[i__2].i = 0.;
                    i__2 = *n;
                    for (i__ = j + 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * a_dim1;
                        i__4 = i__ + j * a_dim1;
                        i__5 = i__;
                        z__3.r = x[i__5].r * temp1.r - x[i__5].i * temp1.i,
                        z__3.i = x[i__5].r * temp1.i + x[i__5].i * temp1.r;
                        z__2.r = a[i__4].r + z__3.r, z__2.i = a[i__4].i + z__3.i;
                        i__6 = i__;
                        z__4.r = y[i__6].r * temp2.r - y[i__6].i * temp2.i,
                        z__4.i = y[i__6].r * temp2.i + y[i__6].i * temp2.r;
                        z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
                        a[i__3].r = z__1.r, a[i__3].i = z__1.i;
                    }
                } else {
                    i__2 = j + j * a_dim1;
                    i__3 = j + j * a_dim1;
                    d__1 = a[i__3].r;
                    a[i__2].r = d__1, a[i__2].i = 0.;
                }
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = jx;
                i__3 = jy;
                if (x[i__2].r != 0. || x[i__2].i != 0. || (y[i__3].r != 0. || y[i__3].i != 0.)) {
                    d_lmp_cnjg(&z__2, &y[jy]);
                    z__1.r = alpha->r * z__2.r - alpha->i * z__2.i,
                    z__1.i = alpha->r * z__2.i + alpha->i * z__2.r;
                    temp1.r = z__1.r, temp1.i = z__1.i;
                    i__2 = jx;
                    z__2.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i,
                    z__2.i = alpha->r * x[i__2].i + alpha->i * x[i__2].r;
                    d_lmp_cnjg(&z__1, &z__2);
                    temp2.r = z__1.r, temp2.i = z__1.i;
                    i__2 = j + j * a_dim1;
                    i__3 = j + j * a_dim1;
                    i__4 = jx;
                    z__2.r = x[i__4].r * temp1.r - x[i__4].i * temp1.i,
                    z__2.i = x[i__4].r * temp1.i + x[i__4].i * temp1.r;
                    i__5 = jy;
                    z__3.r = y[i__5].r * temp2.r - y[i__5].i * temp2.i,
                    z__3.i = y[i__5].r * temp2.i + y[i__5].i * temp2.r;
                    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                    d__1 = a[i__3].r + z__1.r;
                    a[i__2].r = d__1, a[i__2].i = 0.;
                    ix = jx;
                    iy = jy;
                    i__2 = *n;
                    for (i__ = j + 1; i__ <= i__2; ++i__) {
                        ix += *incx;
                        iy += *incy;
                        i__3 = i__ + j * a_dim1;
                        i__4 = i__ + j * a_dim1;
                        i__5 = ix;
                        z__3.r = x[i__5].r * temp1.r - x[i__5].i * temp1.i,
                        z__3.i = x[i__5].r * temp1.i + x[i__5].i * temp1.r;
                        z__2.r = a[i__4].r + z__3.r, z__2.i = a[i__4].i + z__3.i;
                        i__6 = iy;
                        z__4.r = y[i__6].r * temp2.r - y[i__6].i * temp2.i,
                        z__4.i = y[i__6].r * temp2.i + y[i__6].i * temp2.r;
                        z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
                        a[i__3].r = z__1.r, a[i__3].i = z__1.i;
                    }
                } else {
                    i__2 = j + j * a_dim1;
                    i__3 = j + j * a_dim1;
                    d__1 = a[i__3].r;
                    a[i__2].r = d__1, a[i__2].i = 0.;
                }
                jx += *incx;
                jy += *incy;
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
