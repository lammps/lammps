#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int zgemv_(char *trans, integer *m, integer *n, doublecomplex *alpha, doublecomplex *a,
           integer *lda, doublecomplex *x, integer *incx, doublecomplex *beta, doublecomplex *y,
           integer *incy, ftnlen trans_len)
{
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3;
    void d_lmp_cnjg(doublecomplex *, doublecomplex *);
    integer i__, j, ix, iy, jx, jy, kx, ky, info;
    doublecomplex temp;
    integer lenx, leny;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int xerbla_(char *, integer *, ftnlen);
    logical noconj;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    --y;
    info = 0;
    if (!lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1) && !lsame_(trans, (char *)"T", (ftnlen)1, (ftnlen)1) &&
        !lsame_(trans, (char *)"C", (ftnlen)1, (ftnlen)1)) {
        info = 1;
    } else if (*m < 0) {
        info = 2;
    } else if (*n < 0) {
        info = 3;
    } else if (*lda < max(1, *m)) {
        info = 6;
    } else if (*incx == 0) {
        info = 8;
    } else if (*incy == 0) {
        info = 11;
    }
    if (info != 0) {
        xerbla_((char *)"ZGEMV ", &info, (ftnlen)6);
        return 0;
    }
    if (*m == 0 || *n == 0 ||
        alpha->r == 0. && alpha->i == 0. && (beta->r == 1. && beta->i == 0.)) {
        return 0;
    }
    noconj = lsame_(trans, (char *)"T", (ftnlen)1, (ftnlen)1);
    if (lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1)) {
        lenx = *n;
        leny = *m;
    } else {
        lenx = *m;
        leny = *n;
    }
    if (*incx > 0) {
        kx = 1;
    } else {
        kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
        ky = 1;
    } else {
        ky = 1 - (leny - 1) * *incy;
    }
    if (beta->r != 1. || beta->i != 0.) {
        if (*incy == 1) {
            if (beta->r == 0. && beta->i == 0.) {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    i__2 = i__;
                    y[i__2].r = 0., y[i__2].i = 0.;
                }
            } else {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    i__2 = i__;
                    i__3 = i__;
                    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i,
                    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3].r;
                    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
                }
            }
        } else {
            iy = ky;
            if (beta->r == 0. && beta->i == 0.) {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    i__2 = iy;
                    y[i__2].r = 0., y[i__2].i = 0.;
                    iy += *incy;
                }
            } else {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    i__2 = iy;
                    i__3 = iy;
                    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i,
                    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3].r;
                    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
                    iy += *incy;
                }
            }
        }
    }
    if (alpha->r == 0. && alpha->i == 0.) {
        return 0;
    }
    if (lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1)) {
        jx = kx;
        if (*incy == 1) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = jx;
                z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i,
                z__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2].r;
                temp.r = z__1.r, temp.i = z__1.i;
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    i__3 = i__;
                    i__4 = i__;
                    i__5 = i__ + j * a_dim1;
                    z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i,
                    z__2.i = temp.r * a[i__5].i + temp.i * a[i__5].r;
                    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
                    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
                }
                jx += *incx;
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = jx;
                z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i,
                z__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2].r;
                temp.r = z__1.r, temp.i = z__1.i;
                iy = ky;
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    i__3 = iy;
                    i__4 = iy;
                    i__5 = i__ + j * a_dim1;
                    z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i,
                    z__2.i = temp.r * a[i__5].i + temp.i * a[i__5].r;
                    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
                    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
                    iy += *incy;
                }
                jx += *incx;
            }
        }
    } else {
        jy = ky;
        if (*incx == 1) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                temp.r = 0., temp.i = 0.;
                if (noconj) {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * a_dim1;
                        i__4 = i__;
                        z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i,
                        z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4].r;
                        z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                        temp.r = z__1.r, temp.i = z__1.i;
                    }
                } else {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        d_lmp_cnjg(&z__3, &a[i__ + j * a_dim1]);
                        i__3 = i__;
                        z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i,
                        z__2.i = z__3.r * x[i__3].i + z__3.i * x[i__3].r;
                        z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                        temp.r = z__1.r, temp.i = z__1.i;
                    }
                }
                i__2 = jy;
                i__3 = jy;
                z__2.r = alpha->r * temp.r - alpha->i * temp.i,
                z__2.i = alpha->r * temp.i + alpha->i * temp.r;
                z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
                y[i__2].r = z__1.r, y[i__2].i = z__1.i;
                jy += *incy;
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                temp.r = 0., temp.i = 0.;
                ix = kx;
                if (noconj) {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * a_dim1;
                        i__4 = ix;
                        z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i,
                        z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[i__4].r;
                        z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                        temp.r = z__1.r, temp.i = z__1.i;
                        ix += *incx;
                    }
                } else {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        d_lmp_cnjg(&z__3, &a[i__ + j * a_dim1]);
                        i__3 = ix;
                        z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i,
                        z__2.i = z__3.r * x[i__3].i + z__3.i * x[i__3].r;
                        z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                        temp.r = z__1.r, temp.i = z__1.i;
                        ix += *incx;
                    }
                }
                i__2 = jy;
                i__3 = jy;
                z__2.r = alpha->r * temp.r - alpha->i * temp.i,
                z__2.i = alpha->r * temp.i + alpha->i * temp.r;
                z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
                y[i__2].r = z__1.r, y[i__2].i = z__1.i;
                jy += *incy;
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
