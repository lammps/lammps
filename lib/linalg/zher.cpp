#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int zher_(char *uplo, integer *n, doublereal *alpha, doublecomplex *x, integer *incx,
          doublecomplex *a, integer *lda, ftnlen uplo_len)
{
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;
    doublecomplex z__1, z__2;
    void d_lmp_cnjg(doublecomplex *, doublecomplex *);
    integer i__, j, ix, jx, kx, info;
    doublecomplex temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int xerbla_(char *, integer *, ftnlen);
    --x;
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
    } else if (*lda < max(1, *n)) {
        info = 7;
    }
    if (info != 0) {
        xerbla_((char *)"ZHER  ", &info, (ftnlen)6);
        return 0;
    }
    if (*n == 0 || *alpha == 0.) {
        return 0;
    }
    if (*incx <= 0) {
        kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
        kx = 1;
    }
    if (lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1)) {
        if (*incx == 1) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = j;
                if (x[i__2].r != 0. || x[i__2].i != 0.) {
                    d_lmp_cnjg(&z__2, &x[j]);
                    z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
                    temp.r = z__1.r, temp.i = z__1.i;
                    i__2 = j - 1;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * a_dim1;
                        i__4 = i__ + j * a_dim1;
                        i__5 = i__;
                        z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i,
                        z__2.i = x[i__5].r * temp.i + x[i__5].i * temp.r;
                        z__1.r = a[i__4].r + z__2.r, z__1.i = a[i__4].i + z__2.i;
                        a[i__3].r = z__1.r, a[i__3].i = z__1.i;
                    }
                    i__2 = j + j * a_dim1;
                    i__3 = j + j * a_dim1;
                    i__4 = j;
                    z__1.r = x[i__4].r * temp.r - x[i__4].i * temp.i,
                    z__1.i = x[i__4].r * temp.i + x[i__4].i * temp.r;
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
            jx = kx;
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = jx;
                if (x[i__2].r != 0. || x[i__2].i != 0.) {
                    d_lmp_cnjg(&z__2, &x[jx]);
                    z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
                    temp.r = z__1.r, temp.i = z__1.i;
                    ix = kx;
                    i__2 = j - 1;
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
                    i__2 = j + j * a_dim1;
                    i__3 = j + j * a_dim1;
                    i__4 = jx;
                    z__1.r = x[i__4].r * temp.r - x[i__4].i * temp.i,
                    z__1.i = x[i__4].r * temp.i + x[i__4].i * temp.r;
                    d__1 = a[i__3].r + z__1.r;
                    a[i__2].r = d__1, a[i__2].i = 0.;
                } else {
                    i__2 = j + j * a_dim1;
                    i__3 = j + j * a_dim1;
                    d__1 = a[i__3].r;
                    a[i__2].r = d__1, a[i__2].i = 0.;
                }
                jx += *incx;
            }
        }
    } else {
        if (*incx == 1) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = j;
                if (x[i__2].r != 0. || x[i__2].i != 0.) {
                    d_lmp_cnjg(&z__2, &x[j]);
                    z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
                    temp.r = z__1.r, temp.i = z__1.i;
                    i__2 = j + j * a_dim1;
                    i__3 = j + j * a_dim1;
                    i__4 = j;
                    z__1.r = temp.r * x[i__4].r - temp.i * x[i__4].i,
                    z__1.i = temp.r * x[i__4].i + temp.i * x[i__4].r;
                    d__1 = a[i__3].r + z__1.r;
                    a[i__2].r = d__1, a[i__2].i = 0.;
                    i__2 = *n;
                    for (i__ = j + 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * a_dim1;
                        i__4 = i__ + j * a_dim1;
                        i__5 = i__;
                        z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i,
                        z__2.i = x[i__5].r * temp.i + x[i__5].i * temp.r;
                        z__1.r = a[i__4].r + z__2.r, z__1.i = a[i__4].i + z__2.i;
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
            jx = kx;
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = jx;
                if (x[i__2].r != 0. || x[i__2].i != 0.) {
                    d_lmp_cnjg(&z__2, &x[jx]);
                    z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
                    temp.r = z__1.r, temp.i = z__1.i;
                    i__2 = j + j * a_dim1;
                    i__3 = j + j * a_dim1;
                    i__4 = jx;
                    z__1.r = temp.r * x[i__4].r - temp.i * x[i__4].i,
                    z__1.i = temp.r * x[i__4].i + temp.i * x[i__4].r;
                    d__1 = a[i__3].r + z__1.r;
                    a[i__2].r = d__1, a[i__2].i = 0.;
                    ix = jx;
                    i__2 = *n;
                    for (i__ = j + 1; i__ <= i__2; ++i__) {
                        ix += *incx;
                        i__3 = i__ + j * a_dim1;
                        i__4 = i__ + j * a_dim1;
                        i__5 = ix;
                        z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i,
                        z__2.i = x[i__5].r * temp.i + x[i__5].i * temp.r;
                        z__1.r = a[i__4].r + z__2.r, z__1.i = a[i__4].i + z__2.i;
                        a[i__3].r = z__1.r, a[i__3].i = z__1.i;
                    }
                } else {
                    i__2 = j + j * a_dim1;
                    i__3 = j + j * a_dim1;
                    d__1 = a[i__3].r;
                    a[i__2].r = d__1, a[i__2].i = 0.;
                }
                jx += *incx;
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
