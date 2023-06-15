#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int zhpr_(char *uplo, integer *n, doublereal *alpha, doublecomplex *x, integer *incx,
          doublecomplex *ap, ftnlen uplo_len)
{
    integer i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;
    doublecomplex z__1, z__2;
    void d_lmp_cnjg(doublecomplex *, doublecomplex *);
    integer i__, j, k, kk, ix, jx, kx, info;
    doublecomplex temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int xerbla_(char *, integer *, ftnlen);
    --ap;
    --x;
    info = 0;
    if (!lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1) && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        info = 1;
    } else if (*n < 0) {
        info = 2;
    } else if (*incx == 0) {
        info = 5;
    }
    if (info != 0) {
        xerbla_((char *)"ZHPR  ", &info, (ftnlen)6);
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
    kk = 1;
    if (lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1)) {
        if (*incx == 1) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = j;
                if (x[i__2].r != 0. || x[i__2].i != 0.) {
                    d_lmp_cnjg(&z__2, &x[j]);
                    z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
                    temp.r = z__1.r, temp.i = z__1.i;
                    k = kk;
                    i__2 = j - 1;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = k;
                        i__4 = k;
                        i__5 = i__;
                        z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i,
                        z__2.i = x[i__5].r * temp.i + x[i__5].i * temp.r;
                        z__1.r = ap[i__4].r + z__2.r, z__1.i = ap[i__4].i + z__2.i;
                        ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
                        ++k;
                    }
                    i__2 = kk + j - 1;
                    i__3 = kk + j - 1;
                    i__4 = j;
                    z__1.r = x[i__4].r * temp.r - x[i__4].i * temp.i,
                    z__1.i = x[i__4].r * temp.i + x[i__4].i * temp.r;
                    d__1 = ap[i__3].r + z__1.r;
                    ap[i__2].r = d__1, ap[i__2].i = 0.;
                } else {
                    i__2 = kk + j - 1;
                    i__3 = kk + j - 1;
                    d__1 = ap[i__3].r;
                    ap[i__2].r = d__1, ap[i__2].i = 0.;
                }
                kk += j;
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
                    i__2 = kk + j - 2;
                    for (k = kk; k <= i__2; ++k) {
                        i__3 = k;
                        i__4 = k;
                        i__5 = ix;
                        z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i,
                        z__2.i = x[i__5].r * temp.i + x[i__5].i * temp.r;
                        z__1.r = ap[i__4].r + z__2.r, z__1.i = ap[i__4].i + z__2.i;
                        ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
                        ix += *incx;
                    }
                    i__2 = kk + j - 1;
                    i__3 = kk + j - 1;
                    i__4 = jx;
                    z__1.r = x[i__4].r * temp.r - x[i__4].i * temp.i,
                    z__1.i = x[i__4].r * temp.i + x[i__4].i * temp.r;
                    d__1 = ap[i__3].r + z__1.r;
                    ap[i__2].r = d__1, ap[i__2].i = 0.;
                } else {
                    i__2 = kk + j - 1;
                    i__3 = kk + j - 1;
                    d__1 = ap[i__3].r;
                    ap[i__2].r = d__1, ap[i__2].i = 0.;
                }
                jx += *incx;
                kk += j;
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
                    i__2 = kk;
                    i__3 = kk;
                    i__4 = j;
                    z__1.r = temp.r * x[i__4].r - temp.i * x[i__4].i,
                    z__1.i = temp.r * x[i__4].i + temp.i * x[i__4].r;
                    d__1 = ap[i__3].r + z__1.r;
                    ap[i__2].r = d__1, ap[i__2].i = 0.;
                    k = kk + 1;
                    i__2 = *n;
                    for (i__ = j + 1; i__ <= i__2; ++i__) {
                        i__3 = k;
                        i__4 = k;
                        i__5 = i__;
                        z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i,
                        z__2.i = x[i__5].r * temp.i + x[i__5].i * temp.r;
                        z__1.r = ap[i__4].r + z__2.r, z__1.i = ap[i__4].i + z__2.i;
                        ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
                        ++k;
                    }
                } else {
                    i__2 = kk;
                    i__3 = kk;
                    d__1 = ap[i__3].r;
                    ap[i__2].r = d__1, ap[i__2].i = 0.;
                }
                kk = kk + *n - j + 1;
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
                    i__2 = kk;
                    i__3 = kk;
                    i__4 = jx;
                    z__1.r = temp.r * x[i__4].r - temp.i * x[i__4].i,
                    z__1.i = temp.r * x[i__4].i + temp.i * x[i__4].r;
                    d__1 = ap[i__3].r + z__1.r;
                    ap[i__2].r = d__1, ap[i__2].i = 0.;
                    ix = jx;
                    i__2 = kk + *n - j;
                    for (k = kk + 1; k <= i__2; ++k) {
                        ix += *incx;
                        i__3 = k;
                        i__4 = k;
                        i__5 = ix;
                        z__2.r = x[i__5].r * temp.r - x[i__5].i * temp.i,
                        z__2.i = x[i__5].r * temp.i + x[i__5].i * temp.r;
                        z__1.r = ap[i__4].r + z__2.r, z__1.i = ap[i__4].i + z__2.i;
                        ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
                    }
                } else {
                    i__2 = kk;
                    i__3 = kk;
                    d__1 = ap[i__3].r;
                    ap[i__2].r = d__1, ap[i__2].i = 0.;
                }
                jx += *incx;
                kk = kk + *n - j + 1;
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
