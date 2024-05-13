#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dsyr2_(char *uplo, integer *n, doublereal *alpha, doublereal *x, integer *incx, doublereal *y,
           integer *incy, doublereal *a, integer *lda, ftnlen uplo_len)
{
    integer a_dim1, a_offset, i__1, i__2;
    integer i__, j, ix, iy, jx, jy, kx, ky, info;
    doublereal temp1, temp2;
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
        xerbla_((char *)"DSYR2 ", &info, (ftnlen)6);
        return 0;
    }
    if (*n == 0 || *alpha == 0.) {
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
                if (x[j] != 0. || y[j] != 0.) {
                    temp1 = *alpha * y[j];
                    temp2 = *alpha * x[j];
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[i__] * temp1 + y[i__] * temp2;
                    }
                }
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (x[jx] != 0. || y[jy] != 0.) {
                    temp1 = *alpha * y[jy];
                    temp2 = *alpha * x[jx];
                    ix = kx;
                    iy = ky;
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[ix] * temp1 + y[iy] * temp2;
                        ix += *incx;
                        iy += *incy;
                    }
                }
                jx += *incx;
                jy += *incy;
            }
        }
    } else {
        if (*incx == 1 && *incy == 1) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (x[j] != 0. || y[j] != 0.) {
                    temp1 = *alpha * y[j];
                    temp2 = *alpha * x[j];
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__) {
                        a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[i__] * temp1 + y[i__] * temp2;
                    }
                }
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (x[jx] != 0. || y[jy] != 0.) {
                    temp1 = *alpha * y[jy];
                    temp2 = *alpha * x[jx];
                    ix = jx;
                    iy = jy;
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__) {
                        a[i__ + j * a_dim1] = a[i__ + j * a_dim1] + x[ix] * temp1 + y[iy] * temp2;
                        ix += *incx;
                        iy += *incy;
                    }
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
