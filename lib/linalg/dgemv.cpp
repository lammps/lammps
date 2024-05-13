#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dgemv_(char *trans, integer *m, integer *n, doublereal *alpha, doublereal *a, integer *lda,
           doublereal *x, integer *incx, doublereal *beta, doublereal *y, integer *incy,
           ftnlen trans_len)
{
    integer a_dim1, a_offset, i__1, i__2;
    integer i__, j, ix, iy, jx, jy, kx, ky, info;
    doublereal temp;
    integer lenx, leny;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int xerbla_(char *, integer *, ftnlen);
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
        xerbla_((char *)"DGEMV ", &info, (ftnlen)6);
        return 0;
    }
    if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
        return 0;
    }
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
    if (*beta != 1.) {
        if (*incy == 1) {
            if (*beta == 0.) {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    y[i__] = 0.;
                }
            } else {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    y[i__] = *beta * y[i__];
                }
            }
        } else {
            iy = ky;
            if (*beta == 0.) {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    y[iy] = 0.;
                    iy += *incy;
                }
            } else {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    y[iy] = *beta * y[iy];
                    iy += *incy;
                }
            }
        }
    }
    if (*alpha == 0.) {
        return 0;
    }
    if (lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1)) {
        jx = kx;
        if (*incy == 1) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                temp = *alpha * x[jx];
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    y[i__] += temp * a[i__ + j * a_dim1];
                }
                jx += *incx;
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                temp = *alpha * x[jx];
                iy = ky;
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    y[iy] += temp * a[i__ + j * a_dim1];
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
                temp = 0.;
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp += a[i__ + j * a_dim1] * x[i__];
                }
                y[jy] += *alpha * temp;
                jy += *incy;
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                temp = 0.;
                ix = kx;
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp += a[i__ + j * a_dim1] * x[ix];
                    ix += *incx;
                }
                y[jy] += *alpha * temp;
                jy += *incy;
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
