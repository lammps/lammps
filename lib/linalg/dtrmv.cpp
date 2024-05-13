#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dtrmv_(char *uplo, char *trans, char *diag, integer *n, doublereal *a, integer *lda,
           doublereal *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len)
{
    integer a_dim1, a_offset, i__1, i__2;
    integer i__, j, ix, jx, kx, info;
    doublereal temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int xerbla_(char *, integer *, ftnlen);
    logical nounit;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    info = 0;
    if (!lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1) && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        info = 1;
    } else if (!lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1) &&
               !lsame_(trans, (char *)"T", (ftnlen)1, (ftnlen)1) &&
               !lsame_(trans, (char *)"C", (ftnlen)1, (ftnlen)1)) {
        info = 2;
    } else if (!lsame_(diag, (char *)"U", (ftnlen)1, (ftnlen)1) &&
               !lsame_(diag, (char *)"N", (ftnlen)1, (ftnlen)1)) {
        info = 3;
    } else if (*n < 0) {
        info = 4;
    } else if (*lda < max(1, *n)) {
        info = 6;
    } else if (*incx == 0) {
        info = 8;
    }
    if (info != 0) {
        xerbla_((char *)"DTRMV ", &info, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    nounit = lsame_(diag, (char *)"N", (ftnlen)1, (ftnlen)1);
    if (*incx <= 0) {
        kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
        kx = 1;
    }
    if (lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1)) {
        if (lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1)) {
            if (*incx == 1) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    if (x[j] != 0.) {
                        temp = x[j];
                        i__2 = j - 1;
                        for (i__ = 1; i__ <= i__2; ++i__) {
                            x[i__] += temp * a[i__ + j * a_dim1];
                        }
                        if (nounit) {
                            x[j] *= a[j + j * a_dim1];
                        }
                    }
                }
            } else {
                jx = kx;
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    if (x[jx] != 0.) {
                        temp = x[jx];
                        ix = kx;
                        i__2 = j - 1;
                        for (i__ = 1; i__ <= i__2; ++i__) {
                            x[ix] += temp * a[i__ + j * a_dim1];
                            ix += *incx;
                        }
                        if (nounit) {
                            x[jx] *= a[j + j * a_dim1];
                        }
                    }
                    jx += *incx;
                }
            }
        } else {
            if (*incx == 1) {
                for (j = *n; j >= 1; --j) {
                    if (x[j] != 0.) {
                        temp = x[j];
                        i__1 = j + 1;
                        for (i__ = *n; i__ >= i__1; --i__) {
                            x[i__] += temp * a[i__ + j * a_dim1];
                        }
                        if (nounit) {
                            x[j] *= a[j + j * a_dim1];
                        }
                    }
                }
            } else {
                kx += (*n - 1) * *incx;
                jx = kx;
                for (j = *n; j >= 1; --j) {
                    if (x[jx] != 0.) {
                        temp = x[jx];
                        ix = kx;
                        i__1 = j + 1;
                        for (i__ = *n; i__ >= i__1; --i__) {
                            x[ix] += temp * a[i__ + j * a_dim1];
                            ix -= *incx;
                        }
                        if (nounit) {
                            x[jx] *= a[j + j * a_dim1];
                        }
                    }
                    jx -= *incx;
                }
            }
        }
    } else {
        if (lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1)) {
            if (*incx == 1) {
                for (j = *n; j >= 1; --j) {
                    temp = x[j];
                    if (nounit) {
                        temp *= a[j + j * a_dim1];
                    }
                    for (i__ = j - 1; i__ >= 1; --i__) {
                        temp += a[i__ + j * a_dim1] * x[i__];
                    }
                    x[j] = temp;
                }
            } else {
                jx = kx + (*n - 1) * *incx;
                for (j = *n; j >= 1; --j) {
                    temp = x[jx];
                    ix = jx;
                    if (nounit) {
                        temp *= a[j + j * a_dim1];
                    }
                    for (i__ = j - 1; i__ >= 1; --i__) {
                        ix -= *incx;
                        temp += a[i__ + j * a_dim1] * x[ix];
                    }
                    x[jx] = temp;
                    jx -= *incx;
                }
            }
        } else {
            if (*incx == 1) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    temp = x[j];
                    if (nounit) {
                        temp *= a[j + j * a_dim1];
                    }
                    i__2 = *n;
                    for (i__ = j + 1; i__ <= i__2; ++i__) {
                        temp += a[i__ + j * a_dim1] * x[i__];
                    }
                    x[j] = temp;
                }
            } else {
                jx = kx;
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    temp = x[jx];
                    ix = jx;
                    if (nounit) {
                        temp *= a[j + j * a_dim1];
                    }
                    i__2 = *n;
                    for (i__ = j + 1; i__ <= i__2; ++i__) {
                        ix += *incx;
                        temp += a[i__ + j * a_dim1] * x[ix];
                    }
                    x[jx] = temp;
                    jx += *incx;
                }
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
