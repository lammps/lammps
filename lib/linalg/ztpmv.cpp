#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int ztpmv_(char *uplo, char *trans, char *diag, integer *n, doublecomplex *ap, doublecomplex *x,
           integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len)
{
    integer i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3;
    void d_lmp_cnjg(doublecomplex *, doublecomplex *);
    integer i__, j, k, kk, ix, jx, kx, info;
    doublecomplex temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int xerbla_(char *, integer *, ftnlen);
    logical noconj, nounit;
    --x;
    --ap;
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
    } else if (*incx == 0) {
        info = 7;
    }
    if (info != 0) {
        xerbla_((char *)"ZTPMV ", &info, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    noconj = lsame_(trans, (char *)"T", (ftnlen)1, (ftnlen)1);
    nounit = lsame_(diag, (char *)"N", (ftnlen)1, (ftnlen)1);
    if (*incx <= 0) {
        kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
        kx = 1;
    }
    if (lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1)) {
        if (lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1)) {
            kk = 1;
            if (*incx == 1) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = j;
                    if (x[i__2].r != 0. || x[i__2].i != 0.) {
                        i__2 = j;
                        temp.r = x[i__2].r, temp.i = x[i__2].i;
                        k = kk;
                        i__2 = j - 1;
                        for (i__ = 1; i__ <= i__2; ++i__) {
                            i__3 = i__;
                            i__4 = i__;
                            i__5 = k;
                            z__2.r = temp.r * ap[i__5].r - temp.i * ap[i__5].i,
                            z__2.i = temp.r * ap[i__5].i + temp.i * ap[i__5].r;
                            z__1.r = x[i__4].r + z__2.r, z__1.i = x[i__4].i + z__2.i;
                            x[i__3].r = z__1.r, x[i__3].i = z__1.i;
                            ++k;
                        }
                        if (nounit) {
                            i__2 = j;
                            i__3 = j;
                            i__4 = kk + j - 1;
                            z__1.r = x[i__3].r * ap[i__4].r - x[i__3].i * ap[i__4].i,
                            z__1.i = x[i__3].r * ap[i__4].i + x[i__3].i * ap[i__4].r;
                            x[i__2].r = z__1.r, x[i__2].i = z__1.i;
                        }
                    }
                    kk += j;
                }
            } else {
                jx = kx;
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = jx;
                    if (x[i__2].r != 0. || x[i__2].i != 0.) {
                        i__2 = jx;
                        temp.r = x[i__2].r, temp.i = x[i__2].i;
                        ix = kx;
                        i__2 = kk + j - 2;
                        for (k = kk; k <= i__2; ++k) {
                            i__3 = ix;
                            i__4 = ix;
                            i__5 = k;
                            z__2.r = temp.r * ap[i__5].r - temp.i * ap[i__5].i,
                            z__2.i = temp.r * ap[i__5].i + temp.i * ap[i__5].r;
                            z__1.r = x[i__4].r + z__2.r, z__1.i = x[i__4].i + z__2.i;
                            x[i__3].r = z__1.r, x[i__3].i = z__1.i;
                            ix += *incx;
                        }
                        if (nounit) {
                            i__2 = jx;
                            i__3 = jx;
                            i__4 = kk + j - 1;
                            z__1.r = x[i__3].r * ap[i__4].r - x[i__3].i * ap[i__4].i,
                            z__1.i = x[i__3].r * ap[i__4].i + x[i__3].i * ap[i__4].r;
                            x[i__2].r = z__1.r, x[i__2].i = z__1.i;
                        }
                    }
                    jx += *incx;
                    kk += j;
                }
            }
        } else {
            kk = *n * (*n + 1) / 2;
            if (*incx == 1) {
                for (j = *n; j >= 1; --j) {
                    i__1 = j;
                    if (x[i__1].r != 0. || x[i__1].i != 0.) {
                        i__1 = j;
                        temp.r = x[i__1].r, temp.i = x[i__1].i;
                        k = kk;
                        i__1 = j + 1;
                        for (i__ = *n; i__ >= i__1; --i__) {
                            i__2 = i__;
                            i__3 = i__;
                            i__4 = k;
                            z__2.r = temp.r * ap[i__4].r - temp.i * ap[i__4].i,
                            z__2.i = temp.r * ap[i__4].i + temp.i * ap[i__4].r;
                            z__1.r = x[i__3].r + z__2.r, z__1.i = x[i__3].i + z__2.i;
                            x[i__2].r = z__1.r, x[i__2].i = z__1.i;
                            --k;
                        }
                        if (nounit) {
                            i__1 = j;
                            i__2 = j;
                            i__3 = kk - *n + j;
                            z__1.r = x[i__2].r * ap[i__3].r - x[i__2].i * ap[i__3].i,
                            z__1.i = x[i__2].r * ap[i__3].i + x[i__2].i * ap[i__3].r;
                            x[i__1].r = z__1.r, x[i__1].i = z__1.i;
                        }
                    }
                    kk -= *n - j + 1;
                }
            } else {
                kx += (*n - 1) * *incx;
                jx = kx;
                for (j = *n; j >= 1; --j) {
                    i__1 = jx;
                    if (x[i__1].r != 0. || x[i__1].i != 0.) {
                        i__1 = jx;
                        temp.r = x[i__1].r, temp.i = x[i__1].i;
                        ix = kx;
                        i__1 = kk - (*n - (j + 1));
                        for (k = kk; k >= i__1; --k) {
                            i__2 = ix;
                            i__3 = ix;
                            i__4 = k;
                            z__2.r = temp.r * ap[i__4].r - temp.i * ap[i__4].i,
                            z__2.i = temp.r * ap[i__4].i + temp.i * ap[i__4].r;
                            z__1.r = x[i__3].r + z__2.r, z__1.i = x[i__3].i + z__2.i;
                            x[i__2].r = z__1.r, x[i__2].i = z__1.i;
                            ix -= *incx;
                        }
                        if (nounit) {
                            i__1 = jx;
                            i__2 = jx;
                            i__3 = kk - *n + j;
                            z__1.r = x[i__2].r * ap[i__3].r - x[i__2].i * ap[i__3].i,
                            z__1.i = x[i__2].r * ap[i__3].i + x[i__2].i * ap[i__3].r;
                            x[i__1].r = z__1.r, x[i__1].i = z__1.i;
                        }
                    }
                    jx -= *incx;
                    kk -= *n - j + 1;
                }
            }
        }
    } else {
        if (lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1)) {
            kk = *n * (*n + 1) / 2;
            if (*incx == 1) {
                for (j = *n; j >= 1; --j) {
                    i__1 = j;
                    temp.r = x[i__1].r, temp.i = x[i__1].i;
                    k = kk - 1;
                    if (noconj) {
                        if (nounit) {
                            i__1 = kk;
                            z__1.r = temp.r * ap[i__1].r - temp.i * ap[i__1].i,
                            z__1.i = temp.r * ap[i__1].i + temp.i * ap[i__1].r;
                            temp.r = z__1.r, temp.i = z__1.i;
                        }
                        for (i__ = j - 1; i__ >= 1; --i__) {
                            i__1 = k;
                            i__2 = i__;
                            z__2.r = ap[i__1].r * x[i__2].r - ap[i__1].i * x[i__2].i,
                            z__2.i = ap[i__1].r * x[i__2].i + ap[i__1].i * x[i__2].r;
                            z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                            temp.r = z__1.r, temp.i = z__1.i;
                            --k;
                        }
                    } else {
                        if (nounit) {
                            d_lmp_cnjg(&z__2, &ap[kk]);
                            z__1.r = temp.r * z__2.r - temp.i * z__2.i,
                            z__1.i = temp.r * z__2.i + temp.i * z__2.r;
                            temp.r = z__1.r, temp.i = z__1.i;
                        }
                        for (i__ = j - 1; i__ >= 1; --i__) {
                            d_lmp_cnjg(&z__3, &ap[k]);
                            i__1 = i__;
                            z__2.r = z__3.r * x[i__1].r - z__3.i * x[i__1].i,
                            z__2.i = z__3.r * x[i__1].i + z__3.i * x[i__1].r;
                            z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                            temp.r = z__1.r, temp.i = z__1.i;
                            --k;
                        }
                    }
                    i__1 = j;
                    x[i__1].r = temp.r, x[i__1].i = temp.i;
                    kk -= j;
                }
            } else {
                jx = kx + (*n - 1) * *incx;
                for (j = *n; j >= 1; --j) {
                    i__1 = jx;
                    temp.r = x[i__1].r, temp.i = x[i__1].i;
                    ix = jx;
                    if (noconj) {
                        if (nounit) {
                            i__1 = kk;
                            z__1.r = temp.r * ap[i__1].r - temp.i * ap[i__1].i,
                            z__1.i = temp.r * ap[i__1].i + temp.i * ap[i__1].r;
                            temp.r = z__1.r, temp.i = z__1.i;
                        }
                        i__1 = kk - j + 1;
                        for (k = kk - 1; k >= i__1; --k) {
                            ix -= *incx;
                            i__2 = k;
                            i__3 = ix;
                            z__2.r = ap[i__2].r * x[i__3].r - ap[i__2].i * x[i__3].i,
                            z__2.i = ap[i__2].r * x[i__3].i + ap[i__2].i * x[i__3].r;
                            z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                            temp.r = z__1.r, temp.i = z__1.i;
                        }
                    } else {
                        if (nounit) {
                            d_lmp_cnjg(&z__2, &ap[kk]);
                            z__1.r = temp.r * z__2.r - temp.i * z__2.i,
                            z__1.i = temp.r * z__2.i + temp.i * z__2.r;
                            temp.r = z__1.r, temp.i = z__1.i;
                        }
                        i__1 = kk - j + 1;
                        for (k = kk - 1; k >= i__1; --k) {
                            ix -= *incx;
                            d_lmp_cnjg(&z__3, &ap[k]);
                            i__2 = ix;
                            z__2.r = z__3.r * x[i__2].r - z__3.i * x[i__2].i,
                            z__2.i = z__3.r * x[i__2].i + z__3.i * x[i__2].r;
                            z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                            temp.r = z__1.r, temp.i = z__1.i;
                        }
                    }
                    i__1 = jx;
                    x[i__1].r = temp.r, x[i__1].i = temp.i;
                    jx -= *incx;
                    kk -= j;
                }
            }
        } else {
            kk = 1;
            if (*incx == 1) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = j;
                    temp.r = x[i__2].r, temp.i = x[i__2].i;
                    k = kk + 1;
                    if (noconj) {
                        if (nounit) {
                            i__2 = kk;
                            z__1.r = temp.r * ap[i__2].r - temp.i * ap[i__2].i,
                            z__1.i = temp.r * ap[i__2].i + temp.i * ap[i__2].r;
                            temp.r = z__1.r, temp.i = z__1.i;
                        }
                        i__2 = *n;
                        for (i__ = j + 1; i__ <= i__2; ++i__) {
                            i__3 = k;
                            i__4 = i__;
                            z__2.r = ap[i__3].r * x[i__4].r - ap[i__3].i * x[i__4].i,
                            z__2.i = ap[i__3].r * x[i__4].i + ap[i__3].i * x[i__4].r;
                            z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                            temp.r = z__1.r, temp.i = z__1.i;
                            ++k;
                        }
                    } else {
                        if (nounit) {
                            d_lmp_cnjg(&z__2, &ap[kk]);
                            z__1.r = temp.r * z__2.r - temp.i * z__2.i,
                            z__1.i = temp.r * z__2.i + temp.i * z__2.r;
                            temp.r = z__1.r, temp.i = z__1.i;
                        }
                        i__2 = *n;
                        for (i__ = j + 1; i__ <= i__2; ++i__) {
                            d_lmp_cnjg(&z__3, &ap[k]);
                            i__3 = i__;
                            z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i,
                            z__2.i = z__3.r * x[i__3].i + z__3.i * x[i__3].r;
                            z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                            temp.r = z__1.r, temp.i = z__1.i;
                            ++k;
                        }
                    }
                    i__2 = j;
                    x[i__2].r = temp.r, x[i__2].i = temp.i;
                    kk += *n - j + 1;
                }
            } else {
                jx = kx;
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = jx;
                    temp.r = x[i__2].r, temp.i = x[i__2].i;
                    ix = jx;
                    if (noconj) {
                        if (nounit) {
                            i__2 = kk;
                            z__1.r = temp.r * ap[i__2].r - temp.i * ap[i__2].i,
                            z__1.i = temp.r * ap[i__2].i + temp.i * ap[i__2].r;
                            temp.r = z__1.r, temp.i = z__1.i;
                        }
                        i__2 = kk + *n - j;
                        for (k = kk + 1; k <= i__2; ++k) {
                            ix += *incx;
                            i__3 = k;
                            i__4 = ix;
                            z__2.r = ap[i__3].r * x[i__4].r - ap[i__3].i * x[i__4].i,
                            z__2.i = ap[i__3].r * x[i__4].i + ap[i__3].i * x[i__4].r;
                            z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                            temp.r = z__1.r, temp.i = z__1.i;
                        }
                    } else {
                        if (nounit) {
                            d_lmp_cnjg(&z__2, &ap[kk]);
                            z__1.r = temp.r * z__2.r - temp.i * z__2.i,
                            z__1.i = temp.r * z__2.i + temp.i * z__2.r;
                            temp.r = z__1.r, temp.i = z__1.i;
                        }
                        i__2 = kk + *n - j;
                        for (k = kk + 1; k <= i__2; ++k) {
                            ix += *incx;
                            d_lmp_cnjg(&z__3, &ap[k]);
                            i__3 = ix;
                            z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i,
                            z__2.i = z__3.r * x[i__3].i + z__3.i * x[i__3].r;
                            z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                            temp.r = z__1.r, temp.i = z__1.i;
                        }
                    }
                    i__2 = jx;
                    x[i__2].r = temp.r, x[i__2].i = temp.i;
                    jx += *incx;
                    kk += *n - j + 1;
                }
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
