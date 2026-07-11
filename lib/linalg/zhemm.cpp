#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int zhemm_(char *side, char *uplo, integer *m, integer *n, doublecomplex *alpha, doublecomplex *a,
           integer *lda, doublecomplex *b, integer *ldb, doublecomplex *beta, doublecomplex *c__,
           integer *ldc, ftnlen side_len, ftnlen uplo_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4, i__5,
        i__6;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4, z__5;
    void d_lmp_cnjg(doublecomplex *, doublecomplex *);
    integer i__, j, k, info;
    doublecomplex temp1, temp2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer nrowa;
    logical upper;
    extern int xerbla_(char *, integer *, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    if (lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        nrowa = *m;
    } else {
        nrowa = *n;
    }
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    info = 0;
    if (!lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1) && !lsame_(side, (char *)"R", (ftnlen)1, (ftnlen)1)) {
        info = 1;
    } else if (!upper && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        info = 2;
    } else if (*m < 0) {
        info = 3;
    } else if (*n < 0) {
        info = 4;
    } else if (*lda < max(1, nrowa)) {
        info = 7;
    } else if (*ldb < max(1, *m)) {
        info = 9;
    } else if (*ldc < max(1, *m)) {
        info = 12;
    }
    if (info != 0) {
        xerbla_((char *)"ZHEMM ", &info, (ftnlen)6);
        return 0;
    }
    if (*m == 0 || *n == 0 ||
        alpha->r == 0. && alpha->i == 0. && (beta->r == 1. && beta->i == 0.)) {
        return 0;
    }
    if (alpha->r == 0. && alpha->i == 0.) {
        if (beta->r == 0. && beta->i == 0.) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    i__3 = i__ + j * c_dim1;
                    c__[i__3].r = 0., c__[i__3].i = 0.;
                }
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    i__3 = i__ + j * c_dim1;
                    i__4 = i__ + j * c_dim1;
                    z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4].i,
                    z__1.i = beta->r * c__[i__4].i + beta->i * c__[i__4].r;
                    c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                }
            }
        }
        return 0;
    }
    if (lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        if (upper) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    i__3 = i__ + j * b_dim1;
                    z__1.r = alpha->r * b[i__3].r - alpha->i * b[i__3].i,
                    z__1.i = alpha->r * b[i__3].i + alpha->i * b[i__3].r;
                    temp1.r = z__1.r, temp1.i = z__1.i;
                    temp2.r = 0., temp2.i = 0.;
                    i__3 = i__ - 1;
                    for (k = 1; k <= i__3; ++k) {
                        i__4 = k + j * c_dim1;
                        i__5 = k + j * c_dim1;
                        i__6 = k + i__ * a_dim1;
                        z__2.r = temp1.r * a[i__6].r - temp1.i * a[i__6].i,
                        z__2.i = temp1.r * a[i__6].i + temp1.i * a[i__6].r;
                        z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5].i + z__2.i;
                        c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
                        i__4 = k + j * b_dim1;
                        d_lmp_cnjg(&z__3, &a[k + i__ * a_dim1]);
                        z__2.r = b[i__4].r * z__3.r - b[i__4].i * z__3.i,
                        z__2.i = b[i__4].r * z__3.i + b[i__4].i * z__3.r;
                        z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
                        temp2.r = z__1.r, temp2.i = z__1.i;
                    }
                    if (beta->r == 0. && beta->i == 0.) {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + i__ * a_dim1;
                        d__1 = a[i__4].r;
                        z__2.r = d__1 * temp1.r, z__2.i = d__1 * temp1.i;
                        z__3.r = alpha->r * temp2.r - alpha->i * temp2.i,
                        z__3.i = alpha->r * temp2.i + alpha->i * temp2.r;
                        z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    } else {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + j * c_dim1;
                        z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4].i,
                        z__3.i = beta->r * c__[i__4].i + beta->i * c__[i__4].r;
                        i__5 = i__ + i__ * a_dim1;
                        d__1 = a[i__5].r;
                        z__4.r = d__1 * temp1.r, z__4.i = d__1 * temp1.i;
                        z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
                        z__5.r = alpha->r * temp2.r - alpha->i * temp2.i,
                        z__5.i = alpha->r * temp2.i + alpha->i * temp2.r;
                        z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                }
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                for (i__ = *m; i__ >= 1; --i__) {
                    i__2 = i__ + j * b_dim1;
                    z__1.r = alpha->r * b[i__2].r - alpha->i * b[i__2].i,
                    z__1.i = alpha->r * b[i__2].i + alpha->i * b[i__2].r;
                    temp1.r = z__1.r, temp1.i = z__1.i;
                    temp2.r = 0., temp2.i = 0.;
                    i__2 = *m;
                    for (k = i__ + 1; k <= i__2; ++k) {
                        i__3 = k + j * c_dim1;
                        i__4 = k + j * c_dim1;
                        i__5 = k + i__ * a_dim1;
                        z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i,
                        z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5].r;
                        z__1.r = c__[i__4].r + z__2.r, z__1.i = c__[i__4].i + z__2.i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                        i__3 = k + j * b_dim1;
                        d_lmp_cnjg(&z__3, &a[k + i__ * a_dim1]);
                        z__2.r = b[i__3].r * z__3.r - b[i__3].i * z__3.i,
                        z__2.i = b[i__3].r * z__3.i + b[i__3].i * z__3.r;
                        z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
                        temp2.r = z__1.r, temp2.i = z__1.i;
                    }
                    if (beta->r == 0. && beta->i == 0.) {
                        i__2 = i__ + j * c_dim1;
                        i__3 = i__ + i__ * a_dim1;
                        d__1 = a[i__3].r;
                        z__2.r = d__1 * temp1.r, z__2.i = d__1 * temp1.i;
                        z__3.r = alpha->r * temp2.r - alpha->i * temp2.i,
                        z__3.i = alpha->r * temp2.i + alpha->i * temp2.r;
                        z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                        c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
                    } else {
                        i__2 = i__ + j * c_dim1;
                        i__3 = i__ + j * c_dim1;
                        z__3.r = beta->r * c__[i__3].r - beta->i * c__[i__3].i,
                        z__3.i = beta->r * c__[i__3].i + beta->i * c__[i__3].r;
                        i__4 = i__ + i__ * a_dim1;
                        d__1 = a[i__4].r;
                        z__4.r = d__1 * temp1.r, z__4.i = d__1 * temp1.i;
                        z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
                        z__5.r = alpha->r * temp2.r - alpha->i * temp2.i,
                        z__5.i = alpha->r * temp2.i + alpha->i * temp2.r;
                        z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
                        c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
                    }
                }
            }
        }
    } else {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__2 = j + j * a_dim1;
            d__1 = a[i__2].r;
            z__1.r = d__1 * alpha->r, z__1.i = d__1 * alpha->i;
            temp1.r = z__1.r, temp1.i = z__1.i;
            if (beta->r == 0. && beta->i == 0.) {
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    i__3 = i__ + j * c_dim1;
                    i__4 = i__ + j * b_dim1;
                    z__1.r = temp1.r * b[i__4].r - temp1.i * b[i__4].i,
                    z__1.i = temp1.r * b[i__4].i + temp1.i * b[i__4].r;
                    c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                }
            } else {
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    i__3 = i__ + j * c_dim1;
                    i__4 = i__ + j * c_dim1;
                    z__2.r = beta->r * c__[i__4].r - beta->i * c__[i__4].i,
                    z__2.i = beta->r * c__[i__4].i + beta->i * c__[i__4].r;
                    i__5 = i__ + j * b_dim1;
                    z__3.r = temp1.r * b[i__5].r - temp1.i * b[i__5].i,
                    z__3.i = temp1.r * b[i__5].i + temp1.i * b[i__5].r;
                    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                    c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                }
            }
            i__2 = j - 1;
            for (k = 1; k <= i__2; ++k) {
                if (upper) {
                    i__3 = k + j * a_dim1;
                    z__1.r = alpha->r * a[i__3].r - alpha->i * a[i__3].i,
                    z__1.i = alpha->r * a[i__3].i + alpha->i * a[i__3].r;
                    temp1.r = z__1.r, temp1.i = z__1.i;
                } else {
                    d_lmp_cnjg(&z__2, &a[j + k * a_dim1]);
                    z__1.r = alpha->r * z__2.r - alpha->i * z__2.i,
                    z__1.i = alpha->r * z__2.i + alpha->i * z__2.r;
                    temp1.r = z__1.r, temp1.i = z__1.i;
                }
                i__3 = *m;
                for (i__ = 1; i__ <= i__3; ++i__) {
                    i__4 = i__ + j * c_dim1;
                    i__5 = i__ + j * c_dim1;
                    i__6 = i__ + k * b_dim1;
                    z__2.r = temp1.r * b[i__6].r - temp1.i * b[i__6].i,
                    z__2.i = temp1.r * b[i__6].i + temp1.i * b[i__6].r;
                    z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5].i + z__2.i;
                    c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
                }
            }
            i__2 = *n;
            for (k = j + 1; k <= i__2; ++k) {
                if (upper) {
                    d_lmp_cnjg(&z__2, &a[j + k * a_dim1]);
                    z__1.r = alpha->r * z__2.r - alpha->i * z__2.i,
                    z__1.i = alpha->r * z__2.i + alpha->i * z__2.r;
                    temp1.r = z__1.r, temp1.i = z__1.i;
                } else {
                    i__3 = k + j * a_dim1;
                    z__1.r = alpha->r * a[i__3].r - alpha->i * a[i__3].i,
                    z__1.i = alpha->r * a[i__3].i + alpha->i * a[i__3].r;
                    temp1.r = z__1.r, temp1.i = z__1.i;
                }
                i__3 = *m;
                for (i__ = 1; i__ <= i__3; ++i__) {
                    i__4 = i__ + j * c_dim1;
                    i__5 = i__ + j * c_dim1;
                    i__6 = i__ + k * b_dim1;
                    z__2.r = temp1.r * b[i__6].r - temp1.i * b[i__6].i,
                    z__2.i = temp1.r * b[i__6].i + temp1.i * b[i__6].r;
                    z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5].i + z__2.i;
                    c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
                }
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
