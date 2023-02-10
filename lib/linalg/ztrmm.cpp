#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int ztrmm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n,
           doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb,
           ftnlen side_len, ftnlen uplo_len, ftnlen transa_len, ftnlen diag_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublecomplex z__1, z__2, z__3;
    void d_lmp_cnjg(doublecomplex *, doublecomplex *);
    integer i__, j, k, info;
    doublecomplex temp;
    logical lside;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer nrowa;
    logical upper;
    extern int xerbla_(char *, integer *, ftnlen);
    logical noconj, nounit;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    lside = lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1);
    if (lside) {
        nrowa = *m;
    } else {
        nrowa = *n;
    }
    noconj = lsame_(transa, (char *)"T", (ftnlen)1, (ftnlen)1);
    nounit = lsame_(diag, (char *)"N", (ftnlen)1, (ftnlen)1);
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    info = 0;
    if (!lside && !lsame_(side, (char *)"R", (ftnlen)1, (ftnlen)1)) {
        info = 1;
    } else if (!upper && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        info = 2;
    } else if (!lsame_(transa, (char *)"N", (ftnlen)1, (ftnlen)1) &&
               !lsame_(transa, (char *)"T", (ftnlen)1, (ftnlen)1) &&
               !lsame_(transa, (char *)"C", (ftnlen)1, (ftnlen)1)) {
        info = 3;
    } else if (!lsame_(diag, (char *)"U", (ftnlen)1, (ftnlen)1) &&
               !lsame_(diag, (char *)"N", (ftnlen)1, (ftnlen)1)) {
        info = 4;
    } else if (*m < 0) {
        info = 5;
    } else if (*n < 0) {
        info = 6;
    } else if (*lda < max(1, nrowa)) {
        info = 9;
    } else if (*ldb < max(1, *m)) {
        info = 11;
    }
    if (info != 0) {
        xerbla_((char *)"ZTRMM ", &info, (ftnlen)6);
        return 0;
    }
    if (*m == 0 || *n == 0) {
        return 0;
    }
    if (alpha->r == 0. && alpha->i == 0.) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__2 = *m;
            for (i__ = 1; i__ <= i__2; ++i__) {
                i__3 = i__ + j * b_dim1;
                b[i__3].r = 0., b[i__3].i = 0.;
            }
        }
        return 0;
    }
    if (lside) {
        if (lsame_(transa, (char *)"N", (ftnlen)1, (ftnlen)1)) {
            if (upper) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *m;
                    for (k = 1; k <= i__2; ++k) {
                        i__3 = k + j * b_dim1;
                        if (b[i__3].r != 0. || b[i__3].i != 0.) {
                            i__3 = k + j * b_dim1;
                            z__1.r = alpha->r * b[i__3].r - alpha->i * b[i__3].i,
                            z__1.i = alpha->r * b[i__3].i + alpha->i * b[i__3].r;
                            temp.r = z__1.r, temp.i = z__1.i;
                            i__3 = k - 1;
                            for (i__ = 1; i__ <= i__3; ++i__) {
                                i__4 = i__ + j * b_dim1;
                                i__5 = i__ + j * b_dim1;
                                i__6 = i__ + k * a_dim1;
                                z__2.r = temp.r * a[i__6].r - temp.i * a[i__6].i,
                                z__2.i = temp.r * a[i__6].i + temp.i * a[i__6].r;
                                z__1.r = b[i__5].r + z__2.r, z__1.i = b[i__5].i + z__2.i;
                                b[i__4].r = z__1.r, b[i__4].i = z__1.i;
                            }
                            if (nounit) {
                                i__3 = k + k * a_dim1;
                                z__1.r = temp.r * a[i__3].r - temp.i * a[i__3].i,
                                z__1.i = temp.r * a[i__3].i + temp.i * a[i__3].r;
                                temp.r = z__1.r, temp.i = z__1.i;
                            }
                            i__3 = k + j * b_dim1;
                            b[i__3].r = temp.r, b[i__3].i = temp.i;
                        }
                    }
                }
            } else {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    for (k = *m; k >= 1; --k) {
                        i__2 = k + j * b_dim1;
                        if (b[i__2].r != 0. || b[i__2].i != 0.) {
                            i__2 = k + j * b_dim1;
                            z__1.r = alpha->r * b[i__2].r - alpha->i * b[i__2].i,
                            z__1.i = alpha->r * b[i__2].i + alpha->i * b[i__2].r;
                            temp.r = z__1.r, temp.i = z__1.i;
                            i__2 = k + j * b_dim1;
                            b[i__2].r = temp.r, b[i__2].i = temp.i;
                            if (nounit) {
                                i__2 = k + j * b_dim1;
                                i__3 = k + j * b_dim1;
                                i__4 = k + k * a_dim1;
                                z__1.r = b[i__3].r * a[i__4].r - b[i__3].i * a[i__4].i,
                                z__1.i = b[i__3].r * a[i__4].i + b[i__3].i * a[i__4].r;
                                b[i__2].r = z__1.r, b[i__2].i = z__1.i;
                            }
                            i__2 = *m;
                            for (i__ = k + 1; i__ <= i__2; ++i__) {
                                i__3 = i__ + j * b_dim1;
                                i__4 = i__ + j * b_dim1;
                                i__5 = i__ + k * a_dim1;
                                z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i,
                                z__2.i = temp.r * a[i__5].i + temp.i * a[i__5].r;
                                z__1.r = b[i__4].r + z__2.r, z__1.i = b[i__4].i + z__2.i;
                                b[i__3].r = z__1.r, b[i__3].i = z__1.i;
                            }
                        }
                    }
                }
            }
        } else {
            if (upper) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    for (i__ = *m; i__ >= 1; --i__) {
                        i__2 = i__ + j * b_dim1;
                        temp.r = b[i__2].r, temp.i = b[i__2].i;
                        if (noconj) {
                            if (nounit) {
                                i__2 = i__ + i__ * a_dim1;
                                z__1.r = temp.r * a[i__2].r - temp.i * a[i__2].i,
                                z__1.i = temp.r * a[i__2].i + temp.i * a[i__2].r;
                                temp.r = z__1.r, temp.i = z__1.i;
                            }
                            i__2 = i__ - 1;
                            for (k = 1; k <= i__2; ++k) {
                                i__3 = k + i__ * a_dim1;
                                i__4 = k + j * b_dim1;
                                z__2.r = a[i__3].r * b[i__4].r - a[i__3].i * b[i__4].i,
                                z__2.i = a[i__3].r * b[i__4].i + a[i__3].i * b[i__4].r;
                                z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                                temp.r = z__1.r, temp.i = z__1.i;
                            }
                        } else {
                            if (nounit) {
                                d_lmp_cnjg(&z__2, &a[i__ + i__ * a_dim1]);
                                z__1.r = temp.r * z__2.r - temp.i * z__2.i,
                                z__1.i = temp.r * z__2.i + temp.i * z__2.r;
                                temp.r = z__1.r, temp.i = z__1.i;
                            }
                            i__2 = i__ - 1;
                            for (k = 1; k <= i__2; ++k) {
                                d_lmp_cnjg(&z__3, &a[k + i__ * a_dim1]);
                                i__3 = k + j * b_dim1;
                                z__2.r = z__3.r * b[i__3].r - z__3.i * b[i__3].i,
                                z__2.i = z__3.r * b[i__3].i + z__3.i * b[i__3].r;
                                z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                                temp.r = z__1.r, temp.i = z__1.i;
                            }
                        }
                        i__2 = i__ + j * b_dim1;
                        z__1.r = alpha->r * temp.r - alpha->i * temp.i,
                        z__1.i = alpha->r * temp.i + alpha->i * temp.r;
                        b[i__2].r = z__1.r, b[i__2].i = z__1.i;
                    }
                }
            } else {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * b_dim1;
                        temp.r = b[i__3].r, temp.i = b[i__3].i;
                        if (noconj) {
                            if (nounit) {
                                i__3 = i__ + i__ * a_dim1;
                                z__1.r = temp.r * a[i__3].r - temp.i * a[i__3].i,
                                z__1.i = temp.r * a[i__3].i + temp.i * a[i__3].r;
                                temp.r = z__1.r, temp.i = z__1.i;
                            }
                            i__3 = *m;
                            for (k = i__ + 1; k <= i__3; ++k) {
                                i__4 = k + i__ * a_dim1;
                                i__5 = k + j * b_dim1;
                                z__2.r = a[i__4].r * b[i__5].r - a[i__4].i * b[i__5].i,
                                z__2.i = a[i__4].r * b[i__5].i + a[i__4].i * b[i__5].r;
                                z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                                temp.r = z__1.r, temp.i = z__1.i;
                            }
                        } else {
                            if (nounit) {
                                d_lmp_cnjg(&z__2, &a[i__ + i__ * a_dim1]);
                                z__1.r = temp.r * z__2.r - temp.i * z__2.i,
                                z__1.i = temp.r * z__2.i + temp.i * z__2.r;
                                temp.r = z__1.r, temp.i = z__1.i;
                            }
                            i__3 = *m;
                            for (k = i__ + 1; k <= i__3; ++k) {
                                d_lmp_cnjg(&z__3, &a[k + i__ * a_dim1]);
                                i__4 = k + j * b_dim1;
                                z__2.r = z__3.r * b[i__4].r - z__3.i * b[i__4].i,
                                z__2.i = z__3.r * b[i__4].i + z__3.i * b[i__4].r;
                                z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                                temp.r = z__1.r, temp.i = z__1.i;
                            }
                        }
                        i__3 = i__ + j * b_dim1;
                        z__1.r = alpha->r * temp.r - alpha->i * temp.i,
                        z__1.i = alpha->r * temp.i + alpha->i * temp.r;
                        b[i__3].r = z__1.r, b[i__3].i = z__1.i;
                    }
                }
            }
        }
    } else {
        if (lsame_(transa, (char *)"N", (ftnlen)1, (ftnlen)1)) {
            if (upper) {
                for (j = *n; j >= 1; --j) {
                    temp.r = alpha->r, temp.i = alpha->i;
                    if (nounit) {
                        i__1 = j + j * a_dim1;
                        z__1.r = temp.r * a[i__1].r - temp.i * a[i__1].i,
                        z__1.i = temp.r * a[i__1].i + temp.i * a[i__1].r;
                        temp.r = z__1.r, temp.i = z__1.i;
                    }
                    i__1 = *m;
                    for (i__ = 1; i__ <= i__1; ++i__) {
                        i__2 = i__ + j * b_dim1;
                        i__3 = i__ + j * b_dim1;
                        z__1.r = temp.r * b[i__3].r - temp.i * b[i__3].i,
                        z__1.i = temp.r * b[i__3].i + temp.i * b[i__3].r;
                        b[i__2].r = z__1.r, b[i__2].i = z__1.i;
                    }
                    i__1 = j - 1;
                    for (k = 1; k <= i__1; ++k) {
                        i__2 = k + j * a_dim1;
                        if (a[i__2].r != 0. || a[i__2].i != 0.) {
                            i__2 = k + j * a_dim1;
                            z__1.r = alpha->r * a[i__2].r - alpha->i * a[i__2].i,
                            z__1.i = alpha->r * a[i__2].i + alpha->i * a[i__2].r;
                            temp.r = z__1.r, temp.i = z__1.i;
                            i__2 = *m;
                            for (i__ = 1; i__ <= i__2; ++i__) {
                                i__3 = i__ + j * b_dim1;
                                i__4 = i__ + j * b_dim1;
                                i__5 = i__ + k * b_dim1;
                                z__2.r = temp.r * b[i__5].r - temp.i * b[i__5].i,
                                z__2.i = temp.r * b[i__5].i + temp.i * b[i__5].r;
                                z__1.r = b[i__4].r + z__2.r, z__1.i = b[i__4].i + z__2.i;
                                b[i__3].r = z__1.r, b[i__3].i = z__1.i;
                            }
                        }
                    }
                }
            } else {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    temp.r = alpha->r, temp.i = alpha->i;
                    if (nounit) {
                        i__2 = j + j * a_dim1;
                        z__1.r = temp.r * a[i__2].r - temp.i * a[i__2].i,
                        z__1.i = temp.r * a[i__2].i + temp.i * a[i__2].r;
                        temp.r = z__1.r, temp.i = z__1.i;
                    }
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * b_dim1;
                        i__4 = i__ + j * b_dim1;
                        z__1.r = temp.r * b[i__4].r - temp.i * b[i__4].i,
                        z__1.i = temp.r * b[i__4].i + temp.i * b[i__4].r;
                        b[i__3].r = z__1.r, b[i__3].i = z__1.i;
                    }
                    i__2 = *n;
                    for (k = j + 1; k <= i__2; ++k) {
                        i__3 = k + j * a_dim1;
                        if (a[i__3].r != 0. || a[i__3].i != 0.) {
                            i__3 = k + j * a_dim1;
                            z__1.r = alpha->r * a[i__3].r - alpha->i * a[i__3].i,
                            z__1.i = alpha->r * a[i__3].i + alpha->i * a[i__3].r;
                            temp.r = z__1.r, temp.i = z__1.i;
                            i__3 = *m;
                            for (i__ = 1; i__ <= i__3; ++i__) {
                                i__4 = i__ + j * b_dim1;
                                i__5 = i__ + j * b_dim1;
                                i__6 = i__ + k * b_dim1;
                                z__2.r = temp.r * b[i__6].r - temp.i * b[i__6].i,
                                z__2.i = temp.r * b[i__6].i + temp.i * b[i__6].r;
                                z__1.r = b[i__5].r + z__2.r, z__1.i = b[i__5].i + z__2.i;
                                b[i__4].r = z__1.r, b[i__4].i = z__1.i;
                            }
                        }
                    }
                }
            }
        } else {
            if (upper) {
                i__1 = *n;
                for (k = 1; k <= i__1; ++k) {
                    i__2 = k - 1;
                    for (j = 1; j <= i__2; ++j) {
                        i__3 = j + k * a_dim1;
                        if (a[i__3].r != 0. || a[i__3].i != 0.) {
                            if (noconj) {
                                i__3 = j + k * a_dim1;
                                z__1.r = alpha->r * a[i__3].r - alpha->i * a[i__3].i,
                                z__1.i = alpha->r * a[i__3].i + alpha->i * a[i__3].r;
                                temp.r = z__1.r, temp.i = z__1.i;
                            } else {
                                d_lmp_cnjg(&z__2, &a[j + k * a_dim1]);
                                z__1.r = alpha->r * z__2.r - alpha->i * z__2.i,
                                z__1.i = alpha->r * z__2.i + alpha->i * z__2.r;
                                temp.r = z__1.r, temp.i = z__1.i;
                            }
                            i__3 = *m;
                            for (i__ = 1; i__ <= i__3; ++i__) {
                                i__4 = i__ + j * b_dim1;
                                i__5 = i__ + j * b_dim1;
                                i__6 = i__ + k * b_dim1;
                                z__2.r = temp.r * b[i__6].r - temp.i * b[i__6].i,
                                z__2.i = temp.r * b[i__6].i + temp.i * b[i__6].r;
                                z__1.r = b[i__5].r + z__2.r, z__1.i = b[i__5].i + z__2.i;
                                b[i__4].r = z__1.r, b[i__4].i = z__1.i;
                            }
                        }
                    }
                    temp.r = alpha->r, temp.i = alpha->i;
                    if (nounit) {
                        if (noconj) {
                            i__2 = k + k * a_dim1;
                            z__1.r = temp.r * a[i__2].r - temp.i * a[i__2].i,
                            z__1.i = temp.r * a[i__2].i + temp.i * a[i__2].r;
                            temp.r = z__1.r, temp.i = z__1.i;
                        } else {
                            d_lmp_cnjg(&z__2, &a[k + k * a_dim1]);
                            z__1.r = temp.r * z__2.r - temp.i * z__2.i,
                            z__1.i = temp.r * z__2.i + temp.i * z__2.r;
                            temp.r = z__1.r, temp.i = z__1.i;
                        }
                    }
                    if (temp.r != 1. || temp.i != 0.) {
                        i__2 = *m;
                        for (i__ = 1; i__ <= i__2; ++i__) {
                            i__3 = i__ + k * b_dim1;
                            i__4 = i__ + k * b_dim1;
                            z__1.r = temp.r * b[i__4].r - temp.i * b[i__4].i,
                            z__1.i = temp.r * b[i__4].i + temp.i * b[i__4].r;
                            b[i__3].r = z__1.r, b[i__3].i = z__1.i;
                        }
                    }
                }
            } else {
                for (k = *n; k >= 1; --k) {
                    i__1 = *n;
                    for (j = k + 1; j <= i__1; ++j) {
                        i__2 = j + k * a_dim1;
                        if (a[i__2].r != 0. || a[i__2].i != 0.) {
                            if (noconj) {
                                i__2 = j + k * a_dim1;
                                z__1.r = alpha->r * a[i__2].r - alpha->i * a[i__2].i,
                                z__1.i = alpha->r * a[i__2].i + alpha->i * a[i__2].r;
                                temp.r = z__1.r, temp.i = z__1.i;
                            } else {
                                d_lmp_cnjg(&z__2, &a[j + k * a_dim1]);
                                z__1.r = alpha->r * z__2.r - alpha->i * z__2.i,
                                z__1.i = alpha->r * z__2.i + alpha->i * z__2.r;
                                temp.r = z__1.r, temp.i = z__1.i;
                            }
                            i__2 = *m;
                            for (i__ = 1; i__ <= i__2; ++i__) {
                                i__3 = i__ + j * b_dim1;
                                i__4 = i__ + j * b_dim1;
                                i__5 = i__ + k * b_dim1;
                                z__2.r = temp.r * b[i__5].r - temp.i * b[i__5].i,
                                z__2.i = temp.r * b[i__5].i + temp.i * b[i__5].r;
                                z__1.r = b[i__4].r + z__2.r, z__1.i = b[i__4].i + z__2.i;
                                b[i__3].r = z__1.r, b[i__3].i = z__1.i;
                            }
                        }
                    }
                    temp.r = alpha->r, temp.i = alpha->i;
                    if (nounit) {
                        if (noconj) {
                            i__1 = k + k * a_dim1;
                            z__1.r = temp.r * a[i__1].r - temp.i * a[i__1].i,
                            z__1.i = temp.r * a[i__1].i + temp.i * a[i__1].r;
                            temp.r = z__1.r, temp.i = z__1.i;
                        } else {
                            d_lmp_cnjg(&z__2, &a[k + k * a_dim1]);
                            z__1.r = temp.r * z__2.r - temp.i * z__2.i,
                            z__1.i = temp.r * z__2.i + temp.i * z__2.r;
                            temp.r = z__1.r, temp.i = z__1.i;
                        }
                    }
                    if (temp.r != 1. || temp.i != 0.) {
                        i__1 = *m;
                        for (i__ = 1; i__ <= i__1; ++i__) {
                            i__2 = i__ + k * b_dim1;
                            i__3 = i__ + k * b_dim1;
                            z__1.r = temp.r * b[i__3].r - temp.i * b[i__3].i,
                            z__1.i = temp.r * b[i__3].i + temp.i * b[i__3].r;
                            b[i__2].r = z__1.r, b[i__2].i = z__1.i;
                        }
                    }
                }
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
