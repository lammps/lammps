#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int zgemm_(char *transa, char *transb, integer *m, integer *n, integer *k, doublecomplex *alpha,
           doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublecomplex *beta,
           doublecomplex *c__, integer *ldc, ftnlen transa_len, ftnlen transb_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4, i__5,
        i__6;
    doublecomplex z__1, z__2, z__3, z__4;
    void d_lmp_cnjg(doublecomplex *, doublecomplex *);
    integer i__, j, l, info;
    logical nota, notb;
    doublecomplex temp;
    logical conja, conjb;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer nrowa, nrowb;
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
    nota = lsame_(transa, (char *)"N", (ftnlen)1, (ftnlen)1);
    notb = lsame_(transb, (char *)"N", (ftnlen)1, (ftnlen)1);
    conja = lsame_(transa, (char *)"C", (ftnlen)1, (ftnlen)1);
    conjb = lsame_(transb, (char *)"C", (ftnlen)1, (ftnlen)1);
    if (nota) {
        nrowa = *m;
    } else {
        nrowa = *k;
    }
    if (notb) {
        nrowb = *k;
    } else {
        nrowb = *n;
    }
    info = 0;
    if (!nota && !conja && !lsame_(transa, (char *)"T", (ftnlen)1, (ftnlen)1)) {
        info = 1;
    } else if (!notb && !conjb && !lsame_(transb, (char *)"T", (ftnlen)1, (ftnlen)1)) {
        info = 2;
    } else if (*m < 0) {
        info = 3;
    } else if (*n < 0) {
        info = 4;
    } else if (*k < 0) {
        info = 5;
    } else if (*lda < max(1, nrowa)) {
        info = 8;
    } else if (*ldb < max(1, nrowb)) {
        info = 10;
    } else if (*ldc < max(1, *m)) {
        info = 13;
    }
    if (info != 0) {
        xerbla_((char *)"ZGEMM ", &info, (ftnlen)6);
        return 0;
    }
    if (*m == 0 || *n == 0 ||
        (alpha->r == 0. && alpha->i == 0. || *k == 0) && (beta->r == 1. && beta->i == 0.)) {
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
    if (notb) {
        if (nota) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (beta->r == 0. && beta->i == 0.) {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * c_dim1;
                        c__[i__3].r = 0., c__[i__3].i = 0.;
                    }
                } else if (beta->r != 1. || beta->i != 0.) {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + j * c_dim1;
                        z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4].i,
                        z__1.i = beta->r * c__[i__4].i + beta->i * c__[i__4].r;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l) {
                    i__3 = l + j * b_dim1;
                    z__1.r = alpha->r * b[i__3].r - alpha->i * b[i__3].i,
                    z__1.i = alpha->r * b[i__3].i + alpha->i * b[i__3].r;
                    temp.r = z__1.r, temp.i = z__1.i;
                    i__3 = *m;
                    for (i__ = 1; i__ <= i__3; ++i__) {
                        i__4 = i__ + j * c_dim1;
                        i__5 = i__ + j * c_dim1;
                        i__6 = i__ + l * a_dim1;
                        z__2.r = temp.r * a[i__6].r - temp.i * a[i__6].i,
                        z__2.i = temp.r * a[i__6].i + temp.i * a[i__6].r;
                        z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5].i + z__2.i;
                        c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
                    }
                }
            }
        } else if (conja) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp.r = 0., temp.i = 0.;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l) {
                        d_lmp_cnjg(&z__3, &a[l + i__ * a_dim1]);
                        i__4 = l + j * b_dim1;
                        z__2.r = z__3.r * b[i__4].r - z__3.i * b[i__4].i,
                        z__2.i = z__3.r * b[i__4].i + z__3.i * b[i__4].r;
                        z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                        temp.r = z__1.r, temp.i = z__1.i;
                    }
                    if (beta->r == 0. && beta->i == 0.) {
                        i__3 = i__ + j * c_dim1;
                        z__1.r = alpha->r * temp.r - alpha->i * temp.i,
                        z__1.i = alpha->r * temp.i + alpha->i * temp.r;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    } else {
                        i__3 = i__ + j * c_dim1;
                        z__2.r = alpha->r * temp.r - alpha->i * temp.i,
                        z__2.i = alpha->r * temp.i + alpha->i * temp.r;
                        i__4 = i__ + j * c_dim1;
                        z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4].i,
                        z__3.i = beta->r * c__[i__4].i + beta->i * c__[i__4].r;
                        z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                }
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp.r = 0., temp.i = 0.;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l) {
                        i__4 = l + i__ * a_dim1;
                        i__5 = l + j * b_dim1;
                        z__2.r = a[i__4].r * b[i__5].r - a[i__4].i * b[i__5].i,
                        z__2.i = a[i__4].r * b[i__5].i + a[i__4].i * b[i__5].r;
                        z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                        temp.r = z__1.r, temp.i = z__1.i;
                    }
                    if (beta->r == 0. && beta->i == 0.) {
                        i__3 = i__ + j * c_dim1;
                        z__1.r = alpha->r * temp.r - alpha->i * temp.i,
                        z__1.i = alpha->r * temp.i + alpha->i * temp.r;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    } else {
                        i__3 = i__ + j * c_dim1;
                        z__2.r = alpha->r * temp.r - alpha->i * temp.i,
                        z__2.i = alpha->r * temp.i + alpha->i * temp.r;
                        i__4 = i__ + j * c_dim1;
                        z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4].i,
                        z__3.i = beta->r * c__[i__4].i + beta->i * c__[i__4].r;
                        z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                }
            }
        }
    } else if (nota) {
        if (conjb) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (beta->r == 0. && beta->i == 0.) {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * c_dim1;
                        c__[i__3].r = 0., c__[i__3].i = 0.;
                    }
                } else if (beta->r != 1. || beta->i != 0.) {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + j * c_dim1;
                        z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4].i,
                        z__1.i = beta->r * c__[i__4].i + beta->i * c__[i__4].r;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l) {
                    d_lmp_cnjg(&z__2, &b[j + l * b_dim1]);
                    z__1.r = alpha->r * z__2.r - alpha->i * z__2.i,
                    z__1.i = alpha->r * z__2.i + alpha->i * z__2.r;
                    temp.r = z__1.r, temp.i = z__1.i;
                    i__3 = *m;
                    for (i__ = 1; i__ <= i__3; ++i__) {
                        i__4 = i__ + j * c_dim1;
                        i__5 = i__ + j * c_dim1;
                        i__6 = i__ + l * a_dim1;
                        z__2.r = temp.r * a[i__6].r - temp.i * a[i__6].i,
                        z__2.i = temp.r * a[i__6].i + temp.i * a[i__6].r;
                        z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5].i + z__2.i;
                        c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
                    }
                }
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (beta->r == 0. && beta->i == 0.) {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * c_dim1;
                        c__[i__3].r = 0., c__[i__3].i = 0.;
                    }
                } else if (beta->r != 1. || beta->i != 0.) {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + j * c_dim1;
                        z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4].i,
                        z__1.i = beta->r * c__[i__4].i + beta->i * c__[i__4].r;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l) {
                    i__3 = j + l * b_dim1;
                    z__1.r = alpha->r * b[i__3].r - alpha->i * b[i__3].i,
                    z__1.i = alpha->r * b[i__3].i + alpha->i * b[i__3].r;
                    temp.r = z__1.r, temp.i = z__1.i;
                    i__3 = *m;
                    for (i__ = 1; i__ <= i__3; ++i__) {
                        i__4 = i__ + j * c_dim1;
                        i__5 = i__ + j * c_dim1;
                        i__6 = i__ + l * a_dim1;
                        z__2.r = temp.r * a[i__6].r - temp.i * a[i__6].i,
                        z__2.i = temp.r * a[i__6].i + temp.i * a[i__6].r;
                        z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5].i + z__2.i;
                        c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
                    }
                }
            }
        }
    } else if (conja) {
        if (conjb) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp.r = 0., temp.i = 0.;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l) {
                        d_lmp_cnjg(&z__3, &a[l + i__ * a_dim1]);
                        d_lmp_cnjg(&z__4, &b[j + l * b_dim1]);
                        z__2.r = z__3.r * z__4.r - z__3.i * z__4.i,
                        z__2.i = z__3.r * z__4.i + z__3.i * z__4.r;
                        z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                        temp.r = z__1.r, temp.i = z__1.i;
                    }
                    if (beta->r == 0. && beta->i == 0.) {
                        i__3 = i__ + j * c_dim1;
                        z__1.r = alpha->r * temp.r - alpha->i * temp.i,
                        z__1.i = alpha->r * temp.i + alpha->i * temp.r;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    } else {
                        i__3 = i__ + j * c_dim1;
                        z__2.r = alpha->r * temp.r - alpha->i * temp.i,
                        z__2.i = alpha->r * temp.i + alpha->i * temp.r;
                        i__4 = i__ + j * c_dim1;
                        z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4].i,
                        z__3.i = beta->r * c__[i__4].i + beta->i * c__[i__4].r;
                        z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                }
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp.r = 0., temp.i = 0.;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l) {
                        d_lmp_cnjg(&z__3, &a[l + i__ * a_dim1]);
                        i__4 = j + l * b_dim1;
                        z__2.r = z__3.r * b[i__4].r - z__3.i * b[i__4].i,
                        z__2.i = z__3.r * b[i__4].i + z__3.i * b[i__4].r;
                        z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                        temp.r = z__1.r, temp.i = z__1.i;
                    }
                    if (beta->r == 0. && beta->i == 0.) {
                        i__3 = i__ + j * c_dim1;
                        z__1.r = alpha->r * temp.r - alpha->i * temp.i,
                        z__1.i = alpha->r * temp.i + alpha->i * temp.r;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    } else {
                        i__3 = i__ + j * c_dim1;
                        z__2.r = alpha->r * temp.r - alpha->i * temp.i,
                        z__2.i = alpha->r * temp.i + alpha->i * temp.r;
                        i__4 = i__ + j * c_dim1;
                        z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4].i,
                        z__3.i = beta->r * c__[i__4].i + beta->i * c__[i__4].r;
                        z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                }
            }
        }
    } else {
        if (conjb) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp.r = 0., temp.i = 0.;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l) {
                        i__4 = l + i__ * a_dim1;
                        d_lmp_cnjg(&z__3, &b[j + l * b_dim1]);
                        z__2.r = a[i__4].r * z__3.r - a[i__4].i * z__3.i,
                        z__2.i = a[i__4].r * z__3.i + a[i__4].i * z__3.r;
                        z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                        temp.r = z__1.r, temp.i = z__1.i;
                    }
                    if (beta->r == 0. && beta->i == 0.) {
                        i__3 = i__ + j * c_dim1;
                        z__1.r = alpha->r * temp.r - alpha->i * temp.i,
                        z__1.i = alpha->r * temp.i + alpha->i * temp.r;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    } else {
                        i__3 = i__ + j * c_dim1;
                        z__2.r = alpha->r * temp.r - alpha->i * temp.i,
                        z__2.i = alpha->r * temp.i + alpha->i * temp.r;
                        i__4 = i__ + j * c_dim1;
                        z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4].i,
                        z__3.i = beta->r * c__[i__4].i + beta->i * c__[i__4].r;
                        z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                }
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp.r = 0., temp.i = 0.;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l) {
                        i__4 = l + i__ * a_dim1;
                        i__5 = j + l * b_dim1;
                        z__2.r = a[i__4].r * b[i__5].r - a[i__4].i * b[i__5].i,
                        z__2.i = a[i__4].r * b[i__5].i + a[i__4].i * b[i__5].r;
                        z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                        temp.r = z__1.r, temp.i = z__1.i;
                    }
                    if (beta->r == 0. && beta->i == 0.) {
                        i__3 = i__ + j * c_dim1;
                        z__1.r = alpha->r * temp.r - alpha->i * temp.i,
                        z__1.i = alpha->r * temp.i + alpha->i * temp.r;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    } else {
                        i__3 = i__ + j * c_dim1;
                        z__2.r = alpha->r * temp.r - alpha->i * temp.i,
                        z__2.i = alpha->r * temp.i + alpha->i * temp.r;
                        i__4 = i__ + j * c_dim1;
                        z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4].i,
                        z__3.i = beta->r * c__[i__4].i + beta->i * c__[i__4].r;
                        z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
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
