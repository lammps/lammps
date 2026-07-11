#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int zherk_(char *uplo, char *trans, integer *n, integer *k, doublereal *alpha, doublecomplex *a,
           integer *lda, doublereal *beta, doublecomplex *c__, integer *ldc, ftnlen uplo_len,
           ftnlen trans_len)
{
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3;
    void d_lmp_cnjg(doublecomplex *, doublecomplex *);
    integer i__, j, l, info;
    doublecomplex temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer nrowa;
    doublereal rtemp;
    logical upper;
    extern int xerbla_(char *, integer *, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    if (lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1)) {
        nrowa = *n;
    } else {
        nrowa = *k;
    }
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    info = 0;
    if (!upper && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        info = 1;
    } else if (!lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1) &&
               !lsame_(trans, (char *)"C", (ftnlen)1, (ftnlen)1)) {
        info = 2;
    } else if (*n < 0) {
        info = 3;
    } else if (*k < 0) {
        info = 4;
    } else if (*lda < max(1, nrowa)) {
        info = 7;
    } else if (*ldc < max(1, *n)) {
        info = 10;
    }
    if (info != 0) {
        xerbla_((char *)"ZHERK ", &info, (ftnlen)6);
        return 0;
    }
    if (*n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
        return 0;
    }
    if (*alpha == 0.) {
        if (upper) {
            if (*beta == 0.) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * c_dim1;
                        c__[i__3].r = 0., c__[i__3].i = 0.;
                    }
                }
            } else {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = j - 1;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + j * c_dim1;
                        z__1.r = *beta * c__[i__4].r, z__1.i = *beta * c__[i__4].i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                    i__2 = j + j * c_dim1;
                    i__3 = j + j * c_dim1;
                    d__1 = *beta * c__[i__3].r;
                    c__[i__2].r = d__1, c__[i__2].i = 0.;
                }
            }
        } else {
            if (*beta == 0.) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * c_dim1;
                        c__[i__3].r = 0., c__[i__3].i = 0.;
                    }
                }
            } else {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = j + j * c_dim1;
                    i__3 = j + j * c_dim1;
                    d__1 = *beta * c__[i__3].r;
                    c__[i__2].r = d__1, c__[i__2].i = 0.;
                    i__2 = *n;
                    for (i__ = j + 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + j * c_dim1;
                        z__1.r = *beta * c__[i__4].r, z__1.i = *beta * c__[i__4].i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                }
            }
        }
        return 0;
    }
    if (lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1)) {
        if (upper) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (*beta == 0.) {
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * c_dim1;
                        c__[i__3].r = 0., c__[i__3].i = 0.;
                    }
                } else if (*beta != 1.) {
                    i__2 = j - 1;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + j * c_dim1;
                        z__1.r = *beta * c__[i__4].r, z__1.i = *beta * c__[i__4].i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                    i__2 = j + j * c_dim1;
                    i__3 = j + j * c_dim1;
                    d__1 = *beta * c__[i__3].r;
                    c__[i__2].r = d__1, c__[i__2].i = 0.;
                } else {
                    i__2 = j + j * c_dim1;
                    i__3 = j + j * c_dim1;
                    d__1 = c__[i__3].r;
                    c__[i__2].r = d__1, c__[i__2].i = 0.;
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l) {
                    i__3 = j + l * a_dim1;
                    if (a[i__3].r != 0. || a[i__3].i != 0.) {
                        d_lmp_cnjg(&z__2, &a[j + l * a_dim1]);
                        z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
                        temp.r = z__1.r, temp.i = z__1.i;
                        i__3 = j - 1;
                        for (i__ = 1; i__ <= i__3; ++i__) {
                            i__4 = i__ + j * c_dim1;
                            i__5 = i__ + j * c_dim1;
                            i__6 = i__ + l * a_dim1;
                            z__2.r = temp.r * a[i__6].r - temp.i * a[i__6].i,
                            z__2.i = temp.r * a[i__6].i + temp.i * a[i__6].r;
                            z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5].i + z__2.i;
                            c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
                        }
                        i__3 = j + j * c_dim1;
                        i__4 = j + j * c_dim1;
                        i__5 = i__ + l * a_dim1;
                        z__1.r = temp.r * a[i__5].r - temp.i * a[i__5].i,
                        z__1.i = temp.r * a[i__5].i + temp.i * a[i__5].r;
                        d__1 = c__[i__4].r + z__1.r;
                        c__[i__3].r = d__1, c__[i__3].i = 0.;
                    }
                }
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (*beta == 0.) {
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * c_dim1;
                        c__[i__3].r = 0., c__[i__3].i = 0.;
                    }
                } else if (*beta != 1.) {
                    i__2 = j + j * c_dim1;
                    i__3 = j + j * c_dim1;
                    d__1 = *beta * c__[i__3].r;
                    c__[i__2].r = d__1, c__[i__2].i = 0.;
                    i__2 = *n;
                    for (i__ = j + 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + j * c_dim1;
                        z__1.r = *beta * c__[i__4].r, z__1.i = *beta * c__[i__4].i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                } else {
                    i__2 = j + j * c_dim1;
                    i__3 = j + j * c_dim1;
                    d__1 = c__[i__3].r;
                    c__[i__2].r = d__1, c__[i__2].i = 0.;
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l) {
                    i__3 = j + l * a_dim1;
                    if (a[i__3].r != 0. || a[i__3].i != 0.) {
                        d_lmp_cnjg(&z__2, &a[j + l * a_dim1]);
                        z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
                        temp.r = z__1.r, temp.i = z__1.i;
                        i__3 = j + j * c_dim1;
                        i__4 = j + j * c_dim1;
                        i__5 = j + l * a_dim1;
                        z__1.r = temp.r * a[i__5].r - temp.i * a[i__5].i,
                        z__1.i = temp.r * a[i__5].i + temp.i * a[i__5].r;
                        d__1 = c__[i__4].r + z__1.r;
                        c__[i__3].r = d__1, c__[i__3].i = 0.;
                        i__3 = *n;
                        for (i__ = j + 1; i__ <= i__3; ++i__) {
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
        }
    } else {
        if (upper) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = j - 1;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp.r = 0., temp.i = 0.;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l) {
                        d_lmp_cnjg(&z__3, &a[l + i__ * a_dim1]);
                        i__4 = l + j * a_dim1;
                        z__2.r = z__3.r * a[i__4].r - z__3.i * a[i__4].i,
                        z__2.i = z__3.r * a[i__4].i + z__3.i * a[i__4].r;
                        z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                        temp.r = z__1.r, temp.i = z__1.i;
                    }
                    if (*beta == 0.) {
                        i__3 = i__ + j * c_dim1;
                        z__1.r = *alpha * temp.r, z__1.i = *alpha * temp.i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    } else {
                        i__3 = i__ + j * c_dim1;
                        z__2.r = *alpha * temp.r, z__2.i = *alpha * temp.i;
                        i__4 = i__ + j * c_dim1;
                        z__3.r = *beta * c__[i__4].r, z__3.i = *beta * c__[i__4].i;
                        z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                }
                rtemp = 0.;
                i__2 = *k;
                for (l = 1; l <= i__2; ++l) {
                    d_lmp_cnjg(&z__2, &a[l + j * a_dim1]);
                    i__3 = l + j * a_dim1;
                    z__1.r = z__2.r * a[i__3].r - z__2.i * a[i__3].i,
                    z__1.i = z__2.r * a[i__3].i + z__2.i * a[i__3].r;
                    rtemp += z__1.r;
                }
                if (*beta == 0.) {
                    i__2 = j + j * c_dim1;
                    d__1 = *alpha * rtemp;
                    c__[i__2].r = d__1, c__[i__2].i = 0.;
                } else {
                    i__2 = j + j * c_dim1;
                    i__3 = j + j * c_dim1;
                    d__1 = *alpha * rtemp + *beta * c__[i__3].r;
                    c__[i__2].r = d__1, c__[i__2].i = 0.;
                }
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                rtemp = 0.;
                i__2 = *k;
                for (l = 1; l <= i__2; ++l) {
                    d_lmp_cnjg(&z__2, &a[l + j * a_dim1]);
                    i__3 = l + j * a_dim1;
                    z__1.r = z__2.r * a[i__3].r - z__2.i * a[i__3].i,
                    z__1.i = z__2.r * a[i__3].i + z__2.i * a[i__3].r;
                    rtemp += z__1.r;
                }
                if (*beta == 0.) {
                    i__2 = j + j * c_dim1;
                    d__1 = *alpha * rtemp;
                    c__[i__2].r = d__1, c__[i__2].i = 0.;
                } else {
                    i__2 = j + j * c_dim1;
                    i__3 = j + j * c_dim1;
                    d__1 = *alpha * rtemp + *beta * c__[i__3].r;
                    c__[i__2].r = d__1, c__[i__2].i = 0.;
                }
                i__2 = *n;
                for (i__ = j + 1; i__ <= i__2; ++i__) {
                    temp.r = 0., temp.i = 0.;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l) {
                        d_lmp_cnjg(&z__3, &a[l + i__ * a_dim1]);
                        i__4 = l + j * a_dim1;
                        z__2.r = z__3.r * a[i__4].r - z__3.i * a[i__4].i,
                        z__2.i = z__3.r * a[i__4].i + z__3.i * a[i__4].r;
                        z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                        temp.r = z__1.r, temp.i = z__1.i;
                    }
                    if (*beta == 0.) {
                        i__3 = i__ + j * c_dim1;
                        z__1.r = *alpha * temp.r, z__1.i = *alpha * temp.i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    } else {
                        i__3 = i__ + j * c_dim1;
                        z__2.r = *alpha * temp.r, z__2.i = *alpha * temp.i;
                        i__4 = i__ + j * c_dim1;
                        z__3.r = *beta * c__[i__4].r, z__3.i = *beta * c__[i__4].i;
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
