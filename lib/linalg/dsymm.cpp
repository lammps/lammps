#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dsymm_(char *side, char *uplo, integer *m, integer *n, doublereal *alpha, doublereal *a,
           integer *lda, doublereal *b, integer *ldb, doublereal *beta, doublereal *c__,
           integer *ldc, ftnlen side_len, ftnlen uplo_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, i__3;
    integer i__, j, k, info;
    doublereal temp1, temp2;
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
        xerbla_((char *)"DSYMM ", &info, (ftnlen)6);
        return 0;
    }
    if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
        return 0;
    }
    if (*alpha == 0.) {
        if (*beta == 0.) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    c__[i__ + j * c_dim1] = 0.;
                }
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
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
                    temp1 = *alpha * b[i__ + j * b_dim1];
                    temp2 = 0.;
                    i__3 = i__ - 1;
                    for (k = 1; k <= i__3; ++k) {
                        c__[k + j * c_dim1] += temp1 * a[k + i__ * a_dim1];
                        temp2 += b[k + j * b_dim1] * a[k + i__ * a_dim1];
                    }
                    if (*beta == 0.) {
                        c__[i__ + j * c_dim1] = temp1 * a[i__ + i__ * a_dim1] + *alpha * temp2;
                    } else {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] +
                                                temp1 * a[i__ + i__ * a_dim1] + *alpha * temp2;
                    }
                }
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                for (i__ = *m; i__ >= 1; --i__) {
                    temp1 = *alpha * b[i__ + j * b_dim1];
                    temp2 = 0.;
                    i__2 = *m;
                    for (k = i__ + 1; k <= i__2; ++k) {
                        c__[k + j * c_dim1] += temp1 * a[k + i__ * a_dim1];
                        temp2 += b[k + j * b_dim1] * a[k + i__ * a_dim1];
                    }
                    if (*beta == 0.) {
                        c__[i__ + j * c_dim1] = temp1 * a[i__ + i__ * a_dim1] + *alpha * temp2;
                    } else {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] +
                                                temp1 * a[i__ + i__ * a_dim1] + *alpha * temp2;
                    }
                }
            }
        }
    } else {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            temp1 = *alpha * a[j + j * a_dim1];
            if (*beta == 0.) {
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    c__[i__ + j * c_dim1] = temp1 * b[i__ + j * b_dim1];
                }
            } else {
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    c__[i__ + j * c_dim1] =
                        *beta * c__[i__ + j * c_dim1] + temp1 * b[i__ + j * b_dim1];
                }
            }
            i__2 = j - 1;
            for (k = 1; k <= i__2; ++k) {
                if (upper) {
                    temp1 = *alpha * a[k + j * a_dim1];
                } else {
                    temp1 = *alpha * a[j + k * a_dim1];
                }
                i__3 = *m;
                for (i__ = 1; i__ <= i__3; ++i__) {
                    c__[i__ + j * c_dim1] += temp1 * b[i__ + k * b_dim1];
                }
            }
            i__2 = *n;
            for (k = j + 1; k <= i__2; ++k) {
                if (upper) {
                    temp1 = *alpha * a[j + k * a_dim1];
                } else {
                    temp1 = *alpha * a[k + j * a_dim1];
                }
                i__3 = *m;
                for (i__ = 1; i__ <= i__3; ++i__) {
                    c__[i__ + j * c_dim1] += temp1 * b[i__ + k * b_dim1];
                }
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
