#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dsyr2k_(char *uplo, char *trans, integer *n, integer *k, doublereal *alpha, doublereal *a,
            integer *lda, doublereal *b, integer *ldb, doublereal *beta, doublereal *c__,
            integer *ldc, ftnlen uplo_len, ftnlen trans_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, i__3;
    integer i__, j, l, info;
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
               !lsame_(trans, (char *)"T", (ftnlen)1, (ftnlen)1) &&
               !lsame_(trans, (char *)"C", (ftnlen)1, (ftnlen)1)) {
        info = 2;
    } else if (*n < 0) {
        info = 3;
    } else if (*k < 0) {
        info = 4;
    } else if (*lda < max(1, nrowa)) {
        info = 7;
    } else if (*ldb < max(1, nrowa)) {
        info = 9;
    } else if (*ldc < max(1, *n)) {
        info = 12;
    }
    if (info != 0) {
        xerbla_((char *)"DSYR2K", &info, (ftnlen)6);
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
                        c__[i__ + j * c_dim1] = 0.;
                    }
                }
            } else {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                    }
                }
            }
        } else {
            if (*beta == 0.) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = 0.;
                    }
                }
            } else {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
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
                        c__[i__ + j * c_dim1] = 0.;
                    }
                } else if (*beta != 1.) {
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                    }
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l) {
                    if (a[j + l * a_dim1] != 0. || b[j + l * b_dim1] != 0.) {
                        temp1 = *alpha * b[j + l * b_dim1];
                        temp2 = *alpha * a[j + l * a_dim1];
                        i__3 = j;
                        for (i__ = 1; i__ <= i__3; ++i__) {
                            c__[i__ + j * c_dim1] = c__[i__ + j * c_dim1] +
                                                    a[i__ + l * a_dim1] * temp1 +
                                                    b[i__ + l * b_dim1] * temp2;
                        }
                    }
                }
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (*beta == 0.) {
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = 0.;
                    }
                } else if (*beta != 1.) {
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                    }
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l) {
                    if (a[j + l * a_dim1] != 0. || b[j + l * b_dim1] != 0.) {
                        temp1 = *alpha * b[j + l * b_dim1];
                        temp2 = *alpha * a[j + l * a_dim1];
                        i__3 = *n;
                        for (i__ = j; i__ <= i__3; ++i__) {
                            c__[i__ + j * c_dim1] = c__[i__ + j * c_dim1] +
                                                    a[i__ + l * a_dim1] * temp1 +
                                                    b[i__ + l * b_dim1] * temp2;
                        }
                    }
                }
            }
        }
    } else {
        if (upper) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = j;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp1 = 0.;
                    temp2 = 0.;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l) {
                        temp1 += a[l + i__ * a_dim1] * b[l + j * b_dim1];
                        temp2 += b[l + i__ * b_dim1] * a[l + j * a_dim1];
                    }
                    if (*beta == 0.) {
                        c__[i__ + j * c_dim1] = *alpha * temp1 + *alpha * temp2;
                    } else {
                        c__[i__ + j * c_dim1] =
                            *beta * c__[i__ + j * c_dim1] + *alpha * temp1 + *alpha * temp2;
                    }
                }
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *n;
                for (i__ = j; i__ <= i__2; ++i__) {
                    temp1 = 0.;
                    temp2 = 0.;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l) {
                        temp1 += a[l + i__ * a_dim1] * b[l + j * b_dim1];
                        temp2 += b[l + i__ * b_dim1] * a[l + j * a_dim1];
                    }
                    if (*beta == 0.) {
                        c__[i__ + j * c_dim1] = *alpha * temp1 + *alpha * temp2;
                    } else {
                        c__[i__ + j * c_dim1] =
                            *beta * c__[i__ + j * c_dim1] + *alpha * temp1 + *alpha * temp2;
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
