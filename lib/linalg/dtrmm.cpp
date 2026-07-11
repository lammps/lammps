#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dtrmm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n,
           doublereal *alpha, doublereal *a, integer *lda, doublereal *b, integer *ldb,
           ftnlen side_len, ftnlen uplo_len, ftnlen transa_len, ftnlen diag_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    integer i__, j, k, info;
    doublereal temp;
    logical lside;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer nrowa;
    logical upper;
    extern int xerbla_(char *, integer *, ftnlen);
    logical nounit;
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
        xerbla_((char *)"DTRMM ", &info, (ftnlen)6);
        return 0;
    }
    if (*m == 0 || *n == 0) {
        return 0;
    }
    if (*alpha == 0.) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__2 = *m;
            for (i__ = 1; i__ <= i__2; ++i__) {
                b[i__ + j * b_dim1] = 0.;
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
                        if (b[k + j * b_dim1] != 0.) {
                            temp = *alpha * b[k + j * b_dim1];
                            i__3 = k - 1;
                            for (i__ = 1; i__ <= i__3; ++i__) {
                                b[i__ + j * b_dim1] += temp * a[i__ + k * a_dim1];
                            }
                            if (nounit) {
                                temp *= a[k + k * a_dim1];
                            }
                            b[k + j * b_dim1] = temp;
                        }
                    }
                }
            } else {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    for (k = *m; k >= 1; --k) {
                        if (b[k + j * b_dim1] != 0.) {
                            temp = *alpha * b[k + j * b_dim1];
                            b[k + j * b_dim1] = temp;
                            if (nounit) {
                                b[k + j * b_dim1] *= a[k + k * a_dim1];
                            }
                            i__2 = *m;
                            for (i__ = k + 1; i__ <= i__2; ++i__) {
                                b[i__ + j * b_dim1] += temp * a[i__ + k * a_dim1];
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
                        temp = b[i__ + j * b_dim1];
                        if (nounit) {
                            temp *= a[i__ + i__ * a_dim1];
                        }
                        i__2 = i__ - 1;
                        for (k = 1; k <= i__2; ++k) {
                            temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
                        }
                        b[i__ + j * b_dim1] = *alpha * temp;
                    }
                }
            } else {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        temp = b[i__ + j * b_dim1];
                        if (nounit) {
                            temp *= a[i__ + i__ * a_dim1];
                        }
                        i__3 = *m;
                        for (k = i__ + 1; k <= i__3; ++k) {
                            temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
                        }
                        b[i__ + j * b_dim1] = *alpha * temp;
                    }
                }
            }
        }
    } else {
        if (lsame_(transa, (char *)"N", (ftnlen)1, (ftnlen)1)) {
            if (upper) {
                for (j = *n; j >= 1; --j) {
                    temp = *alpha;
                    if (nounit) {
                        temp *= a[j + j * a_dim1];
                    }
                    i__1 = *m;
                    for (i__ = 1; i__ <= i__1; ++i__) {
                        b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
                    }
                    i__1 = j - 1;
                    for (k = 1; k <= i__1; ++k) {
                        if (a[k + j * a_dim1] != 0.) {
                            temp = *alpha * a[k + j * a_dim1];
                            i__2 = *m;
                            for (i__ = 1; i__ <= i__2; ++i__) {
                                b[i__ + j * b_dim1] += temp * b[i__ + k * b_dim1];
                            }
                        }
                    }
                }
            } else {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    temp = *alpha;
                    if (nounit) {
                        temp *= a[j + j * a_dim1];
                    }
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
                    }
                    i__2 = *n;
                    for (k = j + 1; k <= i__2; ++k) {
                        if (a[k + j * a_dim1] != 0.) {
                            temp = *alpha * a[k + j * a_dim1];
                            i__3 = *m;
                            for (i__ = 1; i__ <= i__3; ++i__) {
                                b[i__ + j * b_dim1] += temp * b[i__ + k * b_dim1];
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
                        if (a[j + k * a_dim1] != 0.) {
                            temp = *alpha * a[j + k * a_dim1];
                            i__3 = *m;
                            for (i__ = 1; i__ <= i__3; ++i__) {
                                b[i__ + j * b_dim1] += temp * b[i__ + k * b_dim1];
                            }
                        }
                    }
                    temp = *alpha;
                    if (nounit) {
                        temp *= a[k + k * a_dim1];
                    }
                    if (temp != 1.) {
                        i__2 = *m;
                        for (i__ = 1; i__ <= i__2; ++i__) {
                            b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
                        }
                    }
                }
            } else {
                for (k = *n; k >= 1; --k) {
                    i__1 = *n;
                    for (j = k + 1; j <= i__1; ++j) {
                        if (a[j + k * a_dim1] != 0.) {
                            temp = *alpha * a[j + k * a_dim1];
                            i__2 = *m;
                            for (i__ = 1; i__ <= i__2; ++i__) {
                                b[i__ + j * b_dim1] += temp * b[i__ + k * b_dim1];
                            }
                        }
                    }
                    temp = *alpha;
                    if (nounit) {
                        temp *= a[k + k * a_dim1];
                    }
                    if (temp != 1.) {
                        i__1 = *m;
                        for (i__ = 1; i__ <= i__1; ++i__) {
                            b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
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
