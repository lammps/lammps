#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dgemm_(char *transa, char *transb, integer *m, integer *n, integer *k, doublereal *alpha,
           doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *beta,
           doublereal *c__, integer *ldc, ftnlen transa_len, ftnlen transb_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, i__3;
    integer i__, j, l, info;
    logical nota, notb;
    doublereal temp;
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
    if (!nota && !lsame_(transa, (char *)"C", (ftnlen)1, (ftnlen)1) &&
        !lsame_(transa, (char *)"T", (ftnlen)1, (ftnlen)1)) {
        info = 1;
    } else if (!notb && !lsame_(transb, (char *)"C", (ftnlen)1, (ftnlen)1) &&
               !lsame_(transb, (char *)"T", (ftnlen)1, (ftnlen)1)) {
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
        xerbla_((char *)"DGEMM ", &info, (ftnlen)6);
        return 0;
    }
    if (*m == 0 || *n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
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
    if (notb) {
        if (nota) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (*beta == 0.) {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = 0.;
                    }
                } else if (*beta != 1.) {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                    }
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l) {
                    temp = *alpha * b[l + j * b_dim1];
                    i__3 = *m;
                    for (i__ = 1; i__ <= i__3; ++i__) {
                        c__[i__ + j * c_dim1] += temp * a[i__ + l * a_dim1];
                    }
                }
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp = 0.;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l) {
                        temp += a[l + i__ * a_dim1] * b[l + j * b_dim1];
                    }
                    if (*beta == 0.) {
                        c__[i__ + j * c_dim1] = *alpha * temp;
                    } else {
                        c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[i__ + j * c_dim1];
                    }
                }
            }
        }
    } else {
        if (nota) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (*beta == 0.) {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = 0.;
                    }
                } else if (*beta != 1.) {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
                    }
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l) {
                    temp = *alpha * b[j + l * b_dim1];
                    i__3 = *m;
                    for (i__ = 1; i__ <= i__3; ++i__) {
                        c__[i__ + j * c_dim1] += temp * a[i__ + l * a_dim1];
                    }
                }
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp = 0.;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l) {
                        temp += a[l + i__ * a_dim1] * b[j + l * b_dim1];
                    }
                    if (*beta == 0.) {
                        c__[i__ + j * c_dim1] = *alpha * temp;
                    } else {
                        c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[i__ + j * c_dim1];
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
