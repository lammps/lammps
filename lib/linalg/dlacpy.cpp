#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlacpy_(char *uplo, integer *m, integer *n, doublereal *a, integer *lda, doublereal *b,
            integer *ldb, ftnlen uplo_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    integer i__, j;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    if (lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1)) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__2 = min(j, *m);
            for (i__ = 1; i__ <= i__2; ++i__) {
                b[i__ + j * b_dim1] = a[i__ + j * a_dim1];
            }
        }
    } else if (lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__2 = *m;
            for (i__ = j; i__ <= i__2; ++i__) {
                b[i__ + j * b_dim1] = a[i__ + j * a_dim1];
            }
        }
    } else {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__2 = *m;
            for (i__ = 1; i__ <= i__2; ++i__) {
                b[i__ + j * b_dim1] = a[i__ + j * a_dim1];
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
