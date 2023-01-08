#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int zlacpy_(char *uplo, integer *m, integer *n, doublecomplex *a, integer *lda, doublecomplex *b,
            integer *ldb, ftnlen uplo_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;
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
                i__3 = i__ + j * b_dim1;
                i__4 = i__ + j * a_dim1;
                b[i__3].r = a[i__4].r, b[i__3].i = a[i__4].i;
            }
        }
    } else if (lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__2 = *m;
            for (i__ = j; i__ <= i__2; ++i__) {
                i__3 = i__ + j * b_dim1;
                i__4 = i__ + j * a_dim1;
                b[i__3].r = a[i__4].r, b[i__3].i = a[i__4].i;
            }
        }
    } else {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__2 = *m;
            for (i__ = 1; i__ <= i__2; ++i__) {
                i__3 = i__ + j * b_dim1;
                i__4 = i__ + j * a_dim1;
                b[i__3].r = a[i__4].r, b[i__3].i = a[i__4].i;
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
