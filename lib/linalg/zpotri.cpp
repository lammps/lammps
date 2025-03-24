#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int zpotri_(char *uplo, integer *n, doublecomplex *a, integer *lda, integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, i__1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int xerbla_(char *, integer *, ftnlen),
        zlauum_(char *, integer *, doublecomplex *, integer *, integer *, ftnlen),
        ztrtri_(char *, char *, integer *, doublecomplex *, integer *, integer *, ftnlen, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    *info = 0;
    if (!lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1) && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*lda < max(1, *n)) {
        *info = -4;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"ZPOTRI", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    ztrtri_(uplo, (char *)"N", n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)1);
    if (*info > 0) {
        return 0;
    }
    zlauum_(uplo, n, &a[a_offset], lda, info, (ftnlen)1);
    return 0;
}
#ifdef __cplusplus
}
#endif
