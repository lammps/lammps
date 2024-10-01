#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dpotri_(char *uplo, integer *n, doublereal *a, integer *lda, integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, i__1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int xerbla_(char *, integer *, ftnlen),
        dlauum_(char *, integer *, doublereal *, integer *, integer *, ftnlen),
        dtrtri_(char *, char *, integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
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
        xerbla_((char *)"DPOTRI", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    dtrtri_(uplo, (char *)"Non-unit", n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)8);
    if (*info > 0) {
        return 0;
    }
    dlauum_(uplo, n, &a[a_offset], lda, info, (ftnlen)1);
    return 0;
}
#ifdef __cplusplus
}
#endif
