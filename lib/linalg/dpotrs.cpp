#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b9 = 1.;
int dpotrs_(char *uplo, integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *b,
            integer *ldb, integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int dtrsm_(char *, char *, char *, char *, integer *, integer *, doublereal *,
                      doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen,
                      ftnlen);
    logical upper;
    extern int xerbla_(char *, integer *, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    *info = 0;
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    if (!upper && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*nrhs < 0) {
        *info = -3;
    } else if (*lda < max(1, *n)) {
        *info = -5;
    } else if (*ldb < max(1, *n)) {
        *info = -7;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DPOTRS", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0 || *nrhs == 0) {
        return 0;
    }
    if (upper) {
        dtrsm_((char *)"Left", (char *)"Upper", (char *)"Transpose", (char *)"Non-unit", n, nrhs, &c_b9, &a[a_offset], lda,
               &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (ftnlen)9, (ftnlen)8);
        dtrsm_((char *)"Left", (char *)"Upper", (char *)"No transpose", (char *)"Non-unit", n, nrhs, &c_b9, &a[a_offset], lda,
               &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)8);
    } else {
        dtrsm_((char *)"Left", (char *)"Lower", (char *)"No transpose", (char *)"Non-unit", n, nrhs, &c_b9, &a[a_offset], lda,
               &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)8);
        dtrsm_((char *)"Left", (char *)"Lower", (char *)"Transpose", (char *)"Non-unit", n, nrhs, &c_b9, &a[a_offset], lda,
               &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (ftnlen)9, (ftnlen)8);
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
