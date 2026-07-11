#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b12 = 1.;
int dtrtrs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, doublereal *a,
            integer *lda, doublereal *b, integer *ldb, integer *info, ftnlen uplo_len,
            ftnlen trans_len, ftnlen diag_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int dtrsm_(char *, char *, char *, char *, integer *, integer *, doublereal *,
                      doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen,
                      ftnlen),
        xerbla_(char *, integer *, ftnlen);
    logical nounit;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    *info = 0;
    nounit = lsame_(diag, (char *)"N", (ftnlen)1, (ftnlen)1);
    if (!lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1) && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (!lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1) &&
               !lsame_(trans, (char *)"T", (ftnlen)1, (ftnlen)1) &&
               !lsame_(trans, (char *)"C", (ftnlen)1, (ftnlen)1)) {
        *info = -2;
    } else if (!nounit && !lsame_(diag, (char *)"U", (ftnlen)1, (ftnlen)1)) {
        *info = -3;
    } else if (*n < 0) {
        *info = -4;
    } else if (*nrhs < 0) {
        *info = -5;
    } else if (*lda < max(1, *n)) {
        *info = -7;
    } else if (*ldb < max(1, *n)) {
        *info = -9;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DTRTRS", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    if (nounit) {
        i__1 = *n;
        for (*info = 1; *info <= i__1; ++(*info)) {
            if (a[*info + *info * a_dim1] == 0.) {
                return 0;
            }
        }
    }
    *info = 0;
    dtrsm_((char *)"L", uplo, trans, diag, n, nrhs, &c_b12, &a[a_offset], lda, &b[b_offset], ldb, (ftnlen)1,
           (ftnlen)1, (ftnlen)1, (ftnlen)1);
    return 0;
}
#ifdef __cplusplus
}
#endif
