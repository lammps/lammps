#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static doublereal c_b12 = 1.;
static integer c_n1 = -1;
int dgetrs_(char *trans, integer *n, integer *nrhs, doublereal *a, integer *lda, integer *ipiv,
            doublereal *b, integer *ldb, integer *info, ftnlen trans_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int dtrsm_(char *, char *, char *, char *, integer *, integer *, doublereal *,
                      doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen,
                      ftnlen),
        xerbla_(char *, integer *, ftnlen),
        dlaswp_(integer *, doublereal *, integer *, integer *, integer *, integer *, integer *);
    logical notran;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    *info = 0;
    notran = lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1);
    if (!notran && !lsame_(trans, (char *)"T", (ftnlen)1, (ftnlen)1) &&
        !lsame_(trans, (char *)"C", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*nrhs < 0) {
        *info = -3;
    } else if (*lda < max(1, *n)) {
        *info = -5;
    } else if (*ldb < max(1, *n)) {
        *info = -8;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DGETRS", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0 || *nrhs == 0) {
        return 0;
    }
    if (notran) {
        dlaswp_(nrhs, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c__1);
        dtrsm_((char *)"Left", (char *)"Lower", (char *)"No transpose", (char *)"Unit", n, nrhs, &c_b12, &a[a_offset], lda,
               &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)4);
        dtrsm_((char *)"Left", (char *)"Upper", (char *)"No transpose", (char *)"Non-unit", n, nrhs, &c_b12, &a[a_offset], lda,
               &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)8);
    } else {
        dtrsm_((char *)"Left", (char *)"Upper", (char *)"Transpose", (char *)"Non-unit", n, nrhs, &c_b12, &a[a_offset], lda,
               &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (ftnlen)9, (ftnlen)8);
        dtrsm_((char *)"Left", (char *)"Lower", (char *)"Transpose", (char *)"Unit", n, nrhs, &c_b12, &a[a_offset], lda,
               &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (ftnlen)9, (ftnlen)4);
        dlaswp_(nrhs, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c_n1);
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
