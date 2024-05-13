#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b9 = 1.;
static doublereal c_b11 = -1.;
int dpotrf2_(char *uplo, integer *n, doublereal *a, integer *lda, integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, i__1;
    double sqrt(doublereal);
    integer n1, n2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer iinfo;
    extern int dtrsm_(char *, char *, char *, char *, integer *, integer *, doublereal *,
                      doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen,
                      ftnlen);
    logical upper;
    extern int dsyrk_(char *, char *, integer *, integer *, doublereal *, doublereal *, integer *,
                      doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical disnan_(doublereal *);
    extern int xerbla_(char *, integer *, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    *info = 0;
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    if (!upper && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*lda < max(1, *n)) {
        *info = -4;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DPOTRF2", &i__1, (ftnlen)7);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    if (*n == 1) {
        if (a[a_dim1 + 1] <= 0. || disnan_(&a[a_dim1 + 1])) {
            *info = 1;
            return 0;
        }
        a[a_dim1 + 1] = sqrt(a[a_dim1 + 1]);
    } else {
        n1 = *n / 2;
        n2 = *n - n1;
        dpotrf2_(uplo, &n1, &a[a_dim1 + 1], lda, &iinfo, (ftnlen)1);
        if (iinfo != 0) {
            *info = iinfo;
            return 0;
        }
        if (upper) {
            dtrsm_((char *)"L", (char *)"U", (char *)"T", (char *)"N", &n1, &n2, &c_b9, &a[a_dim1 + 1], lda,
                   &a[(n1 + 1) * a_dim1 + 1], lda, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            dsyrk_(uplo, (char *)"T", &n2, &n1, &c_b11, &a[(n1 + 1) * a_dim1 + 1], lda, &c_b9,
                   &a[n1 + 1 + (n1 + 1) * a_dim1], lda, (ftnlen)1, (ftnlen)1);
            dpotrf2_(uplo, &n2, &a[n1 + 1 + (n1 + 1) * a_dim1], lda, &iinfo, (ftnlen)1);
            if (iinfo != 0) {
                *info = iinfo + n1;
                return 0;
            }
        } else {
            dtrsm_((char *)"R", (char *)"L", (char *)"T", (char *)"N", &n2, &n1, &c_b9, &a[a_dim1 + 1], lda, &a[n1 + 1 + a_dim1],
                   lda, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            dsyrk_(uplo, (char *)"N", &n2, &n1, &c_b11, &a[n1 + 1 + a_dim1], lda, &c_b9,
                   &a[n1 + 1 + (n1 + 1) * a_dim1], lda, (ftnlen)1, (ftnlen)1);
            dpotrf2_(uplo, &n2, &a[n1 + 1 + (n1 + 1) * a_dim1], lda, &iinfo, (ftnlen)1);
            if (iinfo != 0) {
                *info = iinfo + n1;
                return 0;
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
