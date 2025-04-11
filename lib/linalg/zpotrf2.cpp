#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublecomplex c_b1 = {1., 0.};
static doublereal c_b11 = -1.;
static doublereal c_b12 = 1.;
int zpotrf2_(char *uplo, integer *n, doublecomplex *a, integer *lda, integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, i__1;
    doublereal d__1;
    double sqrt(doublereal);
    integer n1, n2;
    doublereal ajj;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer iinfo;
    extern int zherk_(char *, char *, integer *, integer *, doublereal *, doublecomplex *,
                      integer *, doublereal *, doublecomplex *, integer *, ftnlen, ftnlen);
    logical upper;
    extern int ztrsm_(char *, char *, char *, char *, integer *, integer *, doublecomplex *,
                      doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, ftnlen,
                      ftnlen, ftnlen);
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
        xerbla_((char *)"ZPOTRF2", &i__1, (ftnlen)7);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    if (*n == 1) {
        i__1 = a_dim1 + 1;
        ajj = a[i__1].r;
        if (ajj <= 0. || disnan_(&ajj)) {
            *info = 1;
            return 0;
        }
        i__1 = a_dim1 + 1;
        d__1 = sqrt(ajj);
        a[i__1].r = d__1, a[i__1].i = 0.;
    } else {
        n1 = *n / 2;
        n2 = *n - n1;
        zpotrf2_(uplo, &n1, &a[a_dim1 + 1], lda, &iinfo, (ftnlen)1);
        if (iinfo != 0) {
            *info = iinfo;
            return 0;
        }
        if (upper) {
            ztrsm_((char *)"L", (char *)"U", (char *)"C", (char *)"N", &n1, &n2, &c_b1, &a[a_dim1 + 1], lda,
                   &a[(n1 + 1) * a_dim1 + 1], lda, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            zherk_(uplo, (char *)"C", &n2, &n1, &c_b11, &a[(n1 + 1) * a_dim1 + 1], lda, &c_b12,
                   &a[n1 + 1 + (n1 + 1) * a_dim1], lda, (ftnlen)1, (ftnlen)1);
            zpotrf2_(uplo, &n2, &a[n1 + 1 + (n1 + 1) * a_dim1], lda, &iinfo, (ftnlen)1);
            if (iinfo != 0) {
                *info = iinfo + n1;
                return 0;
            }
        } else {
            ztrsm_((char *)"R", (char *)"L", (char *)"C", (char *)"N", &n2, &n1, &c_b1, &a[a_dim1 + 1], lda, &a[n1 + 1 + a_dim1],
                   lda, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            zherk_(uplo, (char *)"N", &n2, &n1, &c_b11, &a[n1 + 1 + a_dim1], lda, &c_b12,
                   &a[n1 + 1 + (n1 + 1) * a_dim1], lda, (ftnlen)1, (ftnlen)1);
            zpotrf2_(uplo, &n2, &a[n1 + 1 + (n1 + 1) * a_dim1], lda, &iinfo, (ftnlen)1);
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
