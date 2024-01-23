#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static doublereal c_b10 = -1.;
static doublereal c_b12 = 1.;
int dpotf2_(char *uplo, integer *n, doublereal *a, integer *lda, integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;
    double sqrt(doublereal);
    integer j;
    doublereal ajj;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, integer *);
    extern int dscal_(integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *,
                      doublereal *, integer *, doublereal *, doublereal *, integer *, ftnlen);
    logical upper;
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
        xerbla_((char *)"DPOTF2", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    if (upper) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__2 = j - 1;
            ajj = a[j + j * a_dim1] -
                  ddot_(&i__2, &a[j * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &c__1);
            if (ajj <= 0. || disnan_(&ajj)) {
                a[j + j * a_dim1] = ajj;
                goto L30;
            }
            ajj = sqrt(ajj);
            a[j + j * a_dim1] = ajj;
            if (j < *n) {
                i__2 = j - 1;
                i__3 = *n - j;
                dgemv_((char *)"Transpose", &i__2, &i__3, &c_b10, &a[(j + 1) * a_dim1 + 1], lda,
                       &a[j * a_dim1 + 1], &c__1, &c_b12, &a[j + (j + 1) * a_dim1], lda, (ftnlen)9);
                i__2 = *n - j;
                d__1 = 1. / ajj;
                dscal_(&i__2, &d__1, &a[j + (j + 1) * a_dim1], lda);
            }
        }
    } else {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__2 = j - 1;
            ajj = a[j + j * a_dim1] - ddot_(&i__2, &a[j + a_dim1], lda, &a[j + a_dim1], lda);
            if (ajj <= 0. || disnan_(&ajj)) {
                a[j + j * a_dim1] = ajj;
                goto L30;
            }
            ajj = sqrt(ajj);
            a[j + j * a_dim1] = ajj;
            if (j < *n) {
                i__2 = *n - j;
                i__3 = j - 1;
                dgemv_((char *)"No transpose", &i__2, &i__3, &c_b10, &a[j + 1 + a_dim1], lda,
                       &a[j + a_dim1], lda, &c_b12, &a[j + 1 + j * a_dim1], &c__1, (ftnlen)12);
                i__2 = *n - j;
                d__1 = 1. / ajj;
                dscal_(&i__2, &d__1, &a[j + 1 + j * a_dim1], &c__1);
            }
        }
    }
    goto L40;
L30:
    *info = j;
L40:
    return 0;
}
#ifdef __cplusplus
}
#endif
