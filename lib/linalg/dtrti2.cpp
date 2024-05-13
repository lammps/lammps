#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
int dtrti2_(char *uplo, char *diag, integer *n, doublereal *a, integer *lda, integer *info,
            ftnlen uplo_len, ftnlen diag_len)
{
    integer a_dim1, a_offset, i__1, i__2;
    integer j;
    doublereal ajj;
    extern int dscal_(integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    logical upper;
    extern int dtrmv_(char *, char *, char *, integer *, doublereal *, integer *, doublereal *,
                      integer *, ftnlen, ftnlen, ftnlen),
        xerbla_(char *, integer *, ftnlen);
    logical nounit;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    *info = 0;
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    nounit = lsame_(diag, (char *)"N", (ftnlen)1, (ftnlen)1);
    if (!upper && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (!nounit && !lsame_(diag, (char *)"U", (ftnlen)1, (ftnlen)1)) {
        *info = -2;
    } else if (*n < 0) {
        *info = -3;
    } else if (*lda < max(1, *n)) {
        *info = -5;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DTRTI2", &i__1, (ftnlen)6);
        return 0;
    }
    if (upper) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            if (nounit) {
                a[j + j * a_dim1] = 1. / a[j + j * a_dim1];
                ajj = -a[j + j * a_dim1];
            } else {
                ajj = -1.;
            }
            i__2 = j - 1;
            dtrmv_((char *)"Upper", (char *)"No transpose", diag, &i__2, &a[a_offset], lda, &a[j * a_dim1 + 1],
                   &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)1);
            i__2 = j - 1;
            dscal_(&i__2, &ajj, &a[j * a_dim1 + 1], &c__1);
        }
    } else {
        for (j = *n; j >= 1; --j) {
            if (nounit) {
                a[j + j * a_dim1] = 1. / a[j + j * a_dim1];
                ajj = -a[j + j * a_dim1];
            } else {
                ajj = -1.;
            }
            if (j < *n) {
                i__1 = *n - j;
                dtrmv_((char *)"Lower", (char *)"No transpose", diag, &i__1, &a[j + 1 + (j + 1) * a_dim1], lda,
                       &a[j + 1 + j * a_dim1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)1);
                i__1 = *n - j;
                dscal_(&i__1, &ajj, &a[j + 1 + j * a_dim1], &c__1);
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
