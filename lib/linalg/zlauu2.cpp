#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublecomplex c_b1 = {1., 0.};
static integer c__1 = 1;
int zlauu2_(char *uplo, integer *n, doublecomplex *a, integer *lda, integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1;
    integer i__;
    doublereal aii;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern VOID zdotc_(doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *,
                       integer *);
    extern int zgemv_(char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *,
                      doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *,
                      ftnlen);
    logical upper;
    extern int xerbla_(char *, integer *, ftnlen),
        zdscal_(integer *, doublereal *, doublecomplex *, integer *),
        zlacgv_(integer *, doublecomplex *, integer *);
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
        xerbla_((char *)"ZLAUU2", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    if (upper) {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            i__2 = i__ + i__ * a_dim1;
            aii = a[i__2].r;
            if (i__ < *n) {
                i__2 = i__ + i__ * a_dim1;
                i__3 = *n - i__;
                zdotc_(&z__1, &i__3, &a[i__ + (i__ + 1) * a_dim1], lda,
                       &a[i__ + (i__ + 1) * a_dim1], lda);
                d__1 = aii * aii + z__1.r;
                a[i__2].r = d__1, a[i__2].i = 0.;
                i__2 = *n - i__;
                zlacgv_(&i__2, &a[i__ + (i__ + 1) * a_dim1], lda);
                i__2 = i__ - 1;
                i__3 = *n - i__;
                z__1.r = aii, z__1.i = 0.;
                zgemv_((char *)"N", &i__2, &i__3, &c_b1, &a[(i__ + 1) * a_dim1 + 1], lda,
                       &a[i__ + (i__ + 1) * a_dim1], lda, &z__1, &a[i__ * a_dim1 + 1], &c__1,
                       (ftnlen)1);
                i__2 = *n - i__;
                zlacgv_(&i__2, &a[i__ + (i__ + 1) * a_dim1], lda);
            } else {
                zdscal_(&i__, &aii, &a[i__ * a_dim1 + 1], &c__1);
            }
        }
    } else {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            i__2 = i__ + i__ * a_dim1;
            aii = a[i__2].r;
            if (i__ < *n) {
                i__2 = i__ + i__ * a_dim1;
                i__3 = *n - i__;
                zdotc_(&z__1, &i__3, &a[i__ + 1 + i__ * a_dim1], &c__1, &a[i__ + 1 + i__ * a_dim1],
                       &c__1);
                d__1 = aii * aii + z__1.r;
                a[i__2].r = d__1, a[i__2].i = 0.;
                i__2 = i__ - 1;
                zlacgv_(&i__2, &a[i__ + a_dim1], lda);
                i__2 = *n - i__;
                i__3 = i__ - 1;
                z__1.r = aii, z__1.i = 0.;
                zgemv_((char *)"C", &i__2, &i__3, &c_b1, &a[i__ + 1 + a_dim1], lda,
                       &a[i__ + 1 + i__ * a_dim1], &c__1, &z__1, &a[i__ + a_dim1], lda, (ftnlen)1);
                i__2 = i__ - 1;
                zlacgv_(&i__2, &a[i__ + a_dim1], lda);
            } else {
                zdscal_(&i__, &aii, &a[i__ + a_dim1], lda);
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
