#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublecomplex c_b1 = {1., 0.};
static integer c__1 = 1;
int ztrti2_(char *uplo, char *diag, integer *n, doublecomplex *a, integer *lda, integer *info,
            ftnlen uplo_len, ftnlen diag_len)
{
    integer a_dim1, a_offset, i__1, i__2;
    doublecomplex z__1;
    void z_lmp_div(doublecomplex *, doublecomplex *, doublecomplex *);
    integer j;
    doublecomplex ajj;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int zscal_(integer *, doublecomplex *, doublecomplex *, integer *);
    logical upper;
    extern int ztrmv_(char *, char *, char *, integer *, doublecomplex *, integer *,
                      doublecomplex *, integer *, ftnlen, ftnlen, ftnlen),
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
        xerbla_((char *)"ZTRTI2", &i__1, (ftnlen)6);
        return 0;
    }
    if (upper) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            if (nounit) {
                i__2 = j + j * a_dim1;
                z_lmp_div(&z__1, &c_b1, &a[j + j * a_dim1]);
                a[i__2].r = z__1.r, a[i__2].i = z__1.i;
                i__2 = j + j * a_dim1;
                z__1.r = -a[i__2].r, z__1.i = -a[i__2].i;
                ajj.r = z__1.r, ajj.i = z__1.i;
            } else {
                z__1.r = -1., z__1.i = -0.;
                ajj.r = z__1.r, ajj.i = z__1.i;
            }
            i__2 = j - 1;
            ztrmv_((char *)"U", (char *)"N", diag, &i__2, &a[a_offset], lda, &a[j * a_dim1 + 1], &c__1, (ftnlen)1,
                   (ftnlen)1, (ftnlen)1);
            i__2 = j - 1;
            zscal_(&i__2, &ajj, &a[j * a_dim1 + 1], &c__1);
        }
    } else {
        for (j = *n; j >= 1; --j) {
            if (nounit) {
                i__1 = j + j * a_dim1;
                z_lmp_div(&z__1, &c_b1, &a[j + j * a_dim1]);
                a[i__1].r = z__1.r, a[i__1].i = z__1.i;
                i__1 = j + j * a_dim1;
                z__1.r = -a[i__1].r, z__1.i = -a[i__1].i;
                ajj.r = z__1.r, ajj.i = z__1.i;
            } else {
                z__1.r = -1., z__1.i = -0.;
                ajj.r = z__1.r, ajj.i = z__1.i;
            }
            if (j < *n) {
                i__1 = *n - j;
                ztrmv_((char *)"L", (char *)"N", diag, &i__1, &a[j + 1 + (j + 1) * a_dim1], lda,
                       &a[j + 1 + j * a_dim1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                i__1 = *n - j;
                zscal_(&i__1, &ajj, &a[j + 1 + j * a_dim1], &c__1);
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
