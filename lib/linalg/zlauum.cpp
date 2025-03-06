#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublecomplex c_b1 = {1., 0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b21 = 1.;
int zlauum_(char *uplo, integer *n, doublecomplex *a, integer *lda, integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    integer i__, ib, nb;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int zgemm_(char *, char *, integer *, integer *, integer *, doublecomplex *,
                      doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *,
                      doublecomplex *, integer *, ftnlen, ftnlen),
        zherk_(char *, char *, integer *, integer *, doublereal *, doublecomplex *, integer *,
               doublereal *, doublecomplex *, integer *, ftnlen, ftnlen);
    logical upper;
    extern int ztrmm_(char *, char *, char *, char *, integer *, integer *, doublecomplex *,
                      doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, ftnlen,
                      ftnlen, ftnlen),
        zlauu2_(char *, integer *, doublecomplex *, integer *, integer *, ftnlen),
        xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
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
        xerbla_((char *)"ZLAUUM", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    nb = ilaenv_(&c__1, (char *)"ZLAUUM", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
    if (nb <= 1 || nb >= *n) {
        zlauu2_(uplo, n, &a[a_offset], lda, info, (ftnlen)1);
    } else {
        if (upper) {
            i__1 = *n;
            i__2 = nb;
            for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
                i__3 = nb, i__4 = *n - i__ + 1;
                ib = min(i__3, i__4);
                i__3 = i__ - 1;
                ztrmm_((char *)"R", (char *)"U", (char *)"C", (char *)"N", &i__3, &ib, &c_b1, &a[i__ + i__ * a_dim1], lda,
                       &a[i__ * a_dim1 + 1], lda, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                zlauu2_((char *)"U", &ib, &a[i__ + i__ * a_dim1], lda, info, (ftnlen)1);
                if (i__ + ib <= *n) {
                    i__3 = i__ - 1;
                    i__4 = *n - i__ - ib + 1;
                    zgemm_((char *)"N", (char *)"C", &i__3, &ib, &i__4, &c_b1, &a[(i__ + ib) * a_dim1 + 1], lda,
                           &a[i__ + (i__ + ib) * a_dim1], lda, &c_b1, &a[i__ * a_dim1 + 1], lda,
                           (ftnlen)1, (ftnlen)1);
                    i__3 = *n - i__ - ib + 1;
                    zherk_((char *)"U", (char *)"N", &ib, &i__3, &c_b21, &a[i__ + (i__ + ib) * a_dim1], lda, &c_b21,
                           &a[i__ + i__ * a_dim1], lda, (ftnlen)1, (ftnlen)1);
                }
            }
        } else {
            i__2 = *n;
            i__1 = nb;
            for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
                i__3 = nb, i__4 = *n - i__ + 1;
                ib = min(i__3, i__4);
                i__3 = i__ - 1;
                ztrmm_((char *)"L", (char *)"L", (char *)"C", (char *)"N", &ib, &i__3, &c_b1, &a[i__ + i__ * a_dim1], lda,
                       &a[i__ + a_dim1], lda, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                zlauu2_((char *)"L", &ib, &a[i__ + i__ * a_dim1], lda, info, (ftnlen)1);
                if (i__ + ib <= *n) {
                    i__3 = i__ - 1;
                    i__4 = *n - i__ - ib + 1;
                    zgemm_((char *)"C", (char *)"N", &ib, &i__3, &i__4, &c_b1, &a[i__ + ib + i__ * a_dim1], lda,
                           &a[i__ + ib + a_dim1], lda, &c_b1, &a[i__ + a_dim1], lda, (ftnlen)1,
                           (ftnlen)1);
                    i__3 = *n - i__ - ib + 1;
                    zherk_((char *)"L", (char *)"C", &ib, &i__3, &c_b21, &a[i__ + ib + i__ * a_dim1], lda, &c_b21,
                           &a[i__ + i__ * a_dim1], lda, (ftnlen)1, (ftnlen)1);
                }
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
