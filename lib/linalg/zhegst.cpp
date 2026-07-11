#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublecomplex c_b1 = {1., 0.};
static doublecomplex c_b2 = {.5, 0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b18 = 1.;
int zhegst_(integer *itype, char *uplo, integer *n, doublecomplex *a, integer *lda,
            doublecomplex *b, integer *ldb, integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    doublecomplex z__1;
    integer k, kb, nb;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int zhemm_(char *, char *, integer *, integer *, doublecomplex *, doublecomplex *,
                      integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *,
                      integer *, ftnlen, ftnlen);
    logical upper;
    extern int ztrmm_(char *, char *, char *, char *, integer *, integer *, doublecomplex *,
                      doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, ftnlen,
                      ftnlen, ftnlen),
        ztrsm_(char *, char *, char *, char *, integer *, integer *, doublecomplex *,
               doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, ftnlen, ftnlen,
               ftnlen),
        zhegs2_(integer *, char *, integer *, doublecomplex *, integer *, doublecomplex *,
                integer *, integer *, ftnlen),
        zher2k_(char *, char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *,
                doublecomplex *, integer *, doublereal *, doublecomplex *, integer *, ftnlen,
                ftnlen),
        xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    *info = 0;
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    if (*itype < 1 || *itype > 3) {
        *info = -1;
    } else if (!upper && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -2;
    } else if (*n < 0) {
        *info = -3;
    } else if (*lda < max(1, *n)) {
        *info = -5;
    } else if (*ldb < max(1, *n)) {
        *info = -7;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"ZHEGST", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    nb = ilaenv_(&c__1, (char *)"ZHEGST", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
    if (nb <= 1 || nb >= *n) {
        zhegs2_(itype, uplo, n, &a[a_offset], lda, &b[b_offset], ldb, info, (ftnlen)1);
    } else {
        if (*itype == 1) {
            if (upper) {
                i__1 = *n;
                i__2 = nb;
                for (k = 1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2) {
                    i__3 = *n - k + 1;
                    kb = min(i__3, nb);
                    zhegs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + k * b_dim1], ldb,
                            info, (ftnlen)1);
                    if (k + kb <= *n) {
                        i__3 = *n - k - kb + 1;
                        ztrsm_((char *)"L", uplo, (char *)"C", (char *)"N", &kb, &i__3, &c_b1, &b[k + k * b_dim1], ldb,
                               &a[k + (k + kb) * a_dim1], lda, (ftnlen)1, (ftnlen)1, (ftnlen)1,
                               (ftnlen)1);
                        i__3 = *n - k - kb + 1;
                        z__1.r = -.5, z__1.i = -0.;
                        zhemm_((char *)"L", uplo, &kb, &i__3, &z__1, &a[k + k * a_dim1], lda,
                               &b[k + (k + kb) * b_dim1], ldb, &c_b1, &a[k + (k + kb) * a_dim1],
                               lda, (ftnlen)1, (ftnlen)1);
                        i__3 = *n - k - kb + 1;
                        z__1.r = -1., z__1.i = -0.;
                        zher2k_(uplo, (char *)"C", &i__3, &kb, &z__1, &a[k + (k + kb) * a_dim1], lda,
                                &b[k + (k + kb) * b_dim1], ldb, &c_b18,
                                &a[k + kb + (k + kb) * a_dim1], lda, (ftnlen)1, (ftnlen)1);
                        i__3 = *n - k - kb + 1;
                        z__1.r = -.5, z__1.i = -0.;
                        zhemm_((char *)"L", uplo, &kb, &i__3, &z__1, &a[k + k * a_dim1], lda,
                               &b[k + (k + kb) * b_dim1], ldb, &c_b1, &a[k + (k + kb) * a_dim1],
                               lda, (ftnlen)1, (ftnlen)1);
                        i__3 = *n - k - kb + 1;
                        ztrsm_((char *)"R", uplo, (char *)"N", (char *)"N", &kb, &i__3, &c_b1,
                               &b[k + kb + (k + kb) * b_dim1], ldb, &a[k + (k + kb) * a_dim1], lda,
                               (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                    }
                }
            } else {
                i__2 = *n;
                i__1 = nb;
                for (k = 1; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1) {
                    i__3 = *n - k + 1;
                    kb = min(i__3, nb);
                    zhegs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + k * b_dim1], ldb,
                            info, (ftnlen)1);
                    if (k + kb <= *n) {
                        i__3 = *n - k - kb + 1;
                        ztrsm_((char *)"R", uplo, (char *)"C", (char *)"N", &i__3, &kb, &c_b1, &b[k + k * b_dim1], ldb,
                               &a[k + kb + k * a_dim1], lda, (ftnlen)1, (ftnlen)1, (ftnlen)1,
                               (ftnlen)1);
                        i__3 = *n - k - kb + 1;
                        z__1.r = -.5, z__1.i = -0.;
                        zhemm_((char *)"R", uplo, &i__3, &kb, &z__1, &a[k + k * a_dim1], lda,
                               &b[k + kb + k * b_dim1], ldb, &c_b1, &a[k + kb + k * a_dim1], lda,
                               (ftnlen)1, (ftnlen)1);
                        i__3 = *n - k - kb + 1;
                        z__1.r = -1., z__1.i = -0.;
                        zher2k_(uplo, (char *)"N", &i__3, &kb, &z__1, &a[k + kb + k * a_dim1], lda,
                                &b[k + kb + k * b_dim1], ldb, &c_b18,
                                &a[k + kb + (k + kb) * a_dim1], lda, (ftnlen)1, (ftnlen)1);
                        i__3 = *n - k - kb + 1;
                        z__1.r = -.5, z__1.i = -0.;
                        zhemm_((char *)"R", uplo, &i__3, &kb, &z__1, &a[k + k * a_dim1], lda,
                               &b[k + kb + k * b_dim1], ldb, &c_b1, &a[k + kb + k * a_dim1], lda,
                               (ftnlen)1, (ftnlen)1);
                        i__3 = *n - k - kb + 1;
                        ztrsm_((char *)"L", uplo, (char *)"N", (char *)"N", &i__3, &kb, &c_b1,
                               &b[k + kb + (k + kb) * b_dim1], ldb, &a[k + kb + k * a_dim1], lda,
                               (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                    }
                }
            }
        } else {
            if (upper) {
                i__1 = *n;
                i__2 = nb;
                for (k = 1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2) {
                    i__3 = *n - k + 1;
                    kb = min(i__3, nb);
                    i__3 = k - 1;
                    ztrmm_((char *)"L", uplo, (char *)"N", (char *)"N", &i__3, &kb, &c_b1, &b[b_offset], ldb,
                           &a[k * a_dim1 + 1], lda, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                    i__3 = k - 1;
                    zhemm_((char *)"R", uplo, &i__3, &kb, &c_b2, &a[k + k * a_dim1], lda,
                           &b[k * b_dim1 + 1], ldb, &c_b1, &a[k * a_dim1 + 1], lda, (ftnlen)1,
                           (ftnlen)1);
                    i__3 = k - 1;
                    zher2k_(uplo, (char *)"N", &i__3, &kb, &c_b1, &a[k * a_dim1 + 1], lda,
                            &b[k * b_dim1 + 1], ldb, &c_b18, &a[a_offset], lda, (ftnlen)1,
                            (ftnlen)1);
                    i__3 = k - 1;
                    zhemm_((char *)"R", uplo, &i__3, &kb, &c_b2, &a[k + k * a_dim1], lda,
                           &b[k * b_dim1 + 1], ldb, &c_b1, &a[k * a_dim1 + 1], lda, (ftnlen)1,
                           (ftnlen)1);
                    i__3 = k - 1;
                    ztrmm_((char *)"R", uplo, (char *)"C", (char *)"N", &i__3, &kb, &c_b1, &b[k + k * b_dim1], ldb,
                           &a[k * a_dim1 + 1], lda, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                    zhegs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + k * b_dim1], ldb,
                            info, (ftnlen)1);
                }
            } else {
                i__2 = *n;
                i__1 = nb;
                for (k = 1; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1) {
                    i__3 = *n - k + 1;
                    kb = min(i__3, nb);
                    i__3 = k - 1;
                    ztrmm_((char *)"R", uplo, (char *)"N", (char *)"N", &kb, &i__3, &c_b1, &b[b_offset], ldb,
                           &a[k + a_dim1], lda, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                    i__3 = k - 1;
                    zhemm_((char *)"L", uplo, &kb, &i__3, &c_b2, &a[k + k * a_dim1], lda, &b[k + b_dim1],
                           ldb, &c_b1, &a[k + a_dim1], lda, (ftnlen)1, (ftnlen)1);
                    i__3 = k - 1;
                    zher2k_(uplo, (char *)"C", &i__3, &kb, &c_b1, &a[k + a_dim1], lda, &b[k + b_dim1], ldb,
                            &c_b18, &a[a_offset], lda, (ftnlen)1, (ftnlen)1);
                    i__3 = k - 1;
                    zhemm_((char *)"L", uplo, &kb, &i__3, &c_b2, &a[k + k * a_dim1], lda, &b[k + b_dim1],
                           ldb, &c_b1, &a[k + a_dim1], lda, (ftnlen)1, (ftnlen)1);
                    i__3 = k - 1;
                    ztrmm_((char *)"L", uplo, (char *)"C", (char *)"N", &kb, &i__3, &c_b1, &b[k + k * b_dim1], ldb,
                           &a[k + a_dim1], lda, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                    zhegs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + k * b_dim1], ldb,
                            info, (ftnlen)1);
                }
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
