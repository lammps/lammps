#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b14 = 1.;
static doublereal c_b16 = -.5;
static doublereal c_b19 = -1.;
static doublereal c_b52 = .5;
int dsygst_(integer *itype, char *uplo, integer *n, doublereal *a, integer *lda, doublereal *b,
            integer *ldb, integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    integer k, kb, nb;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int dtrmm_(char *, char *, char *, char *, integer *, integer *, doublereal *,
                      doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen,
                      ftnlen),
        dsymm_(char *, char *, integer *, integer *, doublereal *, doublereal *, integer *,
               doublereal *, integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    logical upper;
    extern int dtrsm_(char *, char *, char *, char *, integer *, integer *, doublereal *,
                      doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen,
                      ftnlen),
        dsygs2_(integer *, char *, integer *, doublereal *, integer *, doublereal *, integer *,
                integer *, ftnlen),
        dsyr2k_(char *, char *, integer *, integer *, doublereal *, doublereal *, integer *,
                doublereal *, integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
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
        xerbla_((char *)"DSYGST", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    nb = ilaenv_(&c__1, (char *)"DSYGST", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
    if (nb <= 1 || nb >= *n) {
        dsygs2_(itype, uplo, n, &a[a_offset], lda, &b[b_offset], ldb, info, (ftnlen)1);
    } else {
        if (*itype == 1) {
            if (upper) {
                i__1 = *n;
                i__2 = nb;
                for (k = 1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2) {
                    i__3 = *n - k + 1;
                    kb = min(i__3, nb);
                    dsygs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + k * b_dim1], ldb,
                            info, (ftnlen)1);
                    if (k + kb <= *n) {
                        i__3 = *n - k - kb + 1;
                        dtrsm_((char *)"Left", uplo, (char *)"Transpose", (char *)"Non-unit", &kb, &i__3, &c_b14,
                               &b[k + k * b_dim1], ldb, &a[k + (k + kb) * a_dim1], lda, (ftnlen)4,
                               (ftnlen)1, (ftnlen)9, (ftnlen)8);
                        i__3 = *n - k - kb + 1;
                        dsymm_((char *)"Left", uplo, &kb, &i__3, &c_b16, &a[k + k * a_dim1], lda,
                               &b[k + (k + kb) * b_dim1], ldb, &c_b14, &a[k + (k + kb) * a_dim1],
                               lda, (ftnlen)4, (ftnlen)1);
                        i__3 = *n - k - kb + 1;
                        dsyr2k_(uplo, (char *)"Transpose", &i__3, &kb, &c_b19, &a[k + (k + kb) * a_dim1],
                                lda, &b[k + (k + kb) * b_dim1], ldb, &c_b14,
                                &a[k + kb + (k + kb) * a_dim1], lda, (ftnlen)1, (ftnlen)9);
                        i__3 = *n - k - kb + 1;
                        dsymm_((char *)"Left", uplo, &kb, &i__3, &c_b16, &a[k + k * a_dim1], lda,
                               &b[k + (k + kb) * b_dim1], ldb, &c_b14, &a[k + (k + kb) * a_dim1],
                               lda, (ftnlen)4, (ftnlen)1);
                        i__3 = *n - k - kb + 1;
                        dtrsm_((char *)"Right", uplo, (char *)"No transpose", (char *)"Non-unit", &kb, &i__3, &c_b14,
                               &b[k + kb + (k + kb) * b_dim1], ldb, &a[k + (k + kb) * a_dim1], lda,
                               (ftnlen)5, (ftnlen)1, (ftnlen)12, (ftnlen)8);
                    }
                }
            } else {
                i__2 = *n;
                i__1 = nb;
                for (k = 1; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1) {
                    i__3 = *n - k + 1;
                    kb = min(i__3, nb);
                    dsygs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + k * b_dim1], ldb,
                            info, (ftnlen)1);
                    if (k + kb <= *n) {
                        i__3 = *n - k - kb + 1;
                        dtrsm_((char *)"Right", uplo, (char *)"Transpose", (char *)"Non-unit", &i__3, &kb, &c_b14,
                               &b[k + k * b_dim1], ldb, &a[k + kb + k * a_dim1], lda, (ftnlen)5,
                               (ftnlen)1, (ftnlen)9, (ftnlen)8);
                        i__3 = *n - k - kb + 1;
                        dsymm_((char *)"Right", uplo, &i__3, &kb, &c_b16, &a[k + k * a_dim1], lda,
                               &b[k + kb + k * b_dim1], ldb, &c_b14, &a[k + kb + k * a_dim1], lda,
                               (ftnlen)5, (ftnlen)1);
                        i__3 = *n - k - kb + 1;
                        dsyr2k_(uplo, (char *)"No transpose", &i__3, &kb, &c_b19, &a[k + kb + k * a_dim1],
                                lda, &b[k + kb + k * b_dim1], ldb, &c_b14,
                                &a[k + kb + (k + kb) * a_dim1], lda, (ftnlen)1, (ftnlen)12);
                        i__3 = *n - k - kb + 1;
                        dsymm_((char *)"Right", uplo, &i__3, &kb, &c_b16, &a[k + k * a_dim1], lda,
                               &b[k + kb + k * b_dim1], ldb, &c_b14, &a[k + kb + k * a_dim1], lda,
                               (ftnlen)5, (ftnlen)1);
                        i__3 = *n - k - kb + 1;
                        dtrsm_((char *)"Left", uplo, (char *)"No transpose", (char *)"Non-unit", &i__3, &kb, &c_b14,
                               &b[k + kb + (k + kb) * b_dim1], ldb, &a[k + kb + k * a_dim1], lda,
                               (ftnlen)4, (ftnlen)1, (ftnlen)12, (ftnlen)8);
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
                    dtrmm_((char *)"Left", uplo, (char *)"No transpose", (char *)"Non-unit", &i__3, &kb, &c_b14,
                           &b[b_offset], ldb, &a[k * a_dim1 + 1], lda, (ftnlen)4, (ftnlen)1,
                           (ftnlen)12, (ftnlen)8);
                    i__3 = k - 1;
                    dsymm_((char *)"Right", uplo, &i__3, &kb, &c_b52, &a[k + k * a_dim1], lda,
                           &b[k * b_dim1 + 1], ldb, &c_b14, &a[k * a_dim1 + 1], lda, (ftnlen)5,
                           (ftnlen)1);
                    i__3 = k - 1;
                    dsyr2k_(uplo, (char *)"No transpose", &i__3, &kb, &c_b14, &a[k * a_dim1 + 1], lda,
                            &b[k * b_dim1 + 1], ldb, &c_b14, &a[a_offset], lda, (ftnlen)1,
                            (ftnlen)12);
                    i__3 = k - 1;
                    dsymm_((char *)"Right", uplo, &i__3, &kb, &c_b52, &a[k + k * a_dim1], lda,
                           &b[k * b_dim1 + 1], ldb, &c_b14, &a[k * a_dim1 + 1], lda, (ftnlen)5,
                           (ftnlen)1);
                    i__3 = k - 1;
                    dtrmm_((char *)"Right", uplo, (char *)"Transpose", (char *)"Non-unit", &i__3, &kb, &c_b14,
                           &b[k + k * b_dim1], ldb, &a[k * a_dim1 + 1], lda, (ftnlen)5, (ftnlen)1,
                           (ftnlen)9, (ftnlen)8);
                    dsygs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + k * b_dim1], ldb,
                            info, (ftnlen)1);
                }
            } else {
                i__2 = *n;
                i__1 = nb;
                for (k = 1; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1) {
                    i__3 = *n - k + 1;
                    kb = min(i__3, nb);
                    i__3 = k - 1;
                    dtrmm_((char *)"Right", uplo, (char *)"No transpose", (char *)"Non-unit", &kb, &i__3, &c_b14,
                           &b[b_offset], ldb, &a[k + a_dim1], lda, (ftnlen)5, (ftnlen)1, (ftnlen)12,
                           (ftnlen)8);
                    i__3 = k - 1;
                    dsymm_((char *)"Left", uplo, &kb, &i__3, &c_b52, &a[k + k * a_dim1], lda,
                           &b[k + b_dim1], ldb, &c_b14, &a[k + a_dim1], lda, (ftnlen)4, (ftnlen)1);
                    i__3 = k - 1;
                    dsyr2k_(uplo, (char *)"Transpose", &i__3, &kb, &c_b14, &a[k + a_dim1], lda,
                            &b[k + b_dim1], ldb, &c_b14, &a[a_offset], lda, (ftnlen)1, (ftnlen)9);
                    i__3 = k - 1;
                    dsymm_((char *)"Left", uplo, &kb, &i__3, &c_b52, &a[k + k * a_dim1], lda,
                           &b[k + b_dim1], ldb, &c_b14, &a[k + a_dim1], lda, (ftnlen)4, (ftnlen)1);
                    i__3 = k - 1;
                    dtrmm_((char *)"Left", uplo, (char *)"Transpose", (char *)"Non-unit", &kb, &i__3, &c_b14,
                           &b[k + k * b_dim1], ldb, &a[k + a_dim1], lda, (ftnlen)4, (ftnlen)1,
                           (ftnlen)9, (ftnlen)8);
                    dsygs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + k * b_dim1], ldb,
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
