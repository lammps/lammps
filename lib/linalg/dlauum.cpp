#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b15 = 1.;
int dlauum_(char *uplo, integer *n, doublereal *a, integer *lda, integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    integer i__, ib, nb;
    extern int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
                      integer *, doublereal *, integer *, doublereal *, doublereal *, integer *,
                      ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int dtrmm_(char *, char *, char *, char *, integer *, integer *, doublereal *,
                      doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen,
                      ftnlen);
    logical upper;
    extern int dsyrk_(char *, char *, integer *, integer *, doublereal *, doublereal *, integer *,
                      doublereal *, doublereal *, integer *, ftnlen, ftnlen),
        dlauu2_(char *, integer *, doublereal *, integer *, integer *, ftnlen),
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
        xerbla_((char *)"DLAUUM", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    nb = ilaenv_(&c__1, (char *)"DLAUUM", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
    if (nb <= 1 || nb >= *n) {
        dlauu2_(uplo, n, &a[a_offset], lda, info, (ftnlen)1);
    } else {
        if (upper) {
            i__1 = *n;
            i__2 = nb;
            for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
                i__3 = nb, i__4 = *n - i__ + 1;
                ib = min(i__3, i__4);
                i__3 = i__ - 1;
                dtrmm_((char *)"Right", (char *)"Upper", (char *)"Transpose", (char *)"Non-unit", &i__3, &ib, &c_b15,
                       &a[i__ + i__ * a_dim1], lda, &a[i__ * a_dim1 + 1], lda, (ftnlen)5, (ftnlen)5,
                       (ftnlen)9, (ftnlen)8);
                dlauu2_((char *)"Upper", &ib, &a[i__ + i__ * a_dim1], lda, info, (ftnlen)5);
                if (i__ + ib <= *n) {
                    i__3 = i__ - 1;
                    i__4 = *n - i__ - ib + 1;
                    dgemm_((char *)"No transpose", (char *)"Transpose", &i__3, &ib, &i__4, &c_b15,
                           &a[(i__ + ib) * a_dim1 + 1], lda, &a[i__ + (i__ + ib) * a_dim1], lda,
                           &c_b15, &a[i__ * a_dim1 + 1], lda, (ftnlen)12, (ftnlen)9);
                    i__3 = *n - i__ - ib + 1;
                    dsyrk_((char *)"Upper", (char *)"No transpose", &ib, &i__3, &c_b15,
                           &a[i__ + (i__ + ib) * a_dim1], lda, &c_b15, &a[i__ + i__ * a_dim1], lda,
                           (ftnlen)5, (ftnlen)12);
                }
            }
        } else {
            i__2 = *n;
            i__1 = nb;
            for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
                i__3 = nb, i__4 = *n - i__ + 1;
                ib = min(i__3, i__4);
                i__3 = i__ - 1;
                dtrmm_((char *)"Left", (char *)"Lower", (char *)"Transpose", (char *)"Non-unit", &ib, &i__3, &c_b15,
                       &a[i__ + i__ * a_dim1], lda, &a[i__ + a_dim1], lda, (ftnlen)4, (ftnlen)5,
                       (ftnlen)9, (ftnlen)8);
                dlauu2_((char *)"Lower", &ib, &a[i__ + i__ * a_dim1], lda, info, (ftnlen)5);
                if (i__ + ib <= *n) {
                    i__3 = i__ - 1;
                    i__4 = *n - i__ - ib + 1;
                    dgemm_((char *)"Transpose", (char *)"No transpose", &ib, &i__3, &i__4, &c_b15,
                           &a[i__ + ib + i__ * a_dim1], lda, &a[i__ + ib + a_dim1], lda, &c_b15,
                           &a[i__ + a_dim1], lda, (ftnlen)9, (ftnlen)12);
                    i__3 = *n - i__ - ib + 1;
                    dsyrk_((char *)"Lower", (char *)"Transpose", &ib, &i__3, &c_b15, &a[i__ + ib + i__ * a_dim1],
                           lda, &c_b15, &a[i__ + i__ * a_dim1], lda, (ftnlen)5, (ftnlen)9);
                }
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
