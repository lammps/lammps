#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static doublereal c_b8 = 0.;
static doublereal c_b14 = -1.;
int dsytd2_(char *uplo, integer *n, doublereal *a, integer *lda, doublereal *d__, doublereal *e,
            doublereal *tau, integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, i__1, i__2, i__3;
    integer i__;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, integer *);
    doublereal taui;
    extern int dsyr2_(char *, integer *, doublereal *, doublereal *, integer *, doublereal *,
                      integer *, doublereal *, integer *, ftnlen);
    doublereal alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int daxpy_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *);
    logical upper;
    extern int dsymv_(char *, integer *, doublereal *, doublereal *, integer *, doublereal *,
                      integer *, doublereal *, doublereal *, integer *, ftnlen),
        dlarfg_(integer *, doublereal *, doublereal *, integer *, doublereal *),
        xerbla_(char *, integer *, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tau;
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
        xerbla_((char *)"DSYTD2", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n <= 0) {
        return 0;
    }
    if (upper) {
        for (i__ = *n - 1; i__ >= 1; --i__) {
            dlarfg_(&i__, &a[i__ + (i__ + 1) * a_dim1], &a[(i__ + 1) * a_dim1 + 1], &c__1, &taui);
            e[i__] = a[i__ + (i__ + 1) * a_dim1];
            if (taui != 0.) {
                a[i__ + (i__ + 1) * a_dim1] = 1.;
                dsymv_(uplo, &i__, &taui, &a[a_offset], lda, &a[(i__ + 1) * a_dim1 + 1], &c__1,
                       &c_b8, &tau[1], &c__1, (ftnlen)1);
                alpha = taui * -.5 * ddot_(&i__, &tau[1], &c__1, &a[(i__ + 1) * a_dim1 + 1], &c__1);
                daxpy_(&i__, &alpha, &a[(i__ + 1) * a_dim1 + 1], &c__1, &tau[1], &c__1);
                dsyr2_(uplo, &i__, &c_b14, &a[(i__ + 1) * a_dim1 + 1], &c__1, &tau[1], &c__1,
                       &a[a_offset], lda, (ftnlen)1);
                a[i__ + (i__ + 1) * a_dim1] = e[i__];
            }
            d__[i__ + 1] = a[i__ + 1 + (i__ + 1) * a_dim1];
            tau[i__] = taui;
        }
        d__[1] = a[a_dim1 + 1];
    } else {
        i__1 = *n - 1;
        for (i__ = 1; i__ <= i__1; ++i__) {
            i__2 = *n - i__;
            i__3 = i__ + 2;
            dlarfg_(&i__2, &a[i__ + 1 + i__ * a_dim1], &a[min(i__3, *n) + i__ * a_dim1], &c__1,
                    &taui);
            e[i__] = a[i__ + 1 + i__ * a_dim1];
            if (taui != 0.) {
                a[i__ + 1 + i__ * a_dim1] = 1.;
                i__2 = *n - i__;
                dsymv_(uplo, &i__2, &taui, &a[i__ + 1 + (i__ + 1) * a_dim1], lda,
                       &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b8, &tau[i__], &c__1, (ftnlen)1);
                i__2 = *n - i__;
                alpha =
                    taui * -.5 * ddot_(&i__2, &tau[i__], &c__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
                i__2 = *n - i__;
                daxpy_(&i__2, &alpha, &a[i__ + 1 + i__ * a_dim1], &c__1, &tau[i__], &c__1);
                i__2 = *n - i__;
                dsyr2_(uplo, &i__2, &c_b14, &a[i__ + 1 + i__ * a_dim1], &c__1, &tau[i__], &c__1,
                       &a[i__ + 1 + (i__ + 1) * a_dim1], lda, (ftnlen)1);
                a[i__ + 1 + i__ * a_dim1] = e[i__];
            }
            d__[i__] = a[i__ + i__ * a_dim1];
            tau[i__] = taui;
        }
        d__[*n] = a[*n + *n * a_dim1];
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
