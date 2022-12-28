#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b5 = -1.;
static doublereal c_b6 = 1.;
static integer c__1 = 1;
static doublereal c_b16 = 0.;
int dlatrd_(char *uplo, integer *n, integer *nb, doublereal *a, integer *lda, doublereal *e,
            doublereal *tau, doublereal *w, integer *ldw, ftnlen uplo_len)
{
    integer a_dim1, a_offset, w_dim1, w_offset, i__1, i__2, i__3;
    integer i__, iw;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, integer *);
    doublereal alpha;
    extern int dscal_(integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *,
                      doublereal *, integer *, doublereal *, doublereal *, integer *, ftnlen),
        daxpy_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *),
        dsymv_(char *, integer *, doublereal *, doublereal *, integer *, doublereal *, integer *,
               doublereal *, doublereal *, integer *, ftnlen),
        dlarfg_(integer *, doublereal *, doublereal *, integer *, doublereal *);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --e;
    --tau;
    w_dim1 = *ldw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    if (*n <= 0) {
        return 0;
    }
    if (lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1)) {
        i__1 = *n - *nb + 1;
        for (i__ = *n; i__ >= i__1; --i__) {
            iw = i__ - *n + *nb;
            if (i__ < *n) {
                i__2 = *n - i__;
                dgemv_((char *)"No transpose", &i__, &i__2, &c_b5, &a[(i__ + 1) * a_dim1 + 1], lda,
                       &w[i__ + (iw + 1) * w_dim1], ldw, &c_b6, &a[i__ * a_dim1 + 1], &c__1,
                       (ftnlen)12);
                i__2 = *n - i__;
                dgemv_((char *)"No transpose", &i__, &i__2, &c_b5, &w[(iw + 1) * w_dim1 + 1], ldw,
                       &a[i__ + (i__ + 1) * a_dim1], lda, &c_b6, &a[i__ * a_dim1 + 1], &c__1,
                       (ftnlen)12);
            }
            if (i__ > 1) {
                i__2 = i__ - 1;
                dlarfg_(&i__2, &a[i__ - 1 + i__ * a_dim1], &a[i__ * a_dim1 + 1], &c__1,
                        &tau[i__ - 1]);
                e[i__ - 1] = a[i__ - 1 + i__ * a_dim1];
                a[i__ - 1 + i__ * a_dim1] = 1.;
                i__2 = i__ - 1;
                dsymv_((char *)"Upper", &i__2, &c_b6, &a[a_offset], lda, &a[i__ * a_dim1 + 1], &c__1,
                       &c_b16, &w[iw * w_dim1 + 1], &c__1, (ftnlen)5);
                if (i__ < *n) {
                    i__2 = i__ - 1;
                    i__3 = *n - i__;
                    dgemv_((char *)"Transpose", &i__2, &i__3, &c_b6, &w[(iw + 1) * w_dim1 + 1], ldw,
                           &a[i__ * a_dim1 + 1], &c__1, &c_b16, &w[i__ + 1 + iw * w_dim1], &c__1,
                           (ftnlen)9);
                    i__2 = i__ - 1;
                    i__3 = *n - i__;
                    dgemv_((char *)"No transpose", &i__2, &i__3, &c_b5, &a[(i__ + 1) * a_dim1 + 1], lda,
                           &w[i__ + 1 + iw * w_dim1], &c__1, &c_b6, &w[iw * w_dim1 + 1], &c__1,
                           (ftnlen)12);
                    i__2 = i__ - 1;
                    i__3 = *n - i__;
                    dgemv_((char *)"Transpose", &i__2, &i__3, &c_b6, &a[(i__ + 1) * a_dim1 + 1], lda,
                           &a[i__ * a_dim1 + 1], &c__1, &c_b16, &w[i__ + 1 + iw * w_dim1], &c__1,
                           (ftnlen)9);
                    i__2 = i__ - 1;
                    i__3 = *n - i__;
                    dgemv_((char *)"No transpose", &i__2, &i__3, &c_b5, &w[(iw + 1) * w_dim1 + 1], ldw,
                           &w[i__ + 1 + iw * w_dim1], &c__1, &c_b6, &w[iw * w_dim1 + 1], &c__1,
                           (ftnlen)12);
                }
                i__2 = i__ - 1;
                dscal_(&i__2, &tau[i__ - 1], &w[iw * w_dim1 + 1], &c__1);
                i__2 = i__ - 1;
                alpha = tau[i__ - 1] * -.5 *
                        ddot_(&i__2, &w[iw * w_dim1 + 1], &c__1, &a[i__ * a_dim1 + 1], &c__1);
                i__2 = i__ - 1;
                daxpy_(&i__2, &alpha, &a[i__ * a_dim1 + 1], &c__1, &w[iw * w_dim1 + 1], &c__1);
            }
        }
    } else {
        i__1 = *nb;
        for (i__ = 1; i__ <= i__1; ++i__) {
            i__2 = *n - i__ + 1;
            i__3 = i__ - 1;
            dgemv_((char *)"No transpose", &i__2, &i__3, &c_b5, &a[i__ + a_dim1], lda, &w[i__ + w_dim1],
                   ldw, &c_b6, &a[i__ + i__ * a_dim1], &c__1, (ftnlen)12);
            i__2 = *n - i__ + 1;
            i__3 = i__ - 1;
            dgemv_((char *)"No transpose", &i__2, &i__3, &c_b5, &w[i__ + w_dim1], ldw, &a[i__ + a_dim1],
                   lda, &c_b6, &a[i__ + i__ * a_dim1], &c__1, (ftnlen)12);
            if (i__ < *n) {
                i__2 = *n - i__;
                i__3 = i__ + 2;
                dlarfg_(&i__2, &a[i__ + 1 + i__ * a_dim1], &a[min(i__3, *n) + i__ * a_dim1], &c__1,
                        &tau[i__]);
                e[i__] = a[i__ + 1 + i__ * a_dim1];
                a[i__ + 1 + i__ * a_dim1] = 1.;
                i__2 = *n - i__;
                dsymv_((char *)"Lower", &i__2, &c_b6, &a[i__ + 1 + (i__ + 1) * a_dim1], lda,
                       &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b16, &w[i__ + 1 + i__ * w_dim1], &c__1,
                       (ftnlen)5);
                i__2 = *n - i__;
                i__3 = i__ - 1;
                dgemv_((char *)"Transpose", &i__2, &i__3, &c_b6, &w[i__ + 1 + w_dim1], ldw,
                       &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b16, &w[i__ * w_dim1 + 1], &c__1,
                       (ftnlen)9);
                i__2 = *n - i__;
                i__3 = i__ - 1;
                dgemv_((char *)"No transpose", &i__2, &i__3, &c_b5, &a[i__ + 1 + a_dim1], lda,
                       &w[i__ * w_dim1 + 1], &c__1, &c_b6, &w[i__ + 1 + i__ * w_dim1], &c__1,
                       (ftnlen)12);
                i__2 = *n - i__;
                i__3 = i__ - 1;
                dgemv_((char *)"Transpose", &i__2, &i__3, &c_b6, &a[i__ + 1 + a_dim1], lda,
                       &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b16, &w[i__ * w_dim1 + 1], &c__1,
                       (ftnlen)9);
                i__2 = *n - i__;
                i__3 = i__ - 1;
                dgemv_((char *)"No transpose", &i__2, &i__3, &c_b5, &w[i__ + 1 + w_dim1], ldw,
                       &w[i__ * w_dim1 + 1], &c__1, &c_b6, &w[i__ + 1 + i__ * w_dim1], &c__1,
                       (ftnlen)12);
                i__2 = *n - i__;
                dscal_(&i__2, &tau[i__], &w[i__ + 1 + i__ * w_dim1], &c__1);
                i__2 = *n - i__;
                alpha = tau[i__] * -.5 *
                        ddot_(&i__2, &w[i__ + 1 + i__ * w_dim1], &c__1, &a[i__ + 1 + i__ * a_dim1],
                              &c__1);
                i__2 = *n - i__;
                daxpy_(&i__2, &alpha, &a[i__ + 1 + i__ * a_dim1], &c__1, &w[i__ + 1 + i__ * w_dim1],
                       &c__1);
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
