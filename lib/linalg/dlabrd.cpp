#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b4 = -1.;
static doublereal c_b5 = 1.;
static integer c__1 = 1;
static doublereal c_b16 = 0.;
int dlabrd_(integer *m, integer *n, integer *nb, doublereal *a, integer *lda, doublereal *d__,
            doublereal *e, doublereal *tauq, doublereal *taup, doublereal *x, integer *ldx,
            doublereal *y, integer *ldy)
{
    integer a_dim1, a_offset, x_dim1, x_offset, y_dim1, y_offset, i__1, i__2, i__3;
    integer i__;
    extern int dscal_(integer *, doublereal *, doublereal *, integer *),
        dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *,
               integer *, doublereal *, doublereal *, integer *, ftnlen),
        dlarfg_(integer *, doublereal *, doublereal *, integer *, doublereal *);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tauq;
    --taup;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    if (*m <= 0 || *n <= 0) {
        return 0;
    }
    if (*m >= *n) {
        i__1 = *nb;
        for (i__ = 1; i__ <= i__1; ++i__) {
            i__2 = *m - i__ + 1;
            i__3 = i__ - 1;
            dgemv_((char *)"No transpose", &i__2, &i__3, &c_b4, &a[i__ + a_dim1], lda, &y[i__ + y_dim1],
                   ldy, &c_b5, &a[i__ + i__ * a_dim1], &c__1, (ftnlen)12);
            i__2 = *m - i__ + 1;
            i__3 = i__ - 1;
            dgemv_((char *)"No transpose", &i__2, &i__3, &c_b4, &x[i__ + x_dim1], ldx, &a[i__ * a_dim1 + 1],
                   &c__1, &c_b5, &a[i__ + i__ * a_dim1], &c__1, (ftnlen)12);
            i__2 = *m - i__ + 1;
            i__3 = i__ + 1;
            dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[min(i__3, *m) + i__ * a_dim1], &c__1,
                    &tauq[i__]);
            d__[i__] = a[i__ + i__ * a_dim1];
            if (i__ < *n) {
                a[i__ + i__ * a_dim1] = 1.;
                i__2 = *m - i__ + 1;
                i__3 = *n - i__;
                dgemv_((char *)"Transpose", &i__2, &i__3, &c_b5, &a[i__ + (i__ + 1) * a_dim1], lda,
                       &a[i__ + i__ * a_dim1], &c__1, &c_b16, &y[i__ + 1 + i__ * y_dim1], &c__1,
                       (ftnlen)9);
                i__2 = *m - i__ + 1;
                i__3 = i__ - 1;
                dgemv_((char *)"Transpose", &i__2, &i__3, &c_b5, &a[i__ + a_dim1], lda,
                       &a[i__ + i__ * a_dim1], &c__1, &c_b16, &y[i__ * y_dim1 + 1], &c__1,
                       (ftnlen)9);
                i__2 = *n - i__;
                i__3 = i__ - 1;
                dgemv_((char *)"No transpose", &i__2, &i__3, &c_b4, &y[i__ + 1 + y_dim1], ldy,
                       &y[i__ * y_dim1 + 1], &c__1, &c_b5, &y[i__ + 1 + i__ * y_dim1], &c__1,
                       (ftnlen)12);
                i__2 = *m - i__ + 1;
                i__3 = i__ - 1;
                dgemv_((char *)"Transpose", &i__2, &i__3, &c_b5, &x[i__ + x_dim1], ldx,
                       &a[i__ + i__ * a_dim1], &c__1, &c_b16, &y[i__ * y_dim1 + 1], &c__1,
                       (ftnlen)9);
                i__2 = i__ - 1;
                i__3 = *n - i__;
                dgemv_((char *)"Transpose", &i__2, &i__3, &c_b4, &a[(i__ + 1) * a_dim1 + 1], lda,
                       &y[i__ * y_dim1 + 1], &c__1, &c_b5, &y[i__ + 1 + i__ * y_dim1], &c__1,
                       (ftnlen)9);
                i__2 = *n - i__;
                dscal_(&i__2, &tauq[i__], &y[i__ + 1 + i__ * y_dim1], &c__1);
                i__2 = *n - i__;
                dgemv_((char *)"No transpose", &i__2, &i__, &c_b4, &y[i__ + 1 + y_dim1], ldy,
                       &a[i__ + a_dim1], lda, &c_b5, &a[i__ + (i__ + 1) * a_dim1], lda, (ftnlen)12);
                i__2 = i__ - 1;
                i__3 = *n - i__;
                dgemv_((char *)"Transpose", &i__2, &i__3, &c_b4, &a[(i__ + 1) * a_dim1 + 1], lda,
                       &x[i__ + x_dim1], ldx, &c_b5, &a[i__ + (i__ + 1) * a_dim1], lda, (ftnlen)9);
                i__2 = *n - i__;
                i__3 = i__ + 2;
                dlarfg_(&i__2, &a[i__ + (i__ + 1) * a_dim1], &a[i__ + min(i__3, *n) * a_dim1], lda,
                        &taup[i__]);
                e[i__] = a[i__ + (i__ + 1) * a_dim1];
                a[i__ + (i__ + 1) * a_dim1] = 1.;
                i__2 = *m - i__;
                i__3 = *n - i__;
                dgemv_((char *)"No transpose", &i__2, &i__3, &c_b5, &a[i__ + 1 + (i__ + 1) * a_dim1], lda,
                       &a[i__ + (i__ + 1) * a_dim1], lda, &c_b16, &x[i__ + 1 + i__ * x_dim1], &c__1,
                       (ftnlen)12);
                i__2 = *n - i__;
                dgemv_((char *)"Transpose", &i__2, &i__, &c_b5, &y[i__ + 1 + y_dim1], ldy,
                       &a[i__ + (i__ + 1) * a_dim1], lda, &c_b16, &x[i__ * x_dim1 + 1], &c__1,
                       (ftnlen)9);
                i__2 = *m - i__;
                dgemv_((char *)"No transpose", &i__2, &i__, &c_b4, &a[i__ + 1 + a_dim1], lda,
                       &x[i__ * x_dim1 + 1], &c__1, &c_b5, &x[i__ + 1 + i__ * x_dim1], &c__1,
                       (ftnlen)12);
                i__2 = i__ - 1;
                i__3 = *n - i__;
                dgemv_((char *)"No transpose", &i__2, &i__3, &c_b5, &a[(i__ + 1) * a_dim1 + 1], lda,
                       &a[i__ + (i__ + 1) * a_dim1], lda, &c_b16, &x[i__ * x_dim1 + 1], &c__1,
                       (ftnlen)12);
                i__2 = *m - i__;
                i__3 = i__ - 1;
                dgemv_((char *)"No transpose", &i__2, &i__3, &c_b4, &x[i__ + 1 + x_dim1], ldx,
                       &x[i__ * x_dim1 + 1], &c__1, &c_b5, &x[i__ + 1 + i__ * x_dim1], &c__1,
                       (ftnlen)12);
                i__2 = *m - i__;
                dscal_(&i__2, &taup[i__], &x[i__ + 1 + i__ * x_dim1], &c__1);
            }
        }
    } else {
        i__1 = *nb;
        for (i__ = 1; i__ <= i__1; ++i__) {
            i__2 = *n - i__ + 1;
            i__3 = i__ - 1;
            dgemv_((char *)"No transpose", &i__2, &i__3, &c_b4, &y[i__ + y_dim1], ldy, &a[i__ + a_dim1],
                   lda, &c_b5, &a[i__ + i__ * a_dim1], lda, (ftnlen)12);
            i__2 = i__ - 1;
            i__3 = *n - i__ + 1;
            dgemv_((char *)"Transpose", &i__2, &i__3, &c_b4, &a[i__ * a_dim1 + 1], lda, &x[i__ + x_dim1],
                   ldx, &c_b5, &a[i__ + i__ * a_dim1], lda, (ftnlen)9);
            i__2 = *n - i__ + 1;
            i__3 = i__ + 1;
            dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[i__ + min(i__3, *n) * a_dim1], lda,
                    &taup[i__]);
            d__[i__] = a[i__ + i__ * a_dim1];
            if (i__ < *m) {
                a[i__ + i__ * a_dim1] = 1.;
                i__2 = *m - i__;
                i__3 = *n - i__ + 1;
                dgemv_((char *)"No transpose", &i__2, &i__3, &c_b5, &a[i__ + 1 + i__ * a_dim1], lda,
                       &a[i__ + i__ * a_dim1], lda, &c_b16, &x[i__ + 1 + i__ * x_dim1], &c__1,
                       (ftnlen)12);
                i__2 = *n - i__ + 1;
                i__3 = i__ - 1;
                dgemv_((char *)"Transpose", &i__2, &i__3, &c_b5, &y[i__ + y_dim1], ldy,
                       &a[i__ + i__ * a_dim1], lda, &c_b16, &x[i__ * x_dim1 + 1], &c__1, (ftnlen)9);
                i__2 = *m - i__;
                i__3 = i__ - 1;
                dgemv_((char *)"No transpose", &i__2, &i__3, &c_b4, &a[i__ + 1 + a_dim1], lda,
                       &x[i__ * x_dim1 + 1], &c__1, &c_b5, &x[i__ + 1 + i__ * x_dim1], &c__1,
                       (ftnlen)12);
                i__2 = i__ - 1;
                i__3 = *n - i__ + 1;
                dgemv_((char *)"No transpose", &i__2, &i__3, &c_b5, &a[i__ * a_dim1 + 1], lda,
                       &a[i__ + i__ * a_dim1], lda, &c_b16, &x[i__ * x_dim1 + 1], &c__1,
                       (ftnlen)12);
                i__2 = *m - i__;
                i__3 = i__ - 1;
                dgemv_((char *)"No transpose", &i__2, &i__3, &c_b4, &x[i__ + 1 + x_dim1], ldx,
                       &x[i__ * x_dim1 + 1], &c__1, &c_b5, &x[i__ + 1 + i__ * x_dim1], &c__1,
                       (ftnlen)12);
                i__2 = *m - i__;
                dscal_(&i__2, &taup[i__], &x[i__ + 1 + i__ * x_dim1], &c__1);
                i__2 = *m - i__;
                i__3 = i__ - 1;
                dgemv_((char *)"No transpose", &i__2, &i__3, &c_b4, &a[i__ + 1 + a_dim1], lda,
                       &y[i__ + y_dim1], ldy, &c_b5, &a[i__ + 1 + i__ * a_dim1], &c__1, (ftnlen)12);
                i__2 = *m - i__;
                dgemv_((char *)"No transpose", &i__2, &i__, &c_b4, &x[i__ + 1 + x_dim1], ldx,
                       &a[i__ * a_dim1 + 1], &c__1, &c_b5, &a[i__ + 1 + i__ * a_dim1], &c__1,
                       (ftnlen)12);
                i__2 = *m - i__;
                i__3 = i__ + 2;
                dlarfg_(&i__2, &a[i__ + 1 + i__ * a_dim1], &a[min(i__3, *m) + i__ * a_dim1], &c__1,
                        &tauq[i__]);
                e[i__] = a[i__ + 1 + i__ * a_dim1];
                a[i__ + 1 + i__ * a_dim1] = 1.;
                i__2 = *m - i__;
                i__3 = *n - i__;
                dgemv_((char *)"Transpose", &i__2, &i__3, &c_b5, &a[i__ + 1 + (i__ + 1) * a_dim1], lda,
                       &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b16, &y[i__ + 1 + i__ * y_dim1], &c__1,
                       (ftnlen)9);
                i__2 = *m - i__;
                i__3 = i__ - 1;
                dgemv_((char *)"Transpose", &i__2, &i__3, &c_b5, &a[i__ + 1 + a_dim1], lda,
                       &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b16, &y[i__ * y_dim1 + 1], &c__1,
                       (ftnlen)9);
                i__2 = *n - i__;
                i__3 = i__ - 1;
                dgemv_((char *)"No transpose", &i__2, &i__3, &c_b4, &y[i__ + 1 + y_dim1], ldy,
                       &y[i__ * y_dim1 + 1], &c__1, &c_b5, &y[i__ + 1 + i__ * y_dim1], &c__1,
                       (ftnlen)12);
                i__2 = *m - i__;
                dgemv_((char *)"Transpose", &i__2, &i__, &c_b5, &x[i__ + 1 + x_dim1], ldx,
                       &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b16, &y[i__ * y_dim1 + 1], &c__1,
                       (ftnlen)9);
                i__2 = *n - i__;
                dgemv_((char *)"Transpose", &i__, &i__2, &c_b4, &a[(i__ + 1) * a_dim1 + 1], lda,
                       &y[i__ * y_dim1 + 1], &c__1, &c_b5, &y[i__ + 1 + i__ * y_dim1], &c__1,
                       (ftnlen)9);
                i__2 = *n - i__;
                dscal_(&i__2, &tauq[i__], &y[i__ + 1 + i__ * y_dim1], &c__1);
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
