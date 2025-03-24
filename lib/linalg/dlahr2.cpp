#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b4 = -1.;
static doublereal c_b5 = 1.;
static integer c__1 = 1;
static doublereal c_b38 = 0.;
int dlahr2_(integer *n, integer *k, integer *nb, doublereal *a, integer *lda, doublereal *tau,
            doublereal *t, integer *ldt, doublereal *y, integer *ldy)
{
    integer a_dim1, a_offset, t_dim1, t_offset, y_dim1, y_offset, i__1, i__2, i__3;
    doublereal d__1;
    integer i__;
    doublereal ei;
    extern int dscal_(integer *, doublereal *, doublereal *, integer *),
        dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
               integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, ftnlen,
               ftnlen),
        dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *,
               integer *, doublereal *, doublereal *, integer *, ftnlen),
        dcopy_(integer *, doublereal *, integer *, doublereal *, integer *),
        dtrmm_(char *, char *, char *, char *, integer *, integer *, doublereal *, doublereal *,
               integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen),
        daxpy_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *),
        dtrmv_(char *, char *, char *, integer *, doublereal *, integer *, doublereal *, integer *,
               ftnlen, ftnlen, ftnlen),
        dlarfg_(integer *, doublereal *, doublereal *, integer *, doublereal *),
        dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *,
                ftnlen);
    --tau;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    if (*n <= 1) {
        return 0;
    }
    i__1 = *nb;
    for (i__ = 1; i__ <= i__1; ++i__) {
        if (i__ > 1) {
            i__2 = *n - *k;
            i__3 = i__ - 1;
            dgemv_((char *)"T", &i__2, &i__3, &c_b4, &y[*k + 1 + y_dim1], ldy, &a[*k + i__ - 1 + a_dim1],
                   lda, &c_b5, &a[*k + 1 + i__ * a_dim1], &c__1, (ftnlen)1);
            i__2 = i__ - 1;
            dcopy_(&i__2, &a[*k + 1 + i__ * a_dim1], &c__1, &t[*nb * t_dim1 + 1], &c__1);
            i__2 = i__ - 1;
            dtrmv_((char *)"L", (char *)"T", (char *)"U", &i__2, &a[*k + 1 + a_dim1], lda, &t[*nb * t_dim1 + 1], &c__1,
                   (ftnlen)1, (ftnlen)1, (ftnlen)1);
            i__2 = *n - *k - i__ + 1;
            i__3 = i__ - 1;
            dgemv_((char *)"T", &i__2, &i__3, &c_b5, &a[*k + i__ + a_dim1], lda,
                   &a[*k + i__ + i__ * a_dim1], &c__1, &c_b5, &t[*nb * t_dim1 + 1], &c__1,
                   (ftnlen)1);
            i__2 = i__ - 1;
            dtrmv_((char *)"U", (char *)"T", (char *)"N", &i__2, &t[t_offset], ldt, &t[*nb * t_dim1 + 1], &c__1, (ftnlen)1,
                   (ftnlen)1, (ftnlen)1);
            i__2 = *n - *k - i__ + 1;
            i__3 = i__ - 1;
            dgemv_((char *)"T", &i__2, &i__3, &c_b4, &a[*k + i__ + a_dim1], lda, &t[*nb * t_dim1 + 1],
                   &c__1, &c_b5, &a[*k + i__ + i__ * a_dim1], &c__1, (ftnlen)1);
            i__2 = i__ - 1;
            dtrmv_((char *)"L", (char *)"T", (char *)"U", &i__2, &a[*k + 1 + a_dim1], lda, &t[*nb * t_dim1 + 1], &c__1,
                   (ftnlen)1, (ftnlen)1, (ftnlen)1);
            i__2 = i__ - 1;
            daxpy_(&i__2, &c_b4, &t[*nb * t_dim1 + 1], &c__1, &a[*k + 1 + i__ * a_dim1], &c__1);
            a[*k + i__ - 1 + (i__ - 1) * a_dim1] = ei;
        }
        i__2 = *n - *k - i__ + 1;
        i__3 = *k + i__ + 1;
        dlarfg_(&i__2, &a[*k + i__ + i__ * a_dim1], &a[min(i__3, *n) + i__ * a_dim1], &c__1,
                &tau[i__]);
        ei = a[*k + i__ + i__ * a_dim1];
        a[*k + i__ + i__ * a_dim1] = 1.;
        i__2 = *n - *k;
        i__3 = *n - *k - i__ + 1;
        dgemv_((char *)"T", &i__2, &i__3, &c_b5, &a[*k + 1 + (i__ + 1) * a_dim1], lda,
               &a[*k + i__ + i__ * a_dim1], &c__1, &c_b38, &y[*k + 1 + i__ * y_dim1], &c__1,
               (ftnlen)1);
        i__2 = *n - *k - i__ + 1;
        i__3 = i__ - 1;
        dgemv_((char *)"T", &i__2, &i__3, &c_b5, &a[*k + i__ + a_dim1], lda, &a[*k + i__ + i__ * a_dim1],
               &c__1, &c_b38, &t[i__ * t_dim1 + 1], &c__1, (ftnlen)1);
        i__2 = *n - *k;
        i__3 = i__ - 1;
        dgemv_((char *)"T", &i__2, &i__3, &c_b4, &y[*k + 1 + y_dim1], ldy, &t[i__ * t_dim1 + 1], &c__1,
               &c_b5, &y[*k + 1 + i__ * y_dim1], &c__1, (ftnlen)1);
        i__2 = *n - *k;
        dscal_(&i__2, &tau[i__], &y[*k + 1 + i__ * y_dim1], &c__1);
        i__2 = i__ - 1;
        d__1 = -tau[i__];
        dscal_(&i__2, &d__1, &t[i__ * t_dim1 + 1], &c__1);
        i__2 = i__ - 1;
        dtrmv_((char *)"U", (char *)"N", (char *)"N", &i__2, &t[t_offset], ldt, &t[i__ * t_dim1 + 1], &c__1, (ftnlen)1,
               (ftnlen)1, (ftnlen)1);
        t[i__ + i__ * t_dim1] = tau[i__];
    }
    a[*k + *nb + *nb * a_dim1] = ei;
    dlacpy_((char *)"A", k, nb, &a[(a_dim1 << 1) + 1], lda, &y[y_offset], ldy, (ftnlen)1);
    dtrmm_((char *)"R", (char *)"L", (char *)"T", (char *)"U", k, nb, &c_b5, &a[*k + 1 + a_dim1], lda, &y[y_offset], ldy, (ftnlen)1,
           (ftnlen)1, (ftnlen)1, (ftnlen)1);
    if (*n > *k + *nb) {
        i__1 = *n - *k - *nb;
        dgemm_((char *)"T", (char *)"T", k, nb, &i__1, &c_b5, &a[(*nb + 2) * a_dim1 + 1], lda,
               &a[*k + 1 + *nb + a_dim1], lda, &c_b5, &y[y_offset], ldy, (ftnlen)1, (ftnlen)1);
    }
    dtrmm_((char *)"R", (char *)"U", (char *)"T", (char *)"N", k, nb, &c_b5, &t[t_offset], ldt, &y[y_offset], ldy, (ftnlen)1,
           (ftnlen)1, (ftnlen)1, (ftnlen)1);
    return 0;
}
#ifdef __cplusplus
}
#endif
