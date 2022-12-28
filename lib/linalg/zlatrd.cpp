#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublecomplex c_b1 = {0., 0.};
static doublecomplex c_b2 = {1., 0.};
static integer c__1 = 1;
int zlatrd_(char *uplo, integer *n, integer *nb, doublecomplex *a, integer *lda, doublereal *e,
            doublecomplex *tau, doublecomplex *w, integer *ldw, ftnlen uplo_len)
{
    integer a_dim1, a_offset, w_dim1, w_offset, i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4;
    integer i__, iw;
    doublecomplex alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int zscal_(integer *, doublecomplex *, doublecomplex *, integer *);
    extern VOID zdotc_(doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *,
                       integer *);
    extern int zgemv_(char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *,
                      doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *,
                      ftnlen),
        zhemv_(char *, integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *,
               integer *, doublecomplex *, doublecomplex *, integer *, ftnlen),
        zaxpy_(integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *),
        zlarfg_(integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *),
        zlacgv_(integer *, doublecomplex *, integer *);
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
                i__2 = i__ + i__ * a_dim1;
                i__3 = i__ + i__ * a_dim1;
                d__1 = a[i__3].r;
                a[i__2].r = d__1, a[i__2].i = 0.;
                i__2 = *n - i__;
                zlacgv_(&i__2, &w[i__ + (iw + 1) * w_dim1], ldw);
                i__2 = *n - i__;
                z__1.r = -1., z__1.i = -0.;
                zgemv_((char *)"No transpose", &i__, &i__2, &z__1, &a[(i__ + 1) * a_dim1 + 1], lda,
                       &w[i__ + (iw + 1) * w_dim1], ldw, &c_b2, &a[i__ * a_dim1 + 1], &c__1,
                       (ftnlen)12);
                i__2 = *n - i__;
                zlacgv_(&i__2, &w[i__ + (iw + 1) * w_dim1], ldw);
                i__2 = *n - i__;
                zlacgv_(&i__2, &a[i__ + (i__ + 1) * a_dim1], lda);
                i__2 = *n - i__;
                z__1.r = -1., z__1.i = -0.;
                zgemv_((char *)"No transpose", &i__, &i__2, &z__1, &w[(iw + 1) * w_dim1 + 1], ldw,
                       &a[i__ + (i__ + 1) * a_dim1], lda, &c_b2, &a[i__ * a_dim1 + 1], &c__1,
                       (ftnlen)12);
                i__2 = *n - i__;
                zlacgv_(&i__2, &a[i__ + (i__ + 1) * a_dim1], lda);
                i__2 = i__ + i__ * a_dim1;
                i__3 = i__ + i__ * a_dim1;
                d__1 = a[i__3].r;
                a[i__2].r = d__1, a[i__2].i = 0.;
            }
            if (i__ > 1) {
                i__2 = i__ - 1 + i__ * a_dim1;
                alpha.r = a[i__2].r, alpha.i = a[i__2].i;
                i__2 = i__ - 1;
                zlarfg_(&i__2, &alpha, &a[i__ * a_dim1 + 1], &c__1, &tau[i__ - 1]);
                e[i__ - 1] = alpha.r;
                i__2 = i__ - 1 + i__ * a_dim1;
                a[i__2].r = 1., a[i__2].i = 0.;
                i__2 = i__ - 1;
                zhemv_((char *)"Upper", &i__2, &c_b2, &a[a_offset], lda, &a[i__ * a_dim1 + 1], &c__1, &c_b1,
                       &w[iw * w_dim1 + 1], &c__1, (ftnlen)5);
                if (i__ < *n) {
                    i__2 = i__ - 1;
                    i__3 = *n - i__;
                    zgemv_((char *)"Conjugate transpose", &i__2, &i__3, &c_b2, &w[(iw + 1) * w_dim1 + 1],
                           ldw, &a[i__ * a_dim1 + 1], &c__1, &c_b1, &w[i__ + 1 + iw * w_dim1],
                           &c__1, (ftnlen)19);
                    i__2 = i__ - 1;
                    i__3 = *n - i__;
                    z__1.r = -1., z__1.i = -0.;
                    zgemv_((char *)"No transpose", &i__2, &i__3, &z__1, &a[(i__ + 1) * a_dim1 + 1], lda,
                           &w[i__ + 1 + iw * w_dim1], &c__1, &c_b2, &w[iw * w_dim1 + 1], &c__1,
                           (ftnlen)12);
                    i__2 = i__ - 1;
                    i__3 = *n - i__;
                    zgemv_((char *)"Conjugate transpose", &i__2, &i__3, &c_b2, &a[(i__ + 1) * a_dim1 + 1],
                           lda, &a[i__ * a_dim1 + 1], &c__1, &c_b1, &w[i__ + 1 + iw * w_dim1],
                           &c__1, (ftnlen)19);
                    i__2 = i__ - 1;
                    i__3 = *n - i__;
                    z__1.r = -1., z__1.i = -0.;
                    zgemv_((char *)"No transpose", &i__2, &i__3, &z__1, &w[(iw + 1) * w_dim1 + 1], ldw,
                           &w[i__ + 1 + iw * w_dim1], &c__1, &c_b2, &w[iw * w_dim1 + 1], &c__1,
                           (ftnlen)12);
                }
                i__2 = i__ - 1;
                zscal_(&i__2, &tau[i__ - 1], &w[iw * w_dim1 + 1], &c__1);
                z__3.r = -.5, z__3.i = -0.;
                i__2 = i__ - 1;
                z__2.r = z__3.r * tau[i__2].r - z__3.i * tau[i__2].i,
                z__2.i = z__3.r * tau[i__2].i + z__3.i * tau[i__2].r;
                i__3 = i__ - 1;
                zdotc_(&z__4, &i__3, &w[iw * w_dim1 + 1], &c__1, &a[i__ * a_dim1 + 1], &c__1);
                z__1.r = z__2.r * z__4.r - z__2.i * z__4.i,
                z__1.i = z__2.r * z__4.i + z__2.i * z__4.r;
                alpha.r = z__1.r, alpha.i = z__1.i;
                i__2 = i__ - 1;
                zaxpy_(&i__2, &alpha, &a[i__ * a_dim1 + 1], &c__1, &w[iw * w_dim1 + 1], &c__1);
            }
        }
    } else {
        i__1 = *nb;
        for (i__ = 1; i__ <= i__1; ++i__) {
            i__2 = i__ + i__ * a_dim1;
            i__3 = i__ + i__ * a_dim1;
            d__1 = a[i__3].r;
            a[i__2].r = d__1, a[i__2].i = 0.;
            i__2 = i__ - 1;
            zlacgv_(&i__2, &w[i__ + w_dim1], ldw);
            i__2 = *n - i__ + 1;
            i__3 = i__ - 1;
            z__1.r = -1., z__1.i = -0.;
            zgemv_((char *)"No transpose", &i__2, &i__3, &z__1, &a[i__ + a_dim1], lda, &w[i__ + w_dim1],
                   ldw, &c_b2, &a[i__ + i__ * a_dim1], &c__1, (ftnlen)12);
            i__2 = i__ - 1;
            zlacgv_(&i__2, &w[i__ + w_dim1], ldw);
            i__2 = i__ - 1;
            zlacgv_(&i__2, &a[i__ + a_dim1], lda);
            i__2 = *n - i__ + 1;
            i__3 = i__ - 1;
            z__1.r = -1., z__1.i = -0.;
            zgemv_((char *)"No transpose", &i__2, &i__3, &z__1, &w[i__ + w_dim1], ldw, &a[i__ + a_dim1],
                   lda, &c_b2, &a[i__ + i__ * a_dim1], &c__1, (ftnlen)12);
            i__2 = i__ - 1;
            zlacgv_(&i__2, &a[i__ + a_dim1], lda);
            i__2 = i__ + i__ * a_dim1;
            i__3 = i__ + i__ * a_dim1;
            d__1 = a[i__3].r;
            a[i__2].r = d__1, a[i__2].i = 0.;
            if (i__ < *n) {
                i__2 = i__ + 1 + i__ * a_dim1;
                alpha.r = a[i__2].r, alpha.i = a[i__2].i;
                i__2 = *n - i__;
                i__3 = i__ + 2;
                zlarfg_(&i__2, &alpha, &a[min(i__3, *n) + i__ * a_dim1], &c__1, &tau[i__]);
                e[i__] = alpha.r;
                i__2 = i__ + 1 + i__ * a_dim1;
                a[i__2].r = 1., a[i__2].i = 0.;
                i__2 = *n - i__;
                zhemv_((char *)"Lower", &i__2, &c_b2, &a[i__ + 1 + (i__ + 1) * a_dim1], lda,
                       &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b1, &w[i__ + 1 + i__ * w_dim1], &c__1,
                       (ftnlen)5);
                i__2 = *n - i__;
                i__3 = i__ - 1;
                zgemv_((char *)"Conjugate transpose", &i__2, &i__3, &c_b2, &w[i__ + 1 + w_dim1], ldw,
                       &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b1, &w[i__ * w_dim1 + 1], &c__1,
                       (ftnlen)19);
                i__2 = *n - i__;
                i__3 = i__ - 1;
                z__1.r = -1., z__1.i = -0.;
                zgemv_((char *)"No transpose", &i__2, &i__3, &z__1, &a[i__ + 1 + a_dim1], lda,
                       &w[i__ * w_dim1 + 1], &c__1, &c_b2, &w[i__ + 1 + i__ * w_dim1], &c__1,
                       (ftnlen)12);
                i__2 = *n - i__;
                i__3 = i__ - 1;
                zgemv_((char *)"Conjugate transpose", &i__2, &i__3, &c_b2, &a[i__ + 1 + a_dim1], lda,
                       &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b1, &w[i__ * w_dim1 + 1], &c__1,
                       (ftnlen)19);
                i__2 = *n - i__;
                i__3 = i__ - 1;
                z__1.r = -1., z__1.i = -0.;
                zgemv_((char *)"No transpose", &i__2, &i__3, &z__1, &w[i__ + 1 + w_dim1], ldw,
                       &w[i__ * w_dim1 + 1], &c__1, &c_b2, &w[i__ + 1 + i__ * w_dim1], &c__1,
                       (ftnlen)12);
                i__2 = *n - i__;
                zscal_(&i__2, &tau[i__], &w[i__ + 1 + i__ * w_dim1], &c__1);
                z__3.r = -.5, z__3.i = -0.;
                i__2 = i__;
                z__2.r = z__3.r * tau[i__2].r - z__3.i * tau[i__2].i,
                z__2.i = z__3.r * tau[i__2].i + z__3.i * tau[i__2].r;
                i__3 = *n - i__;
                zdotc_(&z__4, &i__3, &w[i__ + 1 + i__ * w_dim1], &c__1, &a[i__ + 1 + i__ * a_dim1],
                       &c__1);
                z__1.r = z__2.r * z__4.r - z__2.i * z__4.i,
                z__1.i = z__2.r * z__4.i + z__2.i * z__4.r;
                alpha.r = z__1.r, alpha.i = z__1.i;
                i__2 = *n - i__;
                zaxpy_(&i__2, &alpha, &a[i__ + 1 + i__ * a_dim1], &c__1, &w[i__ + 1 + i__ * w_dim1],
                       &c__1);
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
