#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublecomplex c_b1 = {1., 0.};
static doublecomplex c_b2 = {0., 0.};
static integer c__1 = 1;
int zlarf1l_(char *side, integer *m, integer *n, doublecomplex *v, integer *incv,
             doublecomplex *tau, doublecomplex *c__, integer *ldc, doublecomplex *work,
             ftnlen side_len)
{
    integer c_dim1, c_offset, i__1, i__2, i__3;
    doublecomplex z__1, z__2, z__3;
    void d_lmp_cnjg(doublecomplex *, doublecomplex *);
    integer i__, j;
    logical applyleft;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer lastc;
    extern int zgerc_(integer *, integer *, doublecomplex *, doublecomplex *, integer *,
                      doublecomplex *, integer *, doublecomplex *, integer *),
        zscal_(integer *, doublecomplex *, doublecomplex *, integer *),
        zgemv_(char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *,
               doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    integer lastv;
    extern int zaxpy_(integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *,
                      integer *);
    extern integer ilazlc_(integer *, integer *, doublecomplex *, integer *),
        ilazlr_(integer *, integer *, doublecomplex *, integer *);
    integer firstv;
    --v;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    applyleft = lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1);
    firstv = 1;
    lastc = 0;
    if (tau->r != 0. || tau->i != 0.) {
        if (applyleft) {
            lastv = *m;
        } else {
            lastv = *n;
        }
        i__ = 1;
        for (;;) {
            i__1 = i__;
            if (!(lastv > firstv && (v[i__1].r == 0. && v[i__1].i == 0.))) break;
            ++firstv;
            i__ += *incv;
        }
        if (applyleft) {
            lastc = ilazlc_(&lastv, n, &c__[c_offset], ldc);
        } else {
            lastc = ilazlr_(m, &lastv, &c__[c_offset], ldc);
        }
    }
    if (lastc == 0) {
        return 0;
    }
    if (applyleft) {
        if (lastv == firstv) {
            z__1.r = 1. - tau->r, z__1.i = 0. - tau->i;
            zscal_(&lastc, &z__1, &c__[lastv + c_dim1], ldc);
        } else {
            i__1 = lastv - firstv;
            zgemv_((char *)"C", &i__1, &lastc, &c_b1, &c__[firstv + c_dim1], ldc, &v[i__], incv, &c_b2,
                   &work[1], &c__1, (ftnlen)1);
            i__1 = lastc;
            for (j = 1; j <= i__1; ++j) {
                i__2 = j;
                i__3 = j;
                d_lmp_cnjg(&z__2, &c__[lastv + j * c_dim1]);
                z__1.r = work[i__3].r + z__2.r, z__1.i = work[i__3].i + z__2.i;
                work[i__2].r = z__1.r, work[i__2].i = z__1.i;
            }
            i__1 = lastc;
            for (j = 1; j <= i__1; ++j) {
                i__2 = lastv + j * c_dim1;
                i__3 = lastv + j * c_dim1;
                d_lmp_cnjg(&z__3, &work[j]);
                z__2.r = tau->r * z__3.r - tau->i * z__3.i,
                z__2.i = tau->r * z__3.i + tau->i * z__3.r;
                z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
                c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
            }
            i__1 = lastv - firstv;
            z__1.r = -tau->r, z__1.i = -tau->i;
            zgerc_(&i__1, &lastc, &z__1, &v[i__], incv, &work[1], &c__1, &c__[firstv + c_dim1],
                   ldc);
        }
    } else {
        if (lastv == firstv) {
            z__1.r = 1. - tau->r, z__1.i = 0. - tau->i;
            zscal_(&lastc, &z__1, &c__[lastv * c_dim1 + 1], &c__1);
        } else {
            i__1 = lastv - firstv;
            zgemv_((char *)"N", &lastc, &i__1, &c_b1, &c__[firstv * c_dim1 + 1], ldc, &v[i__], incv, &c_b2,
                   &work[1], &c__1, (ftnlen)1);
            zaxpy_(&lastc, &c_b1, &c__[lastv * c_dim1 + 1], &c__1, &work[1], &c__1);
            z__1.r = -tau->r, z__1.i = -tau->i;
            zaxpy_(&lastc, &z__1, &work[1], &c__1, &c__[lastv * c_dim1 + 1], &c__1);
            i__1 = lastv - firstv;
            z__1.r = -tau->r, z__1.i = -tau->i;
            zgerc_(&lastc, &i__1, &z__1, &work[1], &c__1, &v[i__], incv, &c__[firstv * c_dim1 + 1],
                   ldc);
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
