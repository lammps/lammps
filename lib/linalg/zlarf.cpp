#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublecomplex c_b1 = {1., 0.};
static doublecomplex c_b2 = {0., 0.};
static integer c__1 = 1;
int zlarf_(char *side, integer *m, integer *n, doublecomplex *v, integer *incv, doublecomplex *tau,
           doublecomplex *c__, integer *ldc, doublecomplex *work, ftnlen side_len)
{
    integer c_dim1, c_offset, i__1;
    doublecomplex z__1;
    integer i__;
    logical applyleft;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer lastc;
    extern int zgerc_(integer *, integer *, doublecomplex *, doublecomplex *, integer *,
                      doublecomplex *, integer *, doublecomplex *, integer *),
        zgemv_(char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *,
               doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    integer lastv;
    extern integer ilazlc_(integer *, integer *, doublecomplex *, integer *),
        ilazlr_(integer *, integer *, doublecomplex *, integer *);
    --v;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    applyleft = lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1);
    lastv = 0;
    lastc = 0;
    if (tau->r != 0. || tau->i != 0.) {
        if (applyleft) {
            lastv = *m;
        } else {
            lastv = *n;
        }
        if (*incv > 0) {
            i__ = (lastv - 1) * *incv + 1;
        } else {
            i__ = 1;
        }
        for (;;) {
            i__1 = i__;
            if (!(lastv > 0 && (v[i__1].r == 0. && v[i__1].i == 0.))) break;
            --lastv;
            i__ -= *incv;
        }
        if (applyleft) {
            lastc = ilazlc_(&lastv, n, &c__[c_offset], ldc);
        } else {
            lastc = ilazlr_(m, &lastv, &c__[c_offset], ldc);
        }
    }
    if (applyleft) {
        if (lastv > 0) {
            zgemv_((char *)"Conjugate transpose", &lastv, &lastc, &c_b1, &c__[c_offset], ldc, &v[1], incv,
                   &c_b2, &work[1], &c__1, (ftnlen)19);
            z__1.r = -tau->r, z__1.i = -tau->i;
            zgerc_(&lastv, &lastc, &z__1, &v[1], incv, &work[1], &c__1, &c__[c_offset], ldc);
        }
    } else {
        if (lastv > 0) {
            zgemv_((char *)"No transpose", &lastc, &lastv, &c_b1, &c__[c_offset], ldc, &v[1], incv, &c_b2,
                   &work[1], &c__1, (ftnlen)12);
            z__1.r = -tau->r, z__1.i = -tau->i;
            zgerc_(&lastc, &lastv, &z__1, &work[1], &c__1, &v[1], incv, &c__[c_offset], ldc);
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
