#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b4 = 1.;
static doublereal c_b5 = 0.;
static integer c__1 = 1;
int dlarf_(char *side, integer *m, integer *n, doublereal *v, integer *incv, doublereal *tau,
           doublereal *c__, integer *ldc, doublereal *work, ftnlen side_len)
{
    integer c_dim1, c_offset;
    doublereal d__1;
    integer i__;
    logical applyleft;
    extern int dger_(integer *, integer *, doublereal *, doublereal *, integer *, doublereal *,
                     integer *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *,
                      doublereal *, integer *, doublereal *, doublereal *, integer *, ftnlen);
    integer lastc, lastv;
    extern integer iladlc_(integer *, integer *, doublereal *, integer *),
        iladlr_(integer *, integer *, doublereal *, integer *);
    --v;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    applyleft = lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1);
    lastv = 0;
    lastc = 0;
    if (*tau != 0.) {
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
        while (lastv > 0 && v[i__] == 0.) {
            --lastv;
            i__ -= *incv;
        }
        if (applyleft) {
            lastc = iladlc_(&lastv, n, &c__[c_offset], ldc);
        } else {
            lastc = iladlr_(m, &lastv, &c__[c_offset], ldc);
        }
    }
    if (applyleft) {
        if (lastv > 0) {
            dgemv_((char *)"Transpose", &lastv, &lastc, &c_b4, &c__[c_offset], ldc, &v[1], incv, &c_b5,
                   &work[1], &c__1, (ftnlen)9);
            d__1 = -(*tau);
            dger_(&lastv, &lastc, &d__1, &v[1], incv, &work[1], &c__1, &c__[c_offset], ldc);
        }
    } else {
        if (lastv > 0) {
            dgemv_((char *)"No transpose", &lastc, &lastv, &c_b4, &c__[c_offset], ldc, &v[1], incv, &c_b5,
                   &work[1], &c__1, (ftnlen)12);
            d__1 = -(*tau);
            dger_(&lastc, &lastv, &d__1, &work[1], &c__1, &v[1], incv, &c__[c_offset], ldc);
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
