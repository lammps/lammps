#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b4 = 1.;
static doublereal c_b5 = 0.;
static integer c__1 = 1;
int dlarf1l_(char *side, integer *m, integer *n, doublereal *v, integer *incv, doublereal *tau,
             doublereal *c__, integer *ldc, doublereal *work, ftnlen side_len)
{
    integer c_dim1, c_offset, i__1;
    doublereal d__1;
    integer i__;
    logical applyleft;
    extern int dger_(integer *, integer *, doublereal *, doublereal *, integer *, doublereal *,
                     integer *, doublereal *, integer *),
        dscal_(integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *,
                      doublereal *, integer *, doublereal *, doublereal *, integer *, ftnlen);
    integer lastc;
    extern int daxpy_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *);
    integer lastv;
    extern integer iladlc_(integer *, integer *, doublereal *, integer *),
        iladlr_(integer *, integer *, doublereal *, integer *);
    integer firstv;
    --v;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    applyleft = lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1);
    firstv = 1;
    lastc = 0;
    if (*tau != 0.) {
        if (applyleft) {
            lastv = *m;
        } else {
            lastv = *n;
        }
        i__ = 1;
        while (lastv > firstv && v[i__] == 0.) {
            ++firstv;
            i__ += *incv;
        }
        if (applyleft) {
            lastc = iladlc_(&lastv, n, &c__[c_offset], ldc);
        } else {
            lastc = iladlr_(m, &lastv, &c__[c_offset], ldc);
        }
    }
    if (lastc == 0) {
        return 0;
    }
    if (applyleft) {
        if (lastv > 0) {
            if (lastv == firstv) {
                d__1 = 1. - *tau;
                dscal_(&lastc, &d__1, &c__[firstv + c_dim1], ldc);
            } else {
                i__1 = lastv - firstv;
                dgemv_((char *)"T", &i__1, &lastc, &c_b4, &c__[firstv + c_dim1], ldc, &v[i__], incv, &c_b5,
                       &work[1], &c__1, (ftnlen)1);
                daxpy_(&lastc, &c_b4, &c__[lastv + c_dim1], ldc, &work[1], &c__1);
                d__1 = -(*tau);
                daxpy_(&lastc, &d__1, &work[1], &c__1, &c__[lastv + c_dim1], ldc);
                i__1 = lastv - firstv;
                d__1 = -(*tau);
                dger_(&i__1, &lastc, &d__1, &v[i__], incv, &work[1], &c__1, &c__[firstv + c_dim1],
                      ldc);
            }
        }
    } else {
        if (lastv > 0) {
            if (lastv == firstv) {
                d__1 = 1. - *tau;
                dscal_(&lastc, &d__1, &c__[c_offset], &c__1);
            } else {
                i__1 = lastv - firstv;
                dgemv_((char *)"N", &lastc, &i__1, &c_b4, &c__[firstv * c_dim1 + 1], ldc, &v[i__], incv,
                       &c_b5, &work[1], &c__1, (ftnlen)1);
                daxpy_(&lastc, &c_b4, &c__[lastv * c_dim1 + 1], &c__1, &work[1], &c__1);
                d__1 = -(*tau);
                daxpy_(&lastc, &d__1, &work[1], &c__1, &c__[lastv * c_dim1 + 1], &c__1);
                i__1 = lastv - firstv;
                d__1 = -(*tau);
                dger_(&lastc, &i__1, &d__1, &work[1], &c__1, &v[i__], incv,
                      &c__[firstv * c_dim1 + 1], ldc);
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
