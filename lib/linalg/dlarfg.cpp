#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlarfg_(integer *n, doublereal *alpha, doublereal *x, integer *incx, doublereal *tau)
{
    integer i__1;
    doublereal d__1;
    double d_lmp_sign(doublereal *, doublereal *);
    integer j, knt;
    doublereal beta;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern int dscal_(integer *, doublereal *, doublereal *, integer *);
    doublereal xnorm;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, ftnlen);
    doublereal safmin, rsafmn;
    --x;
    if (*n <= 1) {
        *tau = 0.;
        return 0;
    }
    i__1 = *n - 1;
    xnorm = dnrm2_(&i__1, &x[1], incx);
    if (xnorm == 0.) {
        *tau = 0.;
    } else {
        d__1 = dlapy2_(alpha, &xnorm);
        beta = -d_lmp_sign(&d__1, alpha);
        safmin = dlamch_((char *)"S", (ftnlen)1) / dlamch_((char *)"E", (ftnlen)1);
        knt = 0;
        if (abs(beta) < safmin) {
            rsafmn = 1. / safmin;
        L10:
            ++knt;
            i__1 = *n - 1;
            dscal_(&i__1, &rsafmn, &x[1], incx);
            beta *= rsafmn;
            *alpha *= rsafmn;
            if (abs(beta) < safmin && knt < 20) {
                goto L10;
            }
            i__1 = *n - 1;
            xnorm = dnrm2_(&i__1, &x[1], incx);
            d__1 = dlapy2_(alpha, &xnorm);
            beta = -d_lmp_sign(&d__1, alpha);
        }
        *tau = (beta - *alpha) / beta;
        i__1 = *n - 1;
        d__1 = 1. / (*alpha - beta);
        dscal_(&i__1, &d__1, &x[1], incx);
        i__1 = knt;
        for (j = 1; j <= i__1; ++j) {
            beta *= safmin;
        }
        *alpha = beta;
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
