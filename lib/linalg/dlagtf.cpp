#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlagtf_(integer *n, doublereal *a, doublereal *lambda, doublereal *b, doublereal *c__,
            doublereal *tol, doublereal *d__, integer *in, integer *info)
{
    integer i__1;
    doublereal d__1, d__2;
    integer k;
    doublereal tl, eps, piv1, piv2, temp, mult, scale1, scale2;
    extern doublereal dlamch_(char *, ftnlen);
    extern int xerbla_(char *, integer *, ftnlen);
    --in;
    --d__;
    --c__;
    --b;
    --a;
    *info = 0;
    if (*n < 0) {
        *info = -1;
        i__1 = -(*info);
        xerbla_((char *)"DLAGTF", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    a[1] -= *lambda;
    in[*n] = 0;
    if (*n == 1) {
        if (a[1] == 0.) {
            in[1] = 1;
        }
        return 0;
    }
    eps = dlamch_((char *)"Epsilon", (ftnlen)7);
    tl = max(*tol, eps);
    scale1 = abs(a[1]) + abs(b[1]);
    i__1 = *n - 1;
    for (k = 1; k <= i__1; ++k) {
        a[k + 1] -= *lambda;
        scale2 = (d__1 = c__[k], abs(d__1)) + (d__2 = a[k + 1], abs(d__2));
        if (k < *n - 1) {
            scale2 += (d__1 = b[k + 1], abs(d__1));
        }
        if (a[k] == 0.) {
            piv1 = 0.;
        } else {
            piv1 = (d__1 = a[k], abs(d__1)) / scale1;
        }
        if (c__[k] == 0.) {
            in[k] = 0;
            piv2 = 0.;
            scale1 = scale2;
            if (k < *n - 1) {
                d__[k] = 0.;
            }
        } else {
            piv2 = (d__1 = c__[k], abs(d__1)) / scale2;
            if (piv2 <= piv1) {
                in[k] = 0;
                scale1 = scale2;
                c__[k] /= a[k];
                a[k + 1] -= c__[k] * b[k];
                if (k < *n - 1) {
                    d__[k] = 0.;
                }
            } else {
                in[k] = 1;
                mult = a[k] / c__[k];
                a[k] = c__[k];
                temp = a[k + 1];
                a[k + 1] = b[k] - mult * temp;
                if (k < *n - 1) {
                    d__[k] = b[k + 1];
                    b[k + 1] = -mult * d__[k];
                }
                b[k] = temp;
                c__[k] = mult;
            }
        }
        if (max(piv1, piv2) <= tl && in[*n] == 0) {
            in[*n] = k;
        }
    }
    if ((d__1 = a[*n], abs(d__1)) <= scale1 * tl && in[*n] == 0) {
        in[*n] = *n;
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
