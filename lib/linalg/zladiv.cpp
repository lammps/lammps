#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
VOID zladiv_(doublecomplex *ret_val, doublecomplex *x, doublecomplex *y)
{
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1;
    double d_lmp_imag(doublecomplex *);
    doublereal zi, zr;
    extern int dladiv_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                       doublereal *);
    d__1 = x->r;
    d__2 = d_lmp_imag(x);
    d__3 = y->r;
    d__4 = d_lmp_imag(y);
    dladiv_(&d__1, &d__2, &d__3, &d__4, &zr, &zi);
    z__1.r = zr, z__1.i = zi;
    ret_val->r = z__1.r, ret_val->i = z__1.i;
    return;
}
#ifdef __cplusplus
}
#endif
