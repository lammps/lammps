#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
integer izamax_(integer *n, doublecomplex *zx, integer *incx)
{
    integer ret_val, i__1;
    integer i__, ix;
    doublereal dmax__;
    extern doublereal dcabs1_(doublecomplex *);
    --zx;
    ret_val = 0;
    if (*n < 1 || *incx <= 0) {
        return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
        return ret_val;
    }
    if (*incx == 1) {
        dmax__ = dcabs1_(&zx[1]);
        i__1 = *n;
        for (i__ = 2; i__ <= i__1; ++i__) {
            if (dcabs1_(&zx[i__]) > dmax__) {
                ret_val = i__;
                dmax__ = dcabs1_(&zx[i__]);
            }
        }
    } else {
        ix = 1;
        dmax__ = dcabs1_(&zx[1]);
        ix += *incx;
        i__1 = *n;
        for (i__ = 2; i__ <= i__1; ++i__) {
            if (dcabs1_(&zx[ix]) > dmax__) {
                ret_val = i__;
                dmax__ = dcabs1_(&zx[ix]);
            }
            ix += *incx;
        }
    }
    return ret_val;
}
#ifdef __cplusplus
}
#endif
