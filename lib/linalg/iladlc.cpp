#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
integer iladlc_(integer *m, integer *n, doublereal *a, integer *lda)
{
    integer a_dim1, a_offset, ret_val, i__1;
    integer i__;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    if (*n == 0) {
        ret_val = *n;
    } else if (a[*n * a_dim1 + 1] != 0. || a[*m + *n * a_dim1] != 0.) {
        ret_val = *n;
    } else {
        for (ret_val = *n; ret_val >= 1; --ret_val) {
            i__1 = *m;
            for (i__ = 1; i__ <= i__1; ++i__) {
                if (a[i__ + ret_val * a_dim1] != 0.) {
                    return ret_val;
                }
            }
        }
    }
    return ret_val;
}
#ifdef __cplusplus
}
#endif
