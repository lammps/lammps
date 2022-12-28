#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
integer ilazlc_(integer *m, integer *n, doublecomplex *a, integer *lda)
{
    integer a_dim1, a_offset, ret_val, i__1, i__2;
    integer i__;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    if (*n == 0) {
        ret_val = *n;
    } else {
        i__1 = *n * a_dim1 + 1;
        i__2 = *m + *n * a_dim1;
        if (a[i__1].r != 0. || a[i__1].i != 0. || (a[i__2].r != 0. || a[i__2].i != 0.)) {
            ret_val = *n;
        } else {
            for (ret_val = *n; ret_val >= 1; --ret_val) {
                i__1 = *m;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    i__2 = i__ + ret_val * a_dim1;
                    if (a[i__2].r != 0. || a[i__2].i != 0.) {
                        return ret_val;
                    }
                }
            }
        }
    }
    return ret_val;
}
#ifdef __cplusplus
}
#endif
