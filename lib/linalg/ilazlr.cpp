#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
integer ilazlr_(integer *m, integer *n, doublecomplex *a, integer *lda)
{
    integer a_dim1, a_offset, ret_val, i__1, i__2;
    integer i__, j;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    if (*m == 0) {
        ret_val = *m;
    } else {
        i__1 = *m + a_dim1;
        i__2 = *m + *n * a_dim1;
        if (a[i__1].r != 0. || a[i__1].i != 0. || (a[i__2].r != 0. || a[i__2].i != 0.)) {
            ret_val = *m;
        } else {
            ret_val = 0;
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__ = *m;
                for (;;) {
                    i__2 = max(i__, 1) + j * a_dim1;
                    if (!(a[i__2].r == 0. && a[i__2].i == 0. && i__ >= 1)) break;
                    --i__;
                }
                ret_val = max(ret_val, i__);
            }
        }
    }
    return ret_val;
}
#ifdef __cplusplus
}
#endif
