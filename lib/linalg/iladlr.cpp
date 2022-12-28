#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
integer iladlr_(integer *m, integer *n, doublereal *a, integer *lda)
{
    integer a_dim1, a_offset, ret_val, i__1;
    integer i__, j;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    if (*m == 0) {
        ret_val = *m;
    } else if (a[*m + a_dim1] != 0. || a[*m + *n * a_dim1] != 0.) {
        ret_val = *m;
    } else {
        ret_val = 0;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__ = *m;
            while (a[max(i__, 1) + j * a_dim1] == 0. && i__ >= 1) {
                --i__;
            }
            ret_val = max(ret_val, i__);
        }
    }
    return ret_val;
}
#ifdef __cplusplus
}
#endif
