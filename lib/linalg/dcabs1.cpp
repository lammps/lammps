#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
doublereal dcabs1_(doublecomplex *z__)
{
    doublereal ret_val, d__1, d__2;
    double d_lmp_imag(doublecomplex *);
    ret_val = (d__1 = z__->r, abs(d__1)) + (d__2 = d_lmp_imag(z__), abs(d__2));
    return ret_val;
}
#ifdef __cplusplus
}
#endif
