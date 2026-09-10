#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
doublereal droundup_lwork__(integer *lwork)
{
    doublereal ret_val;
    doublereal eps;
    extern doublereal dlamch_(char *, ftnlen);
    ret_val = (doublereal)(*lwork);
    if ((integer)ret_val < *lwork) {
        eps = dlamch_((char *)"E", (ftnlen)1);
        ret_val *= eps + 1.;
    }
    return ret_val;
}
#ifdef __cplusplus
}
#endif
