#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
integer ieeeck_(integer *ispec, real *zero, real *one)
{
    integer ret_val;
    real nan1, nan2, nan3, nan4, nan5, nan6, neginf, posinf, negzro, newzro;
    ret_val = 1;
    posinf = *one / *zero;
    if (posinf <= *one) {
        ret_val = 0;
        return ret_val;
    }
    neginf = -(*one) / *zero;
    if (neginf >= *zero) {
        ret_val = 0;
        return ret_val;
    }
    negzro = *one / (neginf + *one);
    if (negzro != *zero) {
        ret_val = 0;
        return ret_val;
    }
    neginf = *one / negzro;
    if (neginf >= *zero) {
        ret_val = 0;
        return ret_val;
    }
    newzro = negzro + *zero;
    if (newzro != *zero) {
        ret_val = 0;
        return ret_val;
    }
    posinf = *one / newzro;
    if (posinf <= *one) {
        ret_val = 0;
        return ret_val;
    }
    neginf *= posinf;
    if (neginf >= *zero) {
        ret_val = 0;
        return ret_val;
    }
    posinf *= posinf;
    if (posinf <= *one) {
        ret_val = 0;
        return ret_val;
    }
    if (*ispec == 0) {
        return ret_val;
    }
    nan1 = posinf + neginf;
    nan2 = posinf / neginf;
    nan3 = posinf / posinf;
    nan4 = posinf * *zero;
    nan5 = neginf * negzro;
    nan6 = nan5 * *zero;
    if (nan1 == nan1) {
        ret_val = 0;
        return ret_val;
    }
    if (nan2 == nan2) {
        ret_val = 0;
        return ret_val;
    }
    if (nan3 == nan3) {
        ret_val = 0;
        return ret_val;
    }
    if (nan4 == nan4) {
        ret_val = 0;
        return ret_val;
    }
    if (nan5 == nan5) {
        ret_val = 0;
        return ret_val;
    }
    if (nan6 == nan6) {
        ret_val = 0;
        return ret_val;
    }
    return ret_val;
}
#ifdef __cplusplus
}
#endif
