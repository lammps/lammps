
#include "lmp_f2c.h"

extern "C" {

double d_lmp_sign(doublereal *a, doublereal *b)
{
    double x;
    x = (*a >= 0 ? *a : -*a);
    return (*b >= 0 ? x : -x);
}
}
