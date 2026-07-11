
#include "lmp_f2c.h"
#undef abs

static constexpr double log10e = 0.43429448190325182765;

#include <cmath>

extern "C" {
double d_lmp_lg10(doublereal *x)
{
    return (log10e * log(*x));
}
}
