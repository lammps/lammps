
#include "lmp_f2c.h"
#undef abs

#include <cmath>

extern "C" {
double pow_lmp_dd(doublereal *ap, doublereal *bp)
{
    return (pow(*ap, *bp));
}
}
