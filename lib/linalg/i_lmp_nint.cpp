
#include "lmp_f2c.h"
#undef abs

#include <cmath>

extern "C" {

integer i_lmp_nint(real *x)
{
    return (integer)(*x >= 0 ? floor(*x + .5) : -floor(.5 - *x));
}
}
