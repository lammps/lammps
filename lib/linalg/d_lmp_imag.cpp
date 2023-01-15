
#include "lmp_f2c.h"

extern "C" {

double d_lmp_imag(doublecomplex *z)
{
    return (z->i);
}
}
