
#include "lmp_f2c.h"

extern "C" {

void d_lmp_cnjg(doublecomplex *r, doublecomplex *z)
{
    doublereal zi = z->i;

    r->r = z->r;
    r->i = -zi;
}
}
