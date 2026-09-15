
#include "lmp_f2c.h"

extern "C" {

double pow_lmp_di(doublereal *ap, integer *bp)
{
    double pow, x;
    integer n;
    unsigned long u;

    pow = 1.0;
    x = *ap;
    n = *bp;

    if (n != 0) {
        if (x == 0.0) return 0.0;
        if (n < 0) {
            n = -n;
            x = 1.0 / x;
        }
        for (u = n;;) {
            if (u & 01) pow *= x;
            if (u >>= 1)
                x *= x;
            else
                break;
        }
    }
    return pow;
}
}
