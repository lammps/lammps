
#include "lmp_f2c.h"

extern "C" {

/* assign strings:  a = b */

void s_lmp_copy(register char *a, register char *b, ftnlen la, ftnlen lb)
{
    register char *aend, *bend;

    aend = a + la;

    if (la <= lb)
        while (a < aend)
            *a++ = *b++;

    else {
        bend = b + lb;
        while (b < bend)
            *a++ = *b++;
        while (a < aend)
            *a++ = ' ';
    }
}
}
