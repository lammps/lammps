
#include "lmp_f2c.h"

extern "C" {

/* assign strings:  a = b */

void s_lmp_copy(char *a, char *b, ftnlen la, ftnlen lb)
{
    char *aend, *bend;

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
