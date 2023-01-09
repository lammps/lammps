
#include "lmp_f2c.h"

// concatenate two strings

extern "C" {
void s_lmp_cat(char *lp, char *rpp[], ftnint rnp[], ftnint *np, ftnlen ll)
{
    ftnlen i, nc;
    char *rp;
    ftnlen n = *np;
    for (i = 0; i < n; ++i) {
        nc = ll;
        if (rnp[i] < nc) nc = rnp[i];
        ll -= nc;
        rp = rpp[i];
        while (--nc >= 0)
            *lp++ = *rp++;
    }
    while (--ll >= 0)
        *lp++ = ' ';
}
}
