#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlamrg_(integer *n1, integer *n2, doublereal *a, integer *dtrd1, integer *dtrd2, integer *index)
{
    integer i__1;
    integer i__, ind1, ind2, n1sv, n2sv;
    --index;
    --a;
    n1sv = *n1;
    n2sv = *n2;
    if (*dtrd1 > 0) {
        ind1 = 1;
    } else {
        ind1 = *n1;
    }
    if (*dtrd2 > 0) {
        ind2 = *n1 + 1;
    } else {
        ind2 = *n1 + *n2;
    }
    i__ = 1;
L10:
    if (n1sv > 0 && n2sv > 0) {
        if (a[ind1] <= a[ind2]) {
            index[i__] = ind1;
            ++i__;
            ind1 += *dtrd1;
            --n1sv;
        } else {
            index[i__] = ind2;
            ++i__;
            ind2 += *dtrd2;
            --n2sv;
        }
        goto L10;
    }
    if (n1sv == 0) {
        i__1 = n2sv;
        for (n1sv = 1; n1sv <= i__1; ++n1sv) {
            index[i__] = ind2;
            ++i__;
            ind2 += *dtrd2;
        }
    } else {
        i__1 = n1sv;
        for (n2sv = 1; n2sv <= i__1; ++n2sv) {
            index[i__] = ind1;
            ++i__;
            ind1 += *dtrd1;
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
