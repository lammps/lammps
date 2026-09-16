#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlarrr_(integer *n, doublereal *d__, doublereal *e, integer *info)
{
    integer i__1;
    doublereal d__1;
    double sqrt(doublereal);
    integer i__;
    doublereal eps, tmp, tmp2, rmin;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal offdig, safmin;
    logical yesrel;
    doublereal smlnum, offdig2;
    --e;
    --d__;
    if (*n <= 0) {
        *info = 0;
        return 0;
    }
    *info = 1;
    safmin = dlamch_((char *)"Safe minimum", (ftnlen)12);
    eps = dlamch_((char *)"Precision", (ftnlen)9);
    smlnum = safmin / eps;
    rmin = sqrt(smlnum);
    yesrel = TRUE_;
    offdig = 0.;
    tmp = sqrt((abs(d__[1])));
    if (tmp < rmin) {
        yesrel = FALSE_;
    }
    if (!yesrel) {
        goto L11;
    }
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
        tmp2 = sqrt((d__1 = d__[i__], abs(d__1)));
        if (tmp2 < rmin) {
            yesrel = FALSE_;
        }
        if (!yesrel) {
            goto L11;
        }
        offdig2 = (d__1 = e[i__ - 1], abs(d__1)) / (tmp * tmp2);
        if (offdig + offdig2 >= .999) {
            yesrel = FALSE_;
        }
        if (!yesrel) {
            goto L11;
        }
        tmp = tmp2;
        offdig = offdig2;
    }
L11:
    if (yesrel) {
        *info = 0;
        return 0;
    } else {
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
