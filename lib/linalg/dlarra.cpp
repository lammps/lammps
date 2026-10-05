#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlarra_(integer *n, doublereal *d__, doublereal *e, doublereal *e2, doublereal *spltol,
            doublereal *tnrm, integer *nsplit, integer *isplit, integer *info)
{
    integer i__1;
    doublereal d__1, d__2;
    double sqrt(doublereal);
    integer i__;
    doublereal tmp1, eabs;
    --isplit;
    --e2;
    --e;
    --d__;
    *info = 0;
    *nsplit = 1;
    if (*n <= 0) {
        return 0;
    }
    if (*spltol < 0.) {
        tmp1 = abs(*spltol) * *tnrm;
        i__1 = *n - 1;
        for (i__ = 1; i__ <= i__1; ++i__) {
            eabs = (d__1 = e[i__], abs(d__1));
            if (eabs <= tmp1) {
                e[i__] = 0.;
                e2[i__] = 0.;
                isplit[*nsplit] = i__;
                ++(*nsplit);
            }
        }
    } else {
        i__1 = *n - 1;
        for (i__ = 1; i__ <= i__1; ++i__) {
            eabs = (d__1 = e[i__], abs(d__1));
            if (eabs <= *spltol * sqrt((d__1 = d__[i__], abs(d__1))) *
                            sqrt((d__2 = d__[i__ + 1], abs(d__2)))) {
                e[i__] = 0.;
                e2[i__] = 0.;
                isplit[*nsplit] = i__;
                ++(*nsplit);
            }
        }
    }
    isplit[*nsplit] = *n;
    return 0;
}
#ifdef __cplusplus
}
#endif
