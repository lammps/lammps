#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int drscl_(integer *n, doublereal *sa, doublereal *sx, integer *incx)
{
    doublereal mul, cden;
    logical done;
    doublereal cnum, cden1, cnum1;
    extern int dscal_(integer *, doublereal *, doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    doublereal bignum, smlnum;
    --sx;
    if (*n <= 0) {
        return 0;
    }
    smlnum = dlamch_((char *)"S", (ftnlen)1);
    bignum = 1. / smlnum;
    cden = *sa;
    cnum = 1.;
L10:
    cden1 = cden * smlnum;
    cnum1 = cnum / bignum;
    if (abs(cden1) > abs(cnum) && cnum != 0.) {
        mul = smlnum;
        done = FALSE_;
        cden = cden1;
    } else if (abs(cnum1) > abs(cden)) {
        mul = bignum;
        done = FALSE_;
        cnum = cnum1;
    } else {
        mul = cnum / cden;
        done = TRUE_;
    }
    dscal_(n, &mul, &sx[1], incx);
    if (!done) {
        goto L10;
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
