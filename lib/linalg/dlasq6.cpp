#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlasq6_(integer *i0, integer *n0, doublereal *z__, integer *pp, doublereal *dmin__,
            doublereal *dmin1, doublereal *dmin2, doublereal *dn, doublereal *dnm1,
            doublereal *dnm2)
{
    integer i__1;
    doublereal d__1, d__2;
    doublereal d__;
    integer j4, j4p2;
    doublereal emin, temp;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal safmin;
    --z__;
    if (*n0 - *i0 - 1 <= 0) {
        return 0;
    }
    safmin = dlamch_((char *)"Safe minimum", (ftnlen)12);
    j4 = (*i0 << 2) + *pp - 3;
    emin = z__[j4 + 4];
    d__ = z__[j4];
    *dmin__ = d__;
    if (*pp == 0) {
        i__1 = *n0 - 3 << 2;
        for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
            z__[j4 - 2] = d__ + z__[j4 - 1];
            if (z__[j4 - 2] == 0.) {
                z__[j4] = 0.;
                d__ = z__[j4 + 1];
                *dmin__ = d__;
                emin = 0.;
            } else if (safmin * z__[j4 + 1] < z__[j4 - 2] && safmin * z__[j4 - 2] < z__[j4 + 1]) {
                temp = z__[j4 + 1] / z__[j4 - 2];
                z__[j4] = z__[j4 - 1] * temp;
                d__ *= temp;
            } else {
                z__[j4] = z__[j4 + 1] * (z__[j4 - 1] / z__[j4 - 2]);
                d__ = z__[j4 + 1] * (d__ / z__[j4 - 2]);
            }
            *dmin__ = min(*dmin__, d__);
            d__1 = emin, d__2 = z__[j4];
            emin = min(d__1, d__2);
        }
    } else {
        i__1 = *n0 - 3 << 2;
        for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
            z__[j4 - 3] = d__ + z__[j4];
            if (z__[j4 - 3] == 0.) {
                z__[j4 - 1] = 0.;
                d__ = z__[j4 + 2];
                *dmin__ = d__;
                emin = 0.;
            } else if (safmin * z__[j4 + 2] < z__[j4 - 3] && safmin * z__[j4 - 3] < z__[j4 + 2]) {
                temp = z__[j4 + 2] / z__[j4 - 3];
                z__[j4 - 1] = z__[j4] * temp;
                d__ *= temp;
            } else {
                z__[j4 - 1] = z__[j4 + 2] * (z__[j4] / z__[j4 - 3]);
                d__ = z__[j4 + 2] * (d__ / z__[j4 - 3]);
            }
            *dmin__ = min(*dmin__, d__);
            d__1 = emin, d__2 = z__[j4 - 1];
            emin = min(d__1, d__2);
        }
    }
    *dnm2 = d__;
    *dmin2 = *dmin__;
    j4 = (*n0 - 2 << 2) - *pp;
    j4p2 = j4 + (*pp << 1) - 1;
    z__[j4 - 2] = *dnm2 + z__[j4p2];
    if (z__[j4 - 2] == 0.) {
        z__[j4] = 0.;
        *dnm1 = z__[j4p2 + 2];
        *dmin__ = *dnm1;
        emin = 0.;
    } else if (safmin * z__[j4p2 + 2] < z__[j4 - 2] && safmin * z__[j4 - 2] < z__[j4p2 + 2]) {
        temp = z__[j4p2 + 2] / z__[j4 - 2];
        z__[j4] = z__[j4p2] * temp;
        *dnm1 = *dnm2 * temp;
    } else {
        z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
        *dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]);
    }
    *dmin__ = min(*dmin__, *dnm1);
    *dmin1 = *dmin__;
    j4 += 4;
    j4p2 = j4 + (*pp << 1) - 1;
    z__[j4 - 2] = *dnm1 + z__[j4p2];
    if (z__[j4 - 2] == 0.) {
        z__[j4] = 0.;
        *dn = z__[j4p2 + 2];
        *dmin__ = *dn;
        emin = 0.;
    } else if (safmin * z__[j4p2 + 2] < z__[j4 - 2] && safmin * z__[j4 - 2] < z__[j4p2 + 2]) {
        temp = z__[j4p2 + 2] / z__[j4 - 2];
        z__[j4] = z__[j4p2] * temp;
        *dn = *dnm1 * temp;
    } else {
        z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
        *dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]);
    }
    *dmin__ = min(*dmin__, *dn);
    z__[j4 + 2] = *dn;
    z__[(*n0 << 2) - *pp] = emin;
    return 0;
}
#ifdef __cplusplus
}
#endif
