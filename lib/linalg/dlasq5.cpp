#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlasq5_(integer *i0, integer *n0, doublereal *z__, integer *pp, doublereal *tau,
            doublereal *sigma, doublereal *dmin__, doublereal *dmin1, doublereal *dmin2,
            doublereal *dn, doublereal *dnm1, doublereal *dnm2, logical *ieee, doublereal *eps)
{
    integer i__1;
    doublereal d__1, d__2;
    doublereal d__;
    integer j4, j4p2;
    doublereal emin, temp, dthresh;
    --z__;
    if (*n0 - *i0 - 1 <= 0) {
        return 0;
    }
    dthresh = *eps * (*sigma + *tau);
    if (*tau < dthresh * .5) {
        *tau = 0.;
    }
    if (*tau != 0.) {
        j4 = (*i0 << 2) + *pp - 3;
        emin = z__[j4 + 4];
        d__ = z__[j4] - *tau;
        *dmin__ = d__;
        *dmin1 = -z__[j4];
        if (*ieee) {
            if (*pp == 0) {
                i__1 = *n0 - 3 << 2;
                for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
                    z__[j4 - 2] = d__ + z__[j4 - 1];
                    temp = z__[j4 + 1] / z__[j4 - 2];
                    d__ = d__ * temp - *tau;
                    *dmin__ = min(*dmin__, d__);
                    z__[j4] = z__[j4 - 1] * temp;
                    d__1 = z__[j4];
                    emin = min(d__1, emin);
                }
            } else {
                i__1 = *n0 - 3 << 2;
                for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
                    z__[j4 - 3] = d__ + z__[j4];
                    temp = z__[j4 + 2] / z__[j4 - 3];
                    d__ = d__ * temp - *tau;
                    *dmin__ = min(*dmin__, d__);
                    z__[j4 - 1] = z__[j4] * temp;
                    d__1 = z__[j4 - 1];
                    emin = min(d__1, emin);
                }
            }
            *dnm2 = d__;
            *dmin2 = *dmin__;
            j4 = (*n0 - 2 << 2) - *pp;
            j4p2 = j4 + (*pp << 1) - 1;
            z__[j4 - 2] = *dnm2 + z__[j4p2];
            z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
            *dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
            *dmin__ = min(*dmin__, *dnm1);
            *dmin1 = *dmin__;
            j4 += 4;
            j4p2 = j4 + (*pp << 1) - 1;
            z__[j4 - 2] = *dnm1 + z__[j4p2];
            z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
            *dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
            *dmin__ = min(*dmin__, *dn);
        } else {
            if (*pp == 0) {
                i__1 = *n0 - 3 << 2;
                for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
                    z__[j4 - 2] = d__ + z__[j4 - 1];
                    if (d__ < 0.) {
                        return 0;
                    } else {
                        z__[j4] = z__[j4 + 1] * (z__[j4 - 1] / z__[j4 - 2]);
                        d__ = z__[j4 + 1] * (d__ / z__[j4 - 2]) - *tau;
                    }
                    *dmin__ = min(*dmin__, d__);
                    d__1 = emin, d__2 = z__[j4];
                    emin = min(d__1, d__2);
                }
            } else {
                i__1 = *n0 - 3 << 2;
                for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
                    z__[j4 - 3] = d__ + z__[j4];
                    if (d__ < 0.) {
                        return 0;
                    } else {
                        z__[j4 - 1] = z__[j4 + 2] * (z__[j4] / z__[j4 - 3]);
                        d__ = z__[j4 + 2] * (d__ / z__[j4 - 3]) - *tau;
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
            if (*dnm2 < 0.) {
                return 0;
            } else {
                z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
                *dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
            }
            *dmin__ = min(*dmin__, *dnm1);
            *dmin1 = *dmin__;
            j4 += 4;
            j4p2 = j4 + (*pp << 1) - 1;
            z__[j4 - 2] = *dnm1 + z__[j4p2];
            if (*dnm1 < 0.) {
                return 0;
            } else {
                z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
                *dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
            }
            *dmin__ = min(*dmin__, *dn);
        }
    } else {
        j4 = (*i0 << 2) + *pp - 3;
        emin = z__[j4 + 4];
        d__ = z__[j4] - *tau;
        *dmin__ = d__;
        *dmin1 = -z__[j4];
        if (*ieee) {
            if (*pp == 0) {
                i__1 = *n0 - 3 << 2;
                for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
                    z__[j4 - 2] = d__ + z__[j4 - 1];
                    temp = z__[j4 + 1] / z__[j4 - 2];
                    d__ = d__ * temp - *tau;
                    if (d__ < dthresh) {
                        d__ = 0.;
                    }
                    *dmin__ = min(*dmin__, d__);
                    z__[j4] = z__[j4 - 1] * temp;
                    d__1 = z__[j4];
                    emin = min(d__1, emin);
                }
            } else {
                i__1 = *n0 - 3 << 2;
                for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
                    z__[j4 - 3] = d__ + z__[j4];
                    temp = z__[j4 + 2] / z__[j4 - 3];
                    d__ = d__ * temp - *tau;
                    if (d__ < dthresh) {
                        d__ = 0.;
                    }
                    *dmin__ = min(*dmin__, d__);
                    z__[j4 - 1] = z__[j4] * temp;
                    d__1 = z__[j4 - 1];
                    emin = min(d__1, emin);
                }
            }
            *dnm2 = d__;
            *dmin2 = *dmin__;
            j4 = (*n0 - 2 << 2) - *pp;
            j4p2 = j4 + (*pp << 1) - 1;
            z__[j4 - 2] = *dnm2 + z__[j4p2];
            z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
            *dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
            *dmin__ = min(*dmin__, *dnm1);
            *dmin1 = *dmin__;
            j4 += 4;
            j4p2 = j4 + (*pp << 1) - 1;
            z__[j4 - 2] = *dnm1 + z__[j4p2];
            z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
            *dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
            *dmin__ = min(*dmin__, *dn);
        } else {
            if (*pp == 0) {
                i__1 = *n0 - 3 << 2;
                for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
                    z__[j4 - 2] = d__ + z__[j4 - 1];
                    if (d__ < 0.) {
                        return 0;
                    } else {
                        z__[j4] = z__[j4 + 1] * (z__[j4 - 1] / z__[j4 - 2]);
                        d__ = z__[j4 + 1] * (d__ / z__[j4 - 2]) - *tau;
                    }
                    if (d__ < dthresh) {
                        d__ = 0.;
                    }
                    *dmin__ = min(*dmin__, d__);
                    d__1 = emin, d__2 = z__[j4];
                    emin = min(d__1, d__2);
                }
            } else {
                i__1 = *n0 - 3 << 2;
                for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
                    z__[j4 - 3] = d__ + z__[j4];
                    if (d__ < 0.) {
                        return 0;
                    } else {
                        z__[j4 - 1] = z__[j4 + 2] * (z__[j4] / z__[j4 - 3]);
                        d__ = z__[j4 + 2] * (d__ / z__[j4 - 3]) - *tau;
                    }
                    if (d__ < dthresh) {
                        d__ = 0.;
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
            if (*dnm2 < 0.) {
                return 0;
            } else {
                z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
                *dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
            }
            *dmin__ = min(*dmin__, *dnm1);
            *dmin1 = *dmin__;
            j4 += 4;
            j4p2 = j4 + (*pp << 1) - 1;
            z__[j4 - 2] = *dnm1 + z__[j4p2];
            if (*dnm1 < 0.) {
                return 0;
            } else {
                z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
                *dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
            }
            *dmin__ = min(*dmin__, *dn);
        }
    }
    z__[j4 + 2] = *dn;
    z__[(*n0 << 2) - *pp] = emin;
    return 0;
}
#ifdef __cplusplus
}
#endif
