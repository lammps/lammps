#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b6 = 1.;
int dlanv2_(doublereal *a, doublereal *b, doublereal *c__, doublereal *d__, doublereal *rt1r,
            doublereal *rt1i, doublereal *rt2r, doublereal *rt2i, doublereal *cs, doublereal *sn)
{
    integer i__1;
    doublereal d__1, d__2;
    double log(doublereal), pow_lmp_di(doublereal *, integer *), d_lmp_sign(doublereal *, doublereal *),
        sqrt(doublereal);
    doublereal p, z__, aa, bb, cc, dd, cs1, sn1, sab, sac, eps, tau, temp, scale, bcmax, bcmis,
        sigma;
    integer count;
    doublereal safmn2;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    doublereal safmx2;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal safmin;
    safmin = dlamch_((char *)"S", (ftnlen)1);
    eps = dlamch_((char *)"P", (ftnlen)1);
    d__1 = dlamch_((char *)"B", (ftnlen)1);
    i__1 = (integer)(log(safmin / eps) / log(dlamch_((char *)"B", (ftnlen)1)) / 2.);
    safmn2 = pow_lmp_di(&d__1, &i__1);
    safmx2 = 1. / safmn2;
    if (*c__ == 0.) {
        *cs = 1.;
        *sn = 0.;
    } else if (*b == 0.) {
        *cs = 0.;
        *sn = 1.;
        temp = *d__;
        *d__ = *a;
        *a = temp;
        *b = -(*c__);
        *c__ = 0.;
    } else if (*a - *d__ == 0. && d_lmp_sign(&c_b6, b) != d_lmp_sign(&c_b6, c__)) {
        *cs = 1.;
        *sn = 0.;
    } else {
        temp = *a - *d__;
        p = temp * .5;
        d__1 = abs(*b), d__2 = abs(*c__);
        bcmax = max(d__1, d__2);
        d__1 = abs(*b), d__2 = abs(*c__);
        bcmis = min(d__1, d__2) * d_lmp_sign(&c_b6, b) * d_lmp_sign(&c_b6, c__);
        d__1 = abs(p);
        scale = max(d__1, bcmax);
        z__ = p / scale * p + bcmax / scale * bcmis;
        if (z__ >= eps * 4.) {
            d__1 = sqrt(scale) * sqrt(z__);
            z__ = p + d_lmp_sign(&d__1, &p);
            *a = *d__ + z__;
            *d__ -= bcmax / z__ * bcmis;
            tau = dlapy2_(c__, &z__);
            *cs = z__ / tau;
            *sn = *c__ / tau;
            *b -= *c__;
            *c__ = 0.;
        } else {
            count = 0;
            sigma = *b + *c__;
        L10:
            ++count;
            d__1 = abs(temp), d__2 = abs(sigma);
            scale = max(d__1, d__2);
            if (scale >= safmx2) {
                sigma *= safmn2;
                temp *= safmn2;
                if (count <= 20) {
                    goto L10;
                }
            }
            if (scale <= safmn2) {
                sigma *= safmx2;
                temp *= safmx2;
                if (count <= 20) {
                    goto L10;
                }
            }
            p = temp * .5;
            tau = dlapy2_(&sigma, &temp);
            *cs = sqrt((abs(sigma) / tau + 1.) * .5);
            *sn = -(p / (tau * *cs)) * d_lmp_sign(&c_b6, &sigma);
            aa = *a * *cs + *b * *sn;
            bb = -(*a) * *sn + *b * *cs;
            cc = *c__ * *cs + *d__ * *sn;
            dd = -(*c__) * *sn + *d__ * *cs;
            *a = aa * *cs + cc * *sn;
            *b = bb * *cs + dd * *sn;
            *c__ = -(aa * *sn) + cc * *cs;
            *d__ = -bb * *sn + dd * *cs;
            temp = (*a + *d__) * .5;
            *a = temp;
            *d__ = temp;
            if (*c__ != 0.) {
                if (*b != 0.) {
                    if (d_lmp_sign(&c_b6, b) == d_lmp_sign(&c_b6, c__)) {
                        sab = sqrt((abs(*b)));
                        sac = sqrt((abs(*c__)));
                        d__1 = sab * sac;
                        p = d_lmp_sign(&d__1, c__);
                        tau = 1. / sqrt((d__1 = *b + *c__, abs(d__1)));
                        *a = temp + p;
                        *d__ = temp - p;
                        *b -= *c__;
                        *c__ = 0.;
                        cs1 = sab * tau;
                        sn1 = sac * tau;
                        temp = *cs * cs1 - *sn * sn1;
                        *sn = *cs * sn1 + *sn * cs1;
                        *cs = temp;
                    }
                } else {
                    *b = -(*c__);
                    *c__ = 0.;
                    temp = *cs;
                    *cs = -(*sn);
                    *sn = temp;
                }
            }
        }
    }
    *rt1r = *a;
    *rt2r = *d__;
    if (*c__ == 0.) {
        *rt1i = 0.;
        *rt2i = 0.;
    } else {
        *rt1i = sqrt((abs(*b))) * sqrt((abs(*c__)));
        *rt2i = -(*rt1i);
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
