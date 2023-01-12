#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlasd5_(integer *i__, doublereal *d__, doublereal *z__, doublereal *delta, doublereal *rho,
            doublereal *dsigma, doublereal *work)
{
    doublereal d__1;
    double sqrt(doublereal);
    doublereal b, c__, w, del, tau, delsq;
    --work;
    --delta;
    --z__;
    --d__;
    del = d__[2] - d__[1];
    delsq = del * (d__[2] + d__[1]);
    if (*i__ == 1) {
        w = *rho * 4. *
                (z__[2] * z__[2] / (d__[1] + d__[2] * 3.) -
                 z__[1] * z__[1] / (d__[1] * 3. + d__[2])) /
                del +
            1.;
        if (w > 0.) {
            b = delsq + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
            c__ = *rho * z__[1] * z__[1] * delsq;
            tau = c__ * 2. / (b + sqrt((d__1 = b * b - c__ * 4., abs(d__1))));
            tau /= d__[1] + sqrt(d__[1] * d__[1] + tau);
            *dsigma = d__[1] + tau;
            delta[1] = -tau;
            delta[2] = del - tau;
            work[1] = d__[1] * 2. + tau;
            work[2] = d__[1] + tau + d__[2];
        } else {
            b = -delsq + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
            c__ = *rho * z__[2] * z__[2] * delsq;
            if (b > 0.) {
                tau = c__ * -2. / (b + sqrt(b * b + c__ * 4.));
            } else {
                tau = (b - sqrt(b * b + c__ * 4.)) / 2.;
            }
            tau /= d__[2] + sqrt((d__1 = d__[2] * d__[2] + tau, abs(d__1)));
            *dsigma = d__[2] + tau;
            delta[1] = -(del + tau);
            delta[2] = -tau;
            work[1] = d__[1] + tau + d__[2];
            work[2] = d__[2] * 2. + tau;
        }
    } else {
        b = -delsq + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
        c__ = *rho * z__[2] * z__[2] * delsq;
        if (b > 0.) {
            tau = (b + sqrt(b * b + c__ * 4.)) / 2.;
        } else {
            tau = c__ * 2. / (-b + sqrt(b * b + c__ * 4.));
        }
        tau /= d__[2] + sqrt(d__[2] * d__[2] + tau);
        *dsigma = d__[2] + tau;
        delta[1] = -(del + tau);
        delta[2] = -tau;
        work[1] = d__[1] + tau + d__[2];
        work[2] = d__[2] * 2. + tau;
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
