#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlaed5_(integer *i__, doublereal *d__, doublereal *z__, doublereal *delta, doublereal *rho,
            doublereal *dlam)
{
    doublereal d__1;
    double sqrt(doublereal);
    doublereal b, c__, w, del, tau, temp;
    --delta;
    --z__;
    --d__;
    del = d__[2] - d__[1];
    if (*i__ == 1) {
        w = *rho * 2. * (z__[2] * z__[2] - z__[1] * z__[1]) / del + 1.;
        if (w > 0.) {
            b = del + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
            c__ = *rho * z__[1] * z__[1] * del;
            tau = c__ * 2. / (b + sqrt((d__1 = b * b - c__ * 4., abs(d__1))));
            *dlam = d__[1] + tau;
            delta[1] = -z__[1] / tau;
            delta[2] = z__[2] / (del - tau);
        } else {
            b = -del + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
            c__ = *rho * z__[2] * z__[2] * del;
            if (b > 0.) {
                tau = c__ * -2. / (b + sqrt(b * b + c__ * 4.));
            } else {
                tau = (b - sqrt(b * b + c__ * 4.)) / 2.;
            }
            *dlam = d__[2] + tau;
            delta[1] = -z__[1] / (del + tau);
            delta[2] = -z__[2] / tau;
        }
        temp = sqrt(delta[1] * delta[1] + delta[2] * delta[2]);
        delta[1] /= temp;
        delta[2] /= temp;
    } else {
        b = -del + *rho * (z__[1] * z__[1] + z__[2] * z__[2]);
        c__ = *rho * z__[2] * z__[2] * del;
        if (b > 0.) {
            tau = (b + sqrt(b * b + c__ * 4.)) / 2.;
        } else {
            tau = c__ * 2. / (-b + sqrt(b * b + c__ * 4.));
        }
        *dlam = d__[2] + tau;
        delta[1] = -z__[1] / (del + tau);
        delta[2] = -z__[2] / tau;
        temp = sqrt(delta[1] * delta[1] + delta[2] * delta[2]);
        delta[1] /= temp;
        delta[2] /= temp;
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
