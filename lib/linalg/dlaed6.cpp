#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlaed6_(integer *kniter, logical *orgati, doublereal *rho, doublereal *d__, doublereal *z__,
            doublereal *finit, doublereal *tau, integer *info)
{
    integer i__1;
    doublereal d__1, d__2, d__3, d__4;
    double sqrt(doublereal), log(doublereal), pow_lmp_di(doublereal *, integer *);
    doublereal a, b, c__, f;
    integer i__;
    doublereal fc, df, ddf, lbd, eta, ubd, eps, base;
    integer iter;
    doublereal temp, temp1, temp2, temp3, temp4;
    logical scale;
    integer niter;
    doublereal small1, small2, sminv1, sminv2;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal dscale[3], sclfac, zscale[3], erretm, sclinv;
    --z__;
    --d__;
    *info = 0;
    if (*orgati) {
        lbd = d__[2];
        ubd = d__[3];
    } else {
        lbd = d__[1];
        ubd = d__[2];
    }
    if (*finit < 0.) {
        lbd = 0.;
    } else {
        ubd = 0.;
    }
    niter = 1;
    *tau = 0.;
    if (*kniter == 2) {
        if (*orgati) {
            temp = (d__[3] - d__[2]) / 2.;
            c__ = *rho + z__[1] / (d__[1] - d__[2] - temp);
            a = c__ * (d__[2] + d__[3]) + z__[2] + z__[3];
            b = c__ * d__[2] * d__[3] + z__[2] * d__[3] + z__[3] * d__[2];
        } else {
            temp = (d__[1] - d__[2]) / 2.;
            c__ = *rho + z__[3] / (d__[3] - d__[2] - temp);
            a = c__ * (d__[1] + d__[2]) + z__[1] + z__[2];
            b = c__ * d__[1] * d__[2] + z__[1] * d__[2] + z__[2] * d__[1];
        }
        d__1 = abs(a), d__2 = abs(b), d__1 = max(d__1, d__2), d__2 = abs(c__);
        temp = max(d__1, d__2);
        a /= temp;
        b /= temp;
        c__ /= temp;
        if (c__ == 0.) {
            *tau = b / a;
        } else if (a <= 0.) {
            *tau = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (c__ * 2.);
        } else {
            *tau = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))));
        }
        if (*tau < lbd || *tau > ubd) {
            *tau = (lbd + ubd) / 2.;
        }
        if (d__[1] == *tau || d__[2] == *tau || d__[3] == *tau) {
            *tau = 0.;
        } else {
            temp = *finit + *tau * z__[1] / (d__[1] * (d__[1] - *tau)) +
                   *tau * z__[2] / (d__[2] * (d__[2] - *tau)) +
                   *tau * z__[3] / (d__[3] * (d__[3] - *tau));
            if (temp <= 0.) {
                lbd = *tau;
            } else {
                ubd = *tau;
            }
            if (abs(*finit) <= abs(temp)) {
                *tau = 0.;
            }
        }
    }
    eps = dlamch_((char *)"Epsilon", (ftnlen)7);
    base = dlamch_((char *)"Base", (ftnlen)4);
    i__1 = (integer)(log(dlamch_((char *)"SafMin", (ftnlen)6)) / log(base) / 3.);
    small1 = pow_lmp_di(&base, &i__1);
    sminv1 = 1. / small1;
    small2 = small1 * small1;
    sminv2 = sminv1 * sminv1;
    if (*orgati) {
        d__3 = (d__1 = d__[2] - *tau, abs(d__1)), d__4 = (d__2 = d__[3] - *tau, abs(d__2));
        temp = min(d__3, d__4);
    } else {
        d__3 = (d__1 = d__[1] - *tau, abs(d__1)), d__4 = (d__2 = d__[2] - *tau, abs(d__2));
        temp = min(d__3, d__4);
    }
    scale = FALSE_;
    if (temp <= small1) {
        scale = TRUE_;
        if (temp <= small2) {
            sclfac = sminv2;
            sclinv = small2;
        } else {
            sclfac = sminv1;
            sclinv = small1;
        }
        for (i__ = 1; i__ <= 3; ++i__) {
            dscale[i__ - 1] = d__[i__] * sclfac;
            zscale[i__ - 1] = z__[i__] * sclfac;
        }
        *tau *= sclfac;
        lbd *= sclfac;
        ubd *= sclfac;
    } else {
        for (i__ = 1; i__ <= 3; ++i__) {
            dscale[i__ - 1] = d__[i__];
            zscale[i__ - 1] = z__[i__];
        }
    }
    fc = 0.;
    df = 0.;
    ddf = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
        temp = 1. / (dscale[i__ - 1] - *tau);
        temp1 = zscale[i__ - 1] * temp;
        temp2 = temp1 * temp;
        temp3 = temp2 * temp;
        fc += temp1 / dscale[i__ - 1];
        df += temp2;
        ddf += temp3;
    }
    f = *finit + *tau * fc;
    if (abs(f) <= 0.) {
        goto L60;
    }
    if (f <= 0.) {
        lbd = *tau;
    } else {
        ubd = *tau;
    }
    iter = niter + 1;
    for (niter = iter; niter <= 40; ++niter) {
        if (*orgati) {
            temp1 = dscale[1] - *tau;
            temp2 = dscale[2] - *tau;
        } else {
            temp1 = dscale[0] - *tau;
            temp2 = dscale[1] - *tau;
        }
        a = (temp1 + temp2) * f - temp1 * temp2 * df;
        b = temp1 * temp2 * f;
        c__ = f - (temp1 + temp2) * df + temp1 * temp2 * ddf;
        d__1 = abs(a), d__2 = abs(b), d__1 = max(d__1, d__2), d__2 = abs(c__);
        temp = max(d__1, d__2);
        a /= temp;
        b /= temp;
        c__ /= temp;
        if (c__ == 0.) {
            eta = b / a;
        } else if (a <= 0.) {
            eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (c__ * 2.);
        } else {
            eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))));
        }
        if (f * eta >= 0.) {
            eta = -f / df;
        }
        *tau += eta;
        if (*tau < lbd || *tau > ubd) {
            *tau = (lbd + ubd) / 2.;
        }
        fc = 0.;
        erretm = 0.;
        df = 0.;
        ddf = 0.;
        for (i__ = 1; i__ <= 3; ++i__) {
            if (dscale[i__ - 1] - *tau != 0.) {
                temp = 1. / (dscale[i__ - 1] - *tau);
                temp1 = zscale[i__ - 1] * temp;
                temp2 = temp1 * temp;
                temp3 = temp2 * temp;
                temp4 = temp1 / dscale[i__ - 1];
                fc += temp4;
                erretm += abs(temp4);
                df += temp2;
                ddf += temp3;
            } else {
                goto L60;
            }
        }
        f = *finit + *tau * fc;
        erretm = (abs(*finit) + abs(*tau) * erretm) * 8. + abs(*tau) * df;
        if (abs(f) <= eps * 4. * erretm || ubd - lbd <= eps * 4. * abs(*tau)) {
            goto L60;
        }
        if (f <= 0.) {
            lbd = *tau;
        } else {
            ubd = *tau;
        }
    }
    *info = 1;
L60:
    if (scale) {
        *tau *= sclinv;
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
