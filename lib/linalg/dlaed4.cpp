#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlaed4_(integer *n, integer *i__, doublereal *d__, doublereal *z__, doublereal *delta,
            doublereal *rho, doublereal *dlam, integer *info)
{
    integer i__1;
    doublereal d__1;
    double sqrt(doublereal);
    doublereal a, b, c__;
    integer j;
    doublereal w;
    integer ii;
    doublereal dw, zz[3];
    integer ip1;
    doublereal del, eta, phi, eps, tau, psi;
    integer iim1, iip1;
    doublereal dphi, dpsi;
    integer iter;
    doublereal temp, prew, temp1, dltlb, dltub, midpt;
    integer niter;
    logical swtch;
    extern int dlaed5_(integer *, doublereal *, doublereal *, doublereal *, doublereal *,
                       doublereal *),
        dlaed6_(integer *, logical *, doublereal *, doublereal *, doublereal *, doublereal *,
                doublereal *, integer *);
    logical swtch3;
    extern doublereal dlamch_(char *, ftnlen);
    logical orgati;
    doublereal erretm, rhoinv;
    --delta;
    --z__;
    --d__;
    *info = 0;
    if (*n == 1) {
        *dlam = d__[1] + *rho * z__[1] * z__[1];
        delta[1] = 1.;
        return 0;
    }
    if (*n == 2) {
        dlaed5_(i__, &d__[1], &z__[1], &delta[1], rho, dlam);
        return 0;
    }
    eps = dlamch_((char *)"Epsilon", (ftnlen)7);
    rhoinv = 1. / *rho;
    if (*i__ == *n) {
        ii = *n - 1;
        niter = 1;
        midpt = *rho / 2.;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            delta[j] = d__[j] - d__[*i__] - midpt;
        }
        psi = 0.;
        i__1 = *n - 2;
        for (j = 1; j <= i__1; ++j) {
            psi += z__[j] * z__[j] / delta[j];
        }
        c__ = rhoinv + psi;
        w = c__ + z__[ii] * z__[ii] / delta[ii] + z__[*n] * z__[*n] / delta[*n];
        if (w <= 0.) {
            temp = z__[*n - 1] * z__[*n - 1] / (d__[*n] - d__[*n - 1] + *rho) +
                   z__[*n] * z__[*n] / *rho;
            if (c__ <= temp) {
                tau = *rho;
            } else {
                del = d__[*n] - d__[*n - 1];
                a = -c__ * del + z__[*n - 1] * z__[*n - 1] + z__[*n] * z__[*n];
                b = z__[*n] * z__[*n] * del;
                if (a < 0.) {
                    tau = b * 2. / (sqrt(a * a + b * 4. * c__) - a);
                } else {
                    tau = (a + sqrt(a * a + b * 4. * c__)) / (c__ * 2.);
                }
            }
            dltlb = midpt;
            dltub = *rho;
        } else {
            del = d__[*n] - d__[*n - 1];
            a = -c__ * del + z__[*n - 1] * z__[*n - 1] + z__[*n] * z__[*n];
            b = z__[*n] * z__[*n] * del;
            if (a < 0.) {
                tau = b * 2. / (sqrt(a * a + b * 4. * c__) - a);
            } else {
                tau = (a + sqrt(a * a + b * 4. * c__)) / (c__ * 2.);
            }
            dltlb = 0.;
            dltub = midpt;
        }
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            delta[j] = d__[j] - d__[*i__] - tau;
        }
        dpsi = 0.;
        psi = 0.;
        erretm = 0.;
        i__1 = ii;
        for (j = 1; j <= i__1; ++j) {
            temp = z__[j] / delta[j];
            psi += z__[j] * temp;
            dpsi += temp * temp;
            erretm += psi;
        }
        erretm = abs(erretm);
        temp = z__[*n] / delta[*n];
        phi = z__[*n] * temp;
        dphi = temp * temp;
        erretm = (-phi - psi) * 8. + erretm - phi + rhoinv + abs(tau) * (dpsi + dphi);
        w = rhoinv + phi + psi;
        if (abs(w) <= eps * erretm) {
            *dlam = d__[*i__] + tau;
            goto L250;
        }
        if (w <= 0.) {
            dltlb = max(dltlb, tau);
        } else {
            dltub = min(dltub, tau);
        }
        ++niter;
        c__ = w - delta[*n - 1] * dpsi - delta[*n] * dphi;
        a = (delta[*n - 1] + delta[*n]) * w - delta[*n - 1] * delta[*n] * (dpsi + dphi);
        b = delta[*n - 1] * delta[*n] * w;
        if (c__ < 0.) {
            c__ = abs(c__);
        }
        if (c__ == 0.) {
            eta = -w / (dpsi + dphi);
        } else if (a >= 0.) {
            eta = (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (c__ * 2.);
        } else {
            eta = b * 2. / (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))));
        }
        if (w * eta > 0.) {
            eta = -w / (dpsi + dphi);
        }
        temp = tau + eta;
        if (temp > dltub || temp < dltlb) {
            if (w < 0.) {
                eta = (dltub - tau) / 2.;
            } else {
                eta = (dltlb - tau) / 2.;
            }
        }
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            delta[j] -= eta;
        }
        tau += eta;
        dpsi = 0.;
        psi = 0.;
        erretm = 0.;
        i__1 = ii;
        for (j = 1; j <= i__1; ++j) {
            temp = z__[j] / delta[j];
            psi += z__[j] * temp;
            dpsi += temp * temp;
            erretm += psi;
        }
        erretm = abs(erretm);
        temp = z__[*n] / delta[*n];
        phi = z__[*n] * temp;
        dphi = temp * temp;
        erretm = (-phi - psi) * 8. + erretm - phi + rhoinv + abs(tau) * (dpsi + dphi);
        w = rhoinv + phi + psi;
        iter = niter + 1;
        for (niter = iter; niter <= 30; ++niter) {
            if (abs(w) <= eps * erretm) {
                *dlam = d__[*i__] + tau;
                goto L250;
            }
            if (w <= 0.) {
                dltlb = max(dltlb, tau);
            } else {
                dltub = min(dltub, tau);
            }
            c__ = w - delta[*n - 1] * dpsi - delta[*n] * dphi;
            a = (delta[*n - 1] + delta[*n]) * w - delta[*n - 1] * delta[*n] * (dpsi + dphi);
            b = delta[*n - 1] * delta[*n] * w;
            if (a >= 0.) {
                eta = (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (c__ * 2.);
            } else {
                eta = b * 2. / (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))));
            }
            if (w * eta > 0.) {
                eta = -w / (dpsi + dphi);
            }
            temp = tau + eta;
            if (temp > dltub || temp < dltlb) {
                if (w < 0.) {
                    eta = (dltub - tau) / 2.;
                } else {
                    eta = (dltlb - tau) / 2.;
                }
            }
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                delta[j] -= eta;
            }
            tau += eta;
            dpsi = 0.;
            psi = 0.;
            erretm = 0.;
            i__1 = ii;
            for (j = 1; j <= i__1; ++j) {
                temp = z__[j] / delta[j];
                psi += z__[j] * temp;
                dpsi += temp * temp;
                erretm += psi;
            }
            erretm = abs(erretm);
            temp = z__[*n] / delta[*n];
            phi = z__[*n] * temp;
            dphi = temp * temp;
            erretm = (-phi - psi) * 8. + erretm - phi + rhoinv + abs(tau) * (dpsi + dphi);
            w = rhoinv + phi + psi;
        }
        *info = 1;
        *dlam = d__[*i__] + tau;
        goto L250;
    } else {
        niter = 1;
        ip1 = *i__ + 1;
        del = d__[ip1] - d__[*i__];
        midpt = del / 2.;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            delta[j] = d__[j] - d__[*i__] - midpt;
        }
        psi = 0.;
        i__1 = *i__ - 1;
        for (j = 1; j <= i__1; ++j) {
            psi += z__[j] * z__[j] / delta[j];
        }
        phi = 0.;
        i__1 = *i__ + 2;
        for (j = *n; j >= i__1; --j) {
            phi += z__[j] * z__[j] / delta[j];
        }
        c__ = rhoinv + psi + phi;
        w = c__ + z__[*i__] * z__[*i__] / delta[*i__] + z__[ip1] * z__[ip1] / delta[ip1];
        if (w > 0.) {
            orgati = TRUE_;
            a = c__ * del + z__[*i__] * z__[*i__] + z__[ip1] * z__[ip1];
            b = z__[*i__] * z__[*i__] * del;
            if (a > 0.) {
                tau = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))));
            } else {
                tau = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (c__ * 2.);
            }
            dltlb = 0.;
            dltub = midpt;
        } else {
            orgati = FALSE_;
            a = c__ * del - z__[*i__] * z__[*i__] - z__[ip1] * z__[ip1];
            b = z__[ip1] * z__[ip1] * del;
            if (a < 0.) {
                tau = b * 2. / (a - sqrt((d__1 = a * a + b * 4. * c__, abs(d__1))));
            } else {
                tau = -(a + sqrt((d__1 = a * a + b * 4. * c__, abs(d__1)))) / (c__ * 2.);
            }
            dltlb = -midpt;
            dltub = 0.;
        }
        if (orgati) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                delta[j] = d__[j] - d__[*i__] - tau;
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                delta[j] = d__[j] - d__[ip1] - tau;
            }
        }
        if (orgati) {
            ii = *i__;
        } else {
            ii = *i__ + 1;
        }
        iim1 = ii - 1;
        iip1 = ii + 1;
        dpsi = 0.;
        psi = 0.;
        erretm = 0.;
        i__1 = iim1;
        for (j = 1; j <= i__1; ++j) {
            temp = z__[j] / delta[j];
            psi += z__[j] * temp;
            dpsi += temp * temp;
            erretm += psi;
        }
        erretm = abs(erretm);
        dphi = 0.;
        phi = 0.;
        i__1 = iip1;
        for (j = *n; j >= i__1; --j) {
            temp = z__[j] / delta[j];
            phi += z__[j] * temp;
            dphi += temp * temp;
            erretm += phi;
        }
        w = rhoinv + phi + psi;
        swtch3 = FALSE_;
        if (orgati) {
            if (w < 0.) {
                swtch3 = TRUE_;
            }
        } else {
            if (w > 0.) {
                swtch3 = TRUE_;
            }
        }
        if (ii == 1 || ii == *n) {
            swtch3 = FALSE_;
        }
        temp = z__[ii] / delta[ii];
        dw = dpsi + dphi + temp * temp;
        temp = z__[ii] * temp;
        w += temp;
        erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + abs(temp) * 3. + abs(tau) * dw;
        if (abs(w) <= eps * erretm) {
            if (orgati) {
                *dlam = d__[*i__] + tau;
            } else {
                *dlam = d__[ip1] + tau;
            }
            goto L250;
        }
        if (w <= 0.) {
            dltlb = max(dltlb, tau);
        } else {
            dltub = min(dltub, tau);
        }
        ++niter;
        if (!swtch3) {
            if (orgati) {
                d__1 = z__[*i__] / delta[*i__];
                c__ = w - delta[ip1] * dw - (d__[*i__] - d__[ip1]) * (d__1 * d__1);
            } else {
                d__1 = z__[ip1] / delta[ip1];
                c__ = w - delta[*i__] * dw - (d__[ip1] - d__[*i__]) * (d__1 * d__1);
            }
            a = (delta[*i__] + delta[ip1]) * w - delta[*i__] * delta[ip1] * dw;
            b = delta[*i__] * delta[ip1] * w;
            if (c__ == 0.) {
                if (a == 0.) {
                    if (orgati) {
                        a = z__[*i__] * z__[*i__] + delta[ip1] * delta[ip1] * (dpsi + dphi);
                    } else {
                        a = z__[ip1] * z__[ip1] + delta[*i__] * delta[*i__] * (dpsi + dphi);
                    }
                }
                eta = b / a;
            } else if (a <= 0.) {
                eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (c__ * 2.);
            } else {
                eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))));
            }
        } else {
            temp = rhoinv + psi + phi;
            if (orgati) {
                temp1 = z__[iim1] / delta[iim1];
                temp1 *= temp1;
                c__ = temp - delta[iip1] * (dpsi + dphi) - (d__[iim1] - d__[iip1]) * temp1;
                zz[0] = z__[iim1] * z__[iim1];
                zz[2] = delta[iip1] * delta[iip1] * (dpsi - temp1 + dphi);
            } else {
                temp1 = z__[iip1] / delta[iip1];
                temp1 *= temp1;
                c__ = temp - delta[iim1] * (dpsi + dphi) - (d__[iip1] - d__[iim1]) * temp1;
                zz[0] = delta[iim1] * delta[iim1] * (dpsi + (dphi - temp1));
                zz[2] = z__[iip1] * z__[iip1];
            }
            zz[1] = z__[ii] * z__[ii];
            dlaed6_(&niter, &orgati, &c__, &delta[iim1], zz, &w, &eta, info);
            if (*info != 0) {
                goto L250;
            }
        }
        if (w * eta >= 0.) {
            eta = -w / dw;
        }
        temp = tau + eta;
        if (temp > dltub || temp < dltlb) {
            if (w < 0.) {
                eta = (dltub - tau) / 2.;
            } else {
                eta = (dltlb - tau) / 2.;
            }
        }
        prew = w;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            delta[j] -= eta;
        }
        dpsi = 0.;
        psi = 0.;
        erretm = 0.;
        i__1 = iim1;
        for (j = 1; j <= i__1; ++j) {
            temp = z__[j] / delta[j];
            psi += z__[j] * temp;
            dpsi += temp * temp;
            erretm += psi;
        }
        erretm = abs(erretm);
        dphi = 0.;
        phi = 0.;
        i__1 = iip1;
        for (j = *n; j >= i__1; --j) {
            temp = z__[j] / delta[j];
            phi += z__[j] * temp;
            dphi += temp * temp;
            erretm += phi;
        }
        temp = z__[ii] / delta[ii];
        dw = dpsi + dphi + temp * temp;
        temp = z__[ii] * temp;
        w = rhoinv + phi + psi + temp;
        erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + abs(temp) * 3. +
                 (d__1 = tau + eta, abs(d__1)) * dw;
        swtch = FALSE_;
        if (orgati) {
            if (-w > abs(prew) / 10.) {
                swtch = TRUE_;
            }
        } else {
            if (w > abs(prew) / 10.) {
                swtch = TRUE_;
            }
        }
        tau += eta;
        iter = niter + 1;
        for (niter = iter; niter <= 30; ++niter) {
            if (abs(w) <= eps * erretm) {
                if (orgati) {
                    *dlam = d__[*i__] + tau;
                } else {
                    *dlam = d__[ip1] + tau;
                }
                goto L250;
            }
            if (w <= 0.) {
                dltlb = max(dltlb, tau);
            } else {
                dltub = min(dltub, tau);
            }
            if (!swtch3) {
                if (!swtch) {
                    if (orgati) {
                        d__1 = z__[*i__] / delta[*i__];
                        c__ = w - delta[ip1] * dw - (d__[*i__] - d__[ip1]) * (d__1 * d__1);
                    } else {
                        d__1 = z__[ip1] / delta[ip1];
                        c__ = w - delta[*i__] * dw - (d__[ip1] - d__[*i__]) * (d__1 * d__1);
                    }
                } else {
                    temp = z__[ii] / delta[ii];
                    if (orgati) {
                        dpsi += temp * temp;
                    } else {
                        dphi += temp * temp;
                    }
                    c__ = w - delta[*i__] * dpsi - delta[ip1] * dphi;
                }
                a = (delta[*i__] + delta[ip1]) * w - delta[*i__] * delta[ip1] * dw;
                b = delta[*i__] * delta[ip1] * w;
                if (c__ == 0.) {
                    if (a == 0.) {
                        if (!swtch) {
                            if (orgati) {
                                a = z__[*i__] * z__[*i__] + delta[ip1] * delta[ip1] * (dpsi + dphi);
                            } else {
                                a = z__[ip1] * z__[ip1] + delta[*i__] * delta[*i__] * (dpsi + dphi);
                            }
                        } else {
                            a = delta[*i__] * delta[*i__] * dpsi + delta[ip1] * delta[ip1] * dphi;
                        }
                    }
                    eta = b / a;
                } else if (a <= 0.) {
                    eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (c__ * 2.);
                } else {
                    eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))));
                }
            } else {
                temp = rhoinv + psi + phi;
                if (swtch) {
                    c__ = temp - delta[iim1] * dpsi - delta[iip1] * dphi;
                    zz[0] = delta[iim1] * delta[iim1] * dpsi;
                    zz[2] = delta[iip1] * delta[iip1] * dphi;
                } else {
                    if (orgati) {
                        temp1 = z__[iim1] / delta[iim1];
                        temp1 *= temp1;
                        c__ = temp - delta[iip1] * (dpsi + dphi) - (d__[iim1] - d__[iip1]) * temp1;
                        zz[0] = z__[iim1] * z__[iim1];
                        zz[2] = delta[iip1] * delta[iip1] * (dpsi - temp1 + dphi);
                    } else {
                        temp1 = z__[iip1] / delta[iip1];
                        temp1 *= temp1;
                        c__ = temp - delta[iim1] * (dpsi + dphi) - (d__[iip1] - d__[iim1]) * temp1;
                        zz[0] = delta[iim1] * delta[iim1] * (dpsi + (dphi - temp1));
                        zz[2] = z__[iip1] * z__[iip1];
                    }
                }
                dlaed6_(&niter, &orgati, &c__, &delta[iim1], zz, &w, &eta, info);
                if (*info != 0) {
                    goto L250;
                }
            }
            if (w * eta >= 0.) {
                eta = -w / dw;
            }
            temp = tau + eta;
            if (temp > dltub || temp < dltlb) {
                if (w < 0.) {
                    eta = (dltub - tau) / 2.;
                } else {
                    eta = (dltlb - tau) / 2.;
                }
            }
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                delta[j] -= eta;
            }
            tau += eta;
            prew = w;
            dpsi = 0.;
            psi = 0.;
            erretm = 0.;
            i__1 = iim1;
            for (j = 1; j <= i__1; ++j) {
                temp = z__[j] / delta[j];
                psi += z__[j] * temp;
                dpsi += temp * temp;
                erretm += psi;
            }
            erretm = abs(erretm);
            dphi = 0.;
            phi = 0.;
            i__1 = iip1;
            for (j = *n; j >= i__1; --j) {
                temp = z__[j] / delta[j];
                phi += z__[j] * temp;
                dphi += temp * temp;
                erretm += phi;
            }
            temp = z__[ii] / delta[ii];
            dw = dpsi + dphi + temp * temp;
            temp = z__[ii] * temp;
            w = rhoinv + phi + psi + temp;
            erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + abs(temp) * 3. + abs(tau) * dw;
            if (w * prew > 0. && abs(w) > abs(prew) / 10.) {
                swtch = !swtch;
            }
        }
        *info = 1;
        if (orgati) {
            *dlam = d__[*i__] + tau;
        } else {
            *dlam = d__[ip1] + tau;
        }
    }
L250:
    return 0;
}
#ifdef __cplusplus
}
#endif
