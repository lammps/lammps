#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlasd4_(integer *n, integer *i__, doublereal *d__, doublereal *z__, doublereal *delta,
            doublereal *rho, doublereal *sigma, doublereal *work, integer *info)
{
    integer i__1;
    doublereal d__1;
    double sqrt(doublereal);
    doublereal a, b, c__;
    integer j;
    doublereal w, dd[3];
    integer ii;
    doublereal dw, zz[3];
    integer ip1;
    doublereal sq2, eta, phi, eps, tau, psi;
    integer iim1, iip1;
    doublereal tau2, dphi, sglb, dpsi, sgub;
    integer iter;
    doublereal temp, prew, temp1, temp2, dtiim, delsq, dtiip;
    integer niter;
    doublereal dtisq;
    logical swtch;
    doublereal dtnsq;
    extern int dlaed6_(integer *, logical *, doublereal *, doublereal *, doublereal *, doublereal *,
                       doublereal *, integer *),
        dlasd5_(integer *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                doublereal *);
    doublereal delsq2, dtnsq1;
    logical swtch3;
    extern doublereal dlamch_(char *, ftnlen);
    logical orgati;
    doublereal erretm, dtipsq, rhoinv;
    logical geomavg;
    --work;
    --delta;
    --z__;
    --d__;
    *info = 0;
    if (*n == 1) {
        *sigma = sqrt(d__[1] * d__[1] + *rho * z__[1] * z__[1]);
        delta[1] = 1.;
        work[1] = 1.;
        return 0;
    }
    if (*n == 2) {
        dlasd5_(i__, &d__[1], &z__[1], &delta[1], rho, sigma, &work[1]);
        return 0;
    }
    eps = dlamch_((char *)"Epsilon", (ftnlen)7);
    rhoinv = 1. / *rho;
    tau2 = 0.;
    if (*i__ == *n) {
        ii = *n - 1;
        niter = 1;
        temp = *rho / 2.;
        temp1 = temp / (d__[*n] + sqrt(d__[*n] * d__[*n] + temp));
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            work[j] = d__[j] + d__[*n] + temp1;
            delta[j] = d__[j] - d__[*n] - temp1;
        }
        psi = 0.;
        i__1 = *n - 2;
        for (j = 1; j <= i__1; ++j) {
            psi += z__[j] * z__[j] / (delta[j] * work[j]);
        }
        c__ = rhoinv + psi;
        w = c__ + z__[ii] * z__[ii] / (delta[ii] * work[ii]) +
            z__[*n] * z__[*n] / (delta[*n] * work[*n]);
        if (w <= 0.) {
            temp1 = sqrt(d__[*n] * d__[*n] + *rho);
            temp =
                z__[*n - 1] * z__[*n - 1] /
                    ((d__[*n - 1] + temp1) * (d__[*n] - d__[*n - 1] + *rho / (d__[*n] + temp1))) +
                z__[*n] * z__[*n] / *rho;
            if (c__ <= temp) {
                tau = *rho;
            } else {
                delsq = (d__[*n] - d__[*n - 1]) * (d__[*n] + d__[*n - 1]);
                a = -c__ * delsq + z__[*n - 1] * z__[*n - 1] + z__[*n] * z__[*n];
                b = z__[*n] * z__[*n] * delsq;
                if (a < 0.) {
                    tau2 = b * 2. / (sqrt(a * a + b * 4. * c__) - a);
                } else {
                    tau2 = (a + sqrt(a * a + b * 4. * c__)) / (c__ * 2.);
                }
                tau = tau2 / (d__[*n] + sqrt(d__[*n] * d__[*n] + tau2));
            }
        } else {
            delsq = (d__[*n] - d__[*n - 1]) * (d__[*n] + d__[*n - 1]);
            a = -c__ * delsq + z__[*n - 1] * z__[*n - 1] + z__[*n] * z__[*n];
            b = z__[*n] * z__[*n] * delsq;
            if (a < 0.) {
                tau2 = b * 2. / (sqrt(a * a + b * 4. * c__) - a);
            } else {
                tau2 = (a + sqrt(a * a + b * 4. * c__)) / (c__ * 2.);
            }
            tau = tau2 / (d__[*n] + sqrt(d__[*n] * d__[*n] + tau2));
        }
        *sigma = d__[*n] + tau;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            delta[j] = d__[j] - d__[*n] - tau;
            work[j] = d__[j] + d__[*n] + tau;
        }
        dpsi = 0.;
        psi = 0.;
        erretm = 0.;
        i__1 = ii;
        for (j = 1; j <= i__1; ++j) {
            temp = z__[j] / (delta[j] * work[j]);
            psi += z__[j] * temp;
            dpsi += temp * temp;
            erretm += psi;
        }
        erretm = abs(erretm);
        temp = z__[*n] / (delta[*n] * work[*n]);
        phi = z__[*n] * temp;
        dphi = temp * temp;
        erretm = (-phi - psi) * 8. + erretm - phi + rhoinv;
        w = rhoinv + phi + psi;
        if (abs(w) <= eps * erretm) {
            goto L240;
        }
        ++niter;
        dtnsq1 = work[*n - 1] * delta[*n - 1];
        dtnsq = work[*n] * delta[*n];
        c__ = w - dtnsq1 * dpsi - dtnsq * dphi;
        a = (dtnsq + dtnsq1) * w - dtnsq * dtnsq1 * (dpsi + dphi);
        b = dtnsq * dtnsq1 * w;
        if (c__ < 0.) {
            c__ = abs(c__);
        }
        if (c__ == 0.) {
            eta = *rho - *sigma * *sigma;
        } else if (a >= 0.) {
            eta = (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (c__ * 2.);
        } else {
            eta = b * 2. / (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))));
        }
        if (w * eta > 0.) {
            eta = -w / (dpsi + dphi);
        }
        temp = eta - dtnsq;
        if (temp > *rho) {
            eta = *rho + dtnsq;
        }
        eta /= *sigma + sqrt(eta + *sigma * *sigma);
        tau += eta;
        *sigma += eta;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            delta[j] -= eta;
            work[j] += eta;
        }
        dpsi = 0.;
        psi = 0.;
        erretm = 0.;
        i__1 = ii;
        for (j = 1; j <= i__1; ++j) {
            temp = z__[j] / (work[j] * delta[j]);
            psi += z__[j] * temp;
            dpsi += temp * temp;
            erretm += psi;
        }
        erretm = abs(erretm);
        tau2 = work[*n] * delta[*n];
        temp = z__[*n] / tau2;
        phi = z__[*n] * temp;
        dphi = temp * temp;
        erretm = (-phi - psi) * 8. + erretm - phi + rhoinv;
        w = rhoinv + phi + psi;
        iter = niter + 1;
        for (niter = iter; niter <= 400; ++niter) {
            if (abs(w) <= eps * erretm) {
                goto L240;
            }
            dtnsq1 = work[*n - 1] * delta[*n - 1];
            dtnsq = work[*n] * delta[*n];
            c__ = w - dtnsq1 * dpsi - dtnsq * dphi;
            a = (dtnsq + dtnsq1) * w - dtnsq1 * dtnsq * (dpsi + dphi);
            b = dtnsq1 * dtnsq * w;
            if (a >= 0.) {
                eta = (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (c__ * 2.);
            } else {
                eta = b * 2. / (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))));
            }
            if (w * eta > 0.) {
                eta = -w / (dpsi + dphi);
            }
            temp = eta - dtnsq;
            if (temp <= 0.) {
                eta /= 2.;
            }
            eta /= *sigma + sqrt(eta + *sigma * *sigma);
            tau += eta;
            *sigma += eta;
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                delta[j] -= eta;
                work[j] += eta;
            }
            dpsi = 0.;
            psi = 0.;
            erretm = 0.;
            i__1 = ii;
            for (j = 1; j <= i__1; ++j) {
                temp = z__[j] / (work[j] * delta[j]);
                psi += z__[j] * temp;
                dpsi += temp * temp;
                erretm += psi;
            }
            erretm = abs(erretm);
            tau2 = work[*n] * delta[*n];
            temp = z__[*n] / tau2;
            phi = z__[*n] * temp;
            dphi = temp * temp;
            erretm = (-phi - psi) * 8. + erretm - phi + rhoinv;
            w = rhoinv + phi + psi;
        }
        *info = 1;
        goto L240;
    } else {
        niter = 1;
        ip1 = *i__ + 1;
        delsq = (d__[ip1] - d__[*i__]) * (d__[ip1] + d__[*i__]);
        delsq2 = delsq / 2.;
        sq2 = sqrt((d__[*i__] * d__[*i__] + d__[ip1] * d__[ip1]) / 2.);
        temp = delsq2 / (d__[*i__] + sq2);
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            work[j] = d__[j] + d__[*i__] + temp;
            delta[j] = d__[j] - d__[*i__] - temp;
        }
        psi = 0.;
        i__1 = *i__ - 1;
        for (j = 1; j <= i__1; ++j) {
            psi += z__[j] * z__[j] / (work[j] * delta[j]);
        }
        phi = 0.;
        i__1 = *i__ + 2;
        for (j = *n; j >= i__1; --j) {
            phi += z__[j] * z__[j] / (work[j] * delta[j]);
        }
        c__ = rhoinv + psi + phi;
        w = c__ + z__[*i__] * z__[*i__] / (work[*i__] * delta[*i__]) +
            z__[ip1] * z__[ip1] / (work[ip1] * delta[ip1]);
        geomavg = FALSE_;
        if (w > 0.) {
            orgati = TRUE_;
            ii = *i__;
            sglb = 0.;
            sgub = delsq2 / (d__[*i__] + sq2);
            a = c__ * delsq + z__[*i__] * z__[*i__] + z__[ip1] * z__[ip1];
            b = z__[*i__] * z__[*i__] * delsq;
            if (a > 0.) {
                tau2 = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))));
            } else {
                tau2 = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (c__ * 2.);
            }
            tau = tau2 / (d__[*i__] + sqrt(d__[*i__] * d__[*i__] + tau2));
            temp = sqrt(eps);
            if (d__[*i__] <= temp * d__[ip1] && (d__1 = z__[*i__], abs(d__1)) <= temp &&
                d__[*i__] > 0.) {
                d__1 = d__[*i__] * 10.;
                tau = min(d__1, sgub);
                geomavg = TRUE_;
            }
        } else {
            orgati = FALSE_;
            ii = ip1;
            sglb = -delsq2 / (d__[ii] + sq2);
            sgub = 0.;
            a = c__ * delsq - z__[*i__] * z__[*i__] - z__[ip1] * z__[ip1];
            b = z__[ip1] * z__[ip1] * delsq;
            if (a < 0.) {
                tau2 = b * 2. / (a - sqrt((d__1 = a * a + b * 4. * c__, abs(d__1))));
            } else {
                tau2 = -(a + sqrt((d__1 = a * a + b * 4. * c__, abs(d__1)))) / (c__ * 2.);
            }
            tau = tau2 / (d__[ip1] + sqrt((d__1 = d__[ip1] * d__[ip1] + tau2, abs(d__1))));
        }
        *sigma = d__[ii] + tau;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            work[j] = d__[j] + d__[ii] + tau;
            delta[j] = d__[j] - d__[ii] - tau;
        }
        iim1 = ii - 1;
        iip1 = ii + 1;
        dpsi = 0.;
        psi = 0.;
        erretm = 0.;
        i__1 = iim1;
        for (j = 1; j <= i__1; ++j) {
            temp = z__[j] / (work[j] * delta[j]);
            psi += z__[j] * temp;
            dpsi += temp * temp;
            erretm += psi;
        }
        erretm = abs(erretm);
        dphi = 0.;
        phi = 0.;
        i__1 = iip1;
        for (j = *n; j >= i__1; --j) {
            temp = z__[j] / (work[j] * delta[j]);
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
        temp = z__[ii] / (work[ii] * delta[ii]);
        dw = dpsi + dphi + temp * temp;
        temp = z__[ii] * temp;
        w += temp;
        erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + abs(temp) * 3.;
        if (abs(w) <= eps * erretm) {
            goto L240;
        }
        if (w <= 0.) {
            sglb = max(sglb, tau);
        } else {
            sgub = min(sgub, tau);
        }
        ++niter;
        if (!swtch3) {
            dtipsq = work[ip1] * delta[ip1];
            dtisq = work[*i__] * delta[*i__];
            if (orgati) {
                d__1 = z__[*i__] / dtisq;
                c__ = w - dtipsq * dw + delsq * (d__1 * d__1);
            } else {
                d__1 = z__[ip1] / dtipsq;
                c__ = w - dtisq * dw - delsq * (d__1 * d__1);
            }
            a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
            b = dtipsq * dtisq * w;
            if (c__ == 0.) {
                if (a == 0.) {
                    if (orgati) {
                        a = z__[*i__] * z__[*i__] + dtipsq * dtipsq * (dpsi + dphi);
                    } else {
                        a = z__[ip1] * z__[ip1] + dtisq * dtisq * (dpsi + dphi);
                    }
                }
                eta = b / a;
            } else if (a <= 0.) {
                eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (c__ * 2.);
            } else {
                eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))));
            }
        } else {
            dtiim = work[iim1] * delta[iim1];
            dtiip = work[iip1] * delta[iip1];
            temp = rhoinv + psi + phi;
            if (orgati) {
                temp1 = z__[iim1] / dtiim;
                temp1 *= temp1;
                c__ = temp - dtiip * (dpsi + dphi) -
                      (d__[iim1] - d__[iip1]) * (d__[iim1] + d__[iip1]) * temp1;
                zz[0] = z__[iim1] * z__[iim1];
                if (dpsi < temp1) {
                    zz[2] = dtiip * dtiip * dphi;
                } else {
                    zz[2] = dtiip * dtiip * (dpsi - temp1 + dphi);
                }
            } else {
                temp1 = z__[iip1] / dtiip;
                temp1 *= temp1;
                c__ = temp - dtiim * (dpsi + dphi) -
                      (d__[iip1] - d__[iim1]) * (d__[iim1] + d__[iip1]) * temp1;
                if (dphi < temp1) {
                    zz[0] = dtiim * dtiim * dpsi;
                } else {
                    zz[0] = dtiim * dtiim * (dpsi + (dphi - temp1));
                }
                zz[2] = z__[iip1] * z__[iip1];
            }
            zz[1] = z__[ii] * z__[ii];
            dd[0] = dtiim;
            dd[1] = delta[ii] * work[ii];
            dd[2] = dtiip;
            dlaed6_(&niter, &orgati, &c__, dd, zz, &w, &eta, info);
            if (*info != 0) {
                swtch3 = FALSE_;
                *info = 0;
                dtipsq = work[ip1] * delta[ip1];
                dtisq = work[*i__] * delta[*i__];
                if (orgati) {
                    d__1 = z__[*i__] / dtisq;
                    c__ = w - dtipsq * dw + delsq * (d__1 * d__1);
                } else {
                    d__1 = z__[ip1] / dtipsq;
                    c__ = w - dtisq * dw - delsq * (d__1 * d__1);
                }
                a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
                b = dtipsq * dtisq * w;
                if (c__ == 0.) {
                    if (a == 0.) {
                        if (orgati) {
                            a = z__[*i__] * z__[*i__] + dtipsq * dtipsq * (dpsi + dphi);
                        } else {
                            a = z__[ip1] * z__[ip1] + dtisq * dtisq * (dpsi + dphi);
                        }
                    }
                    eta = b / a;
                } else if (a <= 0.) {
                    eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (c__ * 2.);
                } else {
                    eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))));
                }
            }
        }
        if (w * eta >= 0.) {
            eta = -w / dw;
        }
        eta /= *sigma + sqrt(*sigma * *sigma + eta);
        temp = tau + eta;
        if (temp > sgub || temp < sglb) {
            if (w < 0.) {
                eta = (sgub - tau) / 2.;
            } else {
                eta = (sglb - tau) / 2.;
            }
            if (geomavg) {
                if (w < 0.) {
                    if (tau > 0.) {
                        eta = sqrt(sgub * tau) - tau;
                    }
                } else {
                    if (sglb > 0.) {
                        eta = sqrt(sglb * tau) - tau;
                    }
                }
            }
        }
        prew = w;
        tau += eta;
        *sigma += eta;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            work[j] += eta;
            delta[j] -= eta;
        }
        dpsi = 0.;
        psi = 0.;
        erretm = 0.;
        i__1 = iim1;
        for (j = 1; j <= i__1; ++j) {
            temp = z__[j] / (work[j] * delta[j]);
            psi += z__[j] * temp;
            dpsi += temp * temp;
            erretm += psi;
        }
        erretm = abs(erretm);
        dphi = 0.;
        phi = 0.;
        i__1 = iip1;
        for (j = *n; j >= i__1; --j) {
            temp = z__[j] / (work[j] * delta[j]);
            phi += z__[j] * temp;
            dphi += temp * temp;
            erretm += phi;
        }
        tau2 = work[ii] * delta[ii];
        temp = z__[ii] / tau2;
        dw = dpsi + dphi + temp * temp;
        temp = z__[ii] * temp;
        w = rhoinv + phi + psi + temp;
        erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + abs(temp) * 3.;
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
        iter = niter + 1;
        for (niter = iter; niter <= 400; ++niter) {
            if (abs(w) <= eps * erretm) {
                goto L240;
            }
            if (w <= 0.) {
                sglb = max(sglb, tau);
            } else {
                sgub = min(sgub, tau);
            }
            if (!swtch3) {
                dtipsq = work[ip1] * delta[ip1];
                dtisq = work[*i__] * delta[*i__];
                if (!swtch) {
                    if (orgati) {
                        d__1 = z__[*i__] / dtisq;
                        c__ = w - dtipsq * dw + delsq * (d__1 * d__1);
                    } else {
                        d__1 = z__[ip1] / dtipsq;
                        c__ = w - dtisq * dw - delsq * (d__1 * d__1);
                    }
                } else {
                    temp = z__[ii] / (work[ii] * delta[ii]);
                    if (orgati) {
                        dpsi += temp * temp;
                    } else {
                        dphi += temp * temp;
                    }
                    c__ = w - dtisq * dpsi - dtipsq * dphi;
                }
                a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
                b = dtipsq * dtisq * w;
                if (c__ == 0.) {
                    if (a == 0.) {
                        if (!swtch) {
                            if (orgati) {
                                a = z__[*i__] * z__[*i__] + dtipsq * dtipsq * (dpsi + dphi);
                            } else {
                                a = z__[ip1] * z__[ip1] + dtisq * dtisq * (dpsi + dphi);
                            }
                        } else {
                            a = dtisq * dtisq * dpsi + dtipsq * dtipsq * dphi;
                        }
                    }
                    eta = b / a;
                } else if (a <= 0.) {
                    eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (c__ * 2.);
                } else {
                    eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))));
                }
            } else {
                dtiim = work[iim1] * delta[iim1];
                dtiip = work[iip1] * delta[iip1];
                temp = rhoinv + psi + phi;
                if (swtch) {
                    c__ = temp - dtiim * dpsi - dtiip * dphi;
                    zz[0] = dtiim * dtiim * dpsi;
                    zz[2] = dtiip * dtiip * dphi;
                } else {
                    if (orgati) {
                        temp1 = z__[iim1] / dtiim;
                        temp1 *= temp1;
                        temp2 = (d__[iim1] - d__[iip1]) * (d__[iim1] + d__[iip1]) * temp1;
                        c__ = temp - dtiip * (dpsi + dphi) - temp2;
                        zz[0] = z__[iim1] * z__[iim1];
                        if (dpsi < temp1) {
                            zz[2] = dtiip * dtiip * dphi;
                        } else {
                            zz[2] = dtiip * dtiip * (dpsi - temp1 + dphi);
                        }
                    } else {
                        temp1 = z__[iip1] / dtiip;
                        temp1 *= temp1;
                        temp2 = (d__[iip1] - d__[iim1]) * (d__[iim1] + d__[iip1]) * temp1;
                        c__ = temp - dtiim * (dpsi + dphi) - temp2;
                        if (dphi < temp1) {
                            zz[0] = dtiim * dtiim * dpsi;
                        } else {
                            zz[0] = dtiim * dtiim * (dpsi + (dphi - temp1));
                        }
                        zz[2] = z__[iip1] * z__[iip1];
                    }
                }
                dd[0] = dtiim;
                dd[1] = delta[ii] * work[ii];
                dd[2] = dtiip;
                dlaed6_(&niter, &orgati, &c__, dd, zz, &w, &eta, info);
                if (*info != 0) {
                    swtch3 = FALSE_;
                    *info = 0;
                    dtipsq = work[ip1] * delta[ip1];
                    dtisq = work[*i__] * delta[*i__];
                    if (!swtch) {
                        if (orgati) {
                            d__1 = z__[*i__] / dtisq;
                            c__ = w - dtipsq * dw + delsq * (d__1 * d__1);
                        } else {
                            d__1 = z__[ip1] / dtipsq;
                            c__ = w - dtisq * dw - delsq * (d__1 * d__1);
                        }
                    } else {
                        temp = z__[ii] / (work[ii] * delta[ii]);
                        if (orgati) {
                            dpsi += temp * temp;
                        } else {
                            dphi += temp * temp;
                        }
                        c__ = w - dtisq * dpsi - dtipsq * dphi;
                    }
                    a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
                    b = dtipsq * dtisq * w;
                    if (c__ == 0.) {
                        if (a == 0.) {
                            if (!swtch) {
                                if (orgati) {
                                    a = z__[*i__] * z__[*i__] + dtipsq * dtipsq * (dpsi + dphi);
                                } else {
                                    a = z__[ip1] * z__[ip1] + dtisq * dtisq * (dpsi + dphi);
                                }
                            } else {
                                a = dtisq * dtisq * dpsi + dtipsq * dtipsq * dphi;
                            }
                        }
                        eta = b / a;
                    } else if (a <= 0.) {
                        eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (c__ * 2.);
                    } else {
                        eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))));
                    }
                }
            }
            if (w * eta >= 0.) {
                eta = -w / dw;
            }
            eta /= *sigma + sqrt(*sigma * *sigma + eta);
            temp = tau + eta;
            if (temp > sgub || temp < sglb) {
                if (w < 0.) {
                    eta = (sgub - tau) / 2.;
                } else {
                    eta = (sglb - tau) / 2.;
                }
                if (geomavg) {
                    if (w < 0.) {
                        if (tau > 0.) {
                            eta = sqrt(sgub * tau) - tau;
                        }
                    } else {
                        if (sglb > 0.) {
                            eta = sqrt(sglb * tau) - tau;
                        }
                    }
                }
            }
            prew = w;
            tau += eta;
            *sigma += eta;
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                work[j] += eta;
                delta[j] -= eta;
            }
            dpsi = 0.;
            psi = 0.;
            erretm = 0.;
            i__1 = iim1;
            for (j = 1; j <= i__1; ++j) {
                temp = z__[j] / (work[j] * delta[j]);
                psi += z__[j] * temp;
                dpsi += temp * temp;
                erretm += psi;
            }
            erretm = abs(erretm);
            dphi = 0.;
            phi = 0.;
            i__1 = iip1;
            for (j = *n; j >= i__1; --j) {
                temp = z__[j] / (work[j] * delta[j]);
                phi += z__[j] * temp;
                dphi += temp * temp;
                erretm += phi;
            }
            tau2 = work[ii] * delta[ii];
            temp = z__[ii] / tau2;
            dw = dpsi + dphi + temp * temp;
            temp = z__[ii] * temp;
            w = rhoinv + phi + psi + temp;
            erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + abs(temp) * 3.;
            if (w * prew > 0. && abs(w) > abs(prew) / 10.) {
                swtch = !swtch;
            }
        }
        *info = 1;
    }
L240:
    return 0;
}
#ifdef __cplusplus
}
#endif
