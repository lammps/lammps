#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlaln2_(logical *ltrans, integer *na, integer *nw, doublereal *smin, doublereal *ca,
            doublereal *a, integer *lda, doublereal *d1, doublereal *d2, doublereal *b,
            integer *ldb, doublereal *wr, doublereal *wi, doublereal *x, integer *ldx,
            doublereal *scale, doublereal *xnorm, integer *info)
{
    static logical zswap[4] = {FALSE_, FALSE_, TRUE_, TRUE_};
    static logical rswap[4] = {FALSE_, TRUE_, FALSE_, TRUE_};
    static integer ipivot[16] = {1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4, 3, 2, 1};
    integer a_dim1, a_offset, b_dim1, b_offset, x_dim1, x_offset;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;
    static doublereal equiv_0[4], equiv_1[4];
    integer j;
#define ci (equiv_0)
#define cr (equiv_1)
    doublereal bi1, bi2, br1, br2, xi1, xi2, xr1, xr2, ci21, ci22, cr21, cr22, li21, csi, ui11,
        lr21, ui12, ui22;
#define civ (equiv_0)
    doublereal csr, ur11, ur12, ur22;
#define crv (equiv_1)
    doublereal bbnd, cmax, ui11r, ui12s, temp, ur11r, ur12s, u22abs;
    integer icmax;
    doublereal bnorm, cnorm, smini;
    extern doublereal dlamch_(char *, ftnlen);
    extern int dladiv_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                       doublereal *);
    doublereal bignum, smlnum;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    smlnum = 2. * dlamch_((char *)"Safe minimum", (ftnlen)12);
    bignum = 1. / smlnum;
    smini = max(*smin, smlnum);
    *info = 0;
    *scale = 1.;
    if (*na == 1) {
        if (*nw == 1) {
            csr = *ca * a[a_dim1 + 1] - *wr * *d1;
            cnorm = abs(csr);
            if (cnorm < smini) {
                csr = smini;
                cnorm = smini;
                *info = 1;
            }
            bnorm = (d__1 = b[b_dim1 + 1], abs(d__1));
            if (cnorm < 1. && bnorm > 1.) {
                if (bnorm > bignum * cnorm) {
                    *scale = 1. / bnorm;
                }
            }
            x[x_dim1 + 1] = b[b_dim1 + 1] * *scale / csr;
            *xnorm = (d__1 = x[x_dim1 + 1], abs(d__1));
        } else {
            csr = *ca * a[a_dim1 + 1] - *wr * *d1;
            csi = -(*wi) * *d1;
            cnorm = abs(csr) + abs(csi);
            if (cnorm < smini) {
                csr = smini;
                csi = 0.;
                cnorm = smini;
                *info = 1;
            }
            bnorm = (d__1 = b[b_dim1 + 1], abs(d__1)) + (d__2 = b[(b_dim1 << 1) + 1], abs(d__2));
            if (cnorm < 1. && bnorm > 1.) {
                if (bnorm > bignum * cnorm) {
                    *scale = 1. / bnorm;
                }
            }
            d__1 = *scale * b[b_dim1 + 1];
            d__2 = *scale * b[(b_dim1 << 1) + 1];
            dladiv_(&d__1, &d__2, &csr, &csi, &x[x_dim1 + 1], &x[(x_dim1 << 1) + 1]);
            *xnorm = (d__1 = x[x_dim1 + 1], abs(d__1)) + (d__2 = x[(x_dim1 << 1) + 1], abs(d__2));
        }
    } else {
        cr[0] = *ca * a[a_dim1 + 1] - *wr * *d1;
        cr[3] = *ca * a[(a_dim1 << 1) + 2] - *wr * *d2;
        if (*ltrans) {
            cr[2] = *ca * a[a_dim1 + 2];
            cr[1] = *ca * a[(a_dim1 << 1) + 1];
        } else {
            cr[1] = *ca * a[a_dim1 + 2];
            cr[2] = *ca * a[(a_dim1 << 1) + 1];
        }
        if (*nw == 1) {
            cmax = 0.;
            icmax = 0;
            for (j = 1; j <= 4; ++j) {
                if ((d__1 = crv[j - 1], abs(d__1)) > cmax) {
                    cmax = (d__1 = crv[j - 1], abs(d__1));
                    icmax = j;
                }
            }
            if (cmax < smini) {
                d__3 = (d__1 = b[b_dim1 + 1], abs(d__1)), d__4 = (d__2 = b[b_dim1 + 2], abs(d__2));
                bnorm = max(d__3, d__4);
                if (smini < 1. && bnorm > 1.) {
                    if (bnorm > bignum * smini) {
                        *scale = 1. / bnorm;
                    }
                }
                temp = *scale / smini;
                x[x_dim1 + 1] = temp * b[b_dim1 + 1];
                x[x_dim1 + 2] = temp * b[b_dim1 + 2];
                *xnorm = temp * bnorm;
                *info = 1;
                return 0;
            }
            ur11 = crv[icmax - 1];
            cr21 = crv[ipivot[(icmax << 2) - 3] - 1];
            ur12 = crv[ipivot[(icmax << 2) - 2] - 1];
            cr22 = crv[ipivot[(icmax << 2) - 1] - 1];
            ur11r = 1. / ur11;
            lr21 = ur11r * cr21;
            ur22 = cr22 - ur12 * lr21;
            if (abs(ur22) < smini) {
                ur22 = smini;
                *info = 1;
            }
            if (rswap[icmax - 1]) {
                br1 = b[b_dim1 + 2];
                br2 = b[b_dim1 + 1];
            } else {
                br1 = b[b_dim1 + 1];
                br2 = b[b_dim1 + 2];
            }
            br2 -= lr21 * br1;
            d__2 = (d__1 = br1 * (ur22 * ur11r), abs(d__1)), d__3 = abs(br2);
            bbnd = max(d__2, d__3);
            if (bbnd > 1. && abs(ur22) < 1.) {
                if (bbnd >= bignum * abs(ur22)) {
                    *scale = 1. / bbnd;
                }
            }
            xr2 = br2 * *scale / ur22;
            xr1 = *scale * br1 * ur11r - xr2 * (ur11r * ur12);
            if (zswap[icmax - 1]) {
                x[x_dim1 + 1] = xr2;
                x[x_dim1 + 2] = xr1;
            } else {
                x[x_dim1 + 1] = xr1;
                x[x_dim1 + 2] = xr2;
            }
            d__1 = abs(xr1), d__2 = abs(xr2);
            *xnorm = max(d__1, d__2);
            if (*xnorm > 1. && cmax > 1.) {
                if (*xnorm > bignum / cmax) {
                    temp = cmax / bignum;
                    x[x_dim1 + 1] = temp * x[x_dim1 + 1];
                    x[x_dim1 + 2] = temp * x[x_dim1 + 2];
                    *xnorm = temp * *xnorm;
                    *scale = temp * *scale;
                }
            }
        } else {
            ci[0] = -(*wi) * *d1;
            ci[1] = 0.;
            ci[2] = 0.;
            ci[3] = -(*wi) * *d2;
            cmax = 0.;
            icmax = 0;
            for (j = 1; j <= 4; ++j) {
                if ((d__1 = crv[j - 1], abs(d__1)) + (d__2 = civ[j - 1], abs(d__2)) > cmax) {
                    cmax = (d__1 = crv[j - 1], abs(d__1)) + (d__2 = civ[j - 1], abs(d__2));
                    icmax = j;
                }
            }
            if (cmax < smini) {
                d__5 = (d__1 = b[b_dim1 + 1], abs(d__1)) + (d__2 = b[(b_dim1 << 1) + 1], abs(d__2)),
                d__6 = (d__3 = b[b_dim1 + 2], abs(d__3)) + (d__4 = b[(b_dim1 << 1) + 2], abs(d__4));
                bnorm = max(d__5, d__6);
                if (smini < 1. && bnorm > 1.) {
                    if (bnorm > bignum * smini) {
                        *scale = 1. / bnorm;
                    }
                }
                temp = *scale / smini;
                x[x_dim1 + 1] = temp * b[b_dim1 + 1];
                x[x_dim1 + 2] = temp * b[b_dim1 + 2];
                x[(x_dim1 << 1) + 1] = temp * b[(b_dim1 << 1) + 1];
                x[(x_dim1 << 1) + 2] = temp * b[(b_dim1 << 1) + 2];
                *xnorm = temp * bnorm;
                *info = 1;
                return 0;
            }
            ur11 = crv[icmax - 1];
            ui11 = civ[icmax - 1];
            cr21 = crv[ipivot[(icmax << 2) - 3] - 1];
            ci21 = civ[ipivot[(icmax << 2) - 3] - 1];
            ur12 = crv[ipivot[(icmax << 2) - 2] - 1];
            ui12 = civ[ipivot[(icmax << 2) - 2] - 1];
            cr22 = crv[ipivot[(icmax << 2) - 1] - 1];
            ci22 = civ[ipivot[(icmax << 2) - 1] - 1];
            if (icmax == 1 || icmax == 4) {
                if (abs(ur11) > abs(ui11)) {
                    temp = ui11 / ur11;
                    d__1 = temp;
                    ur11r = 1. / (ur11 * (d__1 * d__1 + 1.));
                    ui11r = -temp * ur11r;
                } else {
                    temp = ur11 / ui11;
                    d__1 = temp;
                    ui11r = -1. / (ui11 * (d__1 * d__1 + 1.));
                    ur11r = -temp * ui11r;
                }
                lr21 = cr21 * ur11r;
                li21 = cr21 * ui11r;
                ur12s = ur12 * ur11r;
                ui12s = ur12 * ui11r;
                ur22 = cr22 - ur12 * lr21;
                ui22 = ci22 - ur12 * li21;
            } else {
                ur11r = 1. / ur11;
                ui11r = 0.;
                lr21 = cr21 * ur11r;
                li21 = ci21 * ur11r;
                ur12s = ur12 * ur11r;
                ui12s = ui12 * ur11r;
                ur22 = cr22 - ur12 * lr21 + ui12 * li21;
                ui22 = -ur12 * li21 - ui12 * lr21;
            }
            u22abs = abs(ur22) + abs(ui22);
            if (u22abs < smini) {
                ur22 = smini;
                ui22 = 0.;
                *info = 1;
            }
            if (rswap[icmax - 1]) {
                br2 = b[b_dim1 + 1];
                br1 = b[b_dim1 + 2];
                bi2 = b[(b_dim1 << 1) + 1];
                bi1 = b[(b_dim1 << 1) + 2];
            } else {
                br1 = b[b_dim1 + 1];
                br2 = b[b_dim1 + 2];
                bi1 = b[(b_dim1 << 1) + 1];
                bi2 = b[(b_dim1 << 1) + 2];
            }
            br2 = br2 - lr21 * br1 + li21 * bi1;
            bi2 = bi2 - li21 * br1 - lr21 * bi1;
            d__1 = (abs(br1) + abs(bi1)) * (u22abs * (abs(ur11r) + abs(ui11r))),
            d__2 = abs(br2) + abs(bi2);
            bbnd = max(d__1, d__2);
            if (bbnd > 1. && u22abs < 1.) {
                if (bbnd >= bignum * u22abs) {
                    *scale = 1. / bbnd;
                    br1 = *scale * br1;
                    bi1 = *scale * bi1;
                    br2 = *scale * br2;
                    bi2 = *scale * bi2;
                }
            }
            dladiv_(&br2, &bi2, &ur22, &ui22, &xr2, &xi2);
            xr1 = ur11r * br1 - ui11r * bi1 - ur12s * xr2 + ui12s * xi2;
            xi1 = ui11r * br1 + ur11r * bi1 - ui12s * xr2 - ur12s * xi2;
            if (zswap[icmax - 1]) {
                x[x_dim1 + 1] = xr2;
                x[x_dim1 + 2] = xr1;
                x[(x_dim1 << 1) + 1] = xi2;
                x[(x_dim1 << 1) + 2] = xi1;
            } else {
                x[x_dim1 + 1] = xr1;
                x[x_dim1 + 2] = xr2;
                x[(x_dim1 << 1) + 1] = xi1;
                x[(x_dim1 << 1) + 2] = xi2;
            }
            d__1 = abs(xr1) + abs(xi1), d__2 = abs(xr2) + abs(xi2);
            *xnorm = max(d__1, d__2);
            if (*xnorm > 1. && cmax > 1.) {
                if (*xnorm > bignum / cmax) {
                    temp = cmax / bignum;
                    x[x_dim1 + 1] = temp * x[x_dim1 + 1];
                    x[x_dim1 + 2] = temp * x[x_dim1 + 2];
                    x[(x_dim1 << 1) + 1] = temp * x[(x_dim1 << 1) + 1];
                    x[(x_dim1 << 1) + 2] = temp * x[(x_dim1 << 1) + 2];
                    *xnorm = temp * *xnorm;
                    *scale = temp * *scale;
                }
            }
        }
    }
    return 0;
}
#undef crv
#undef civ
#undef cr
#undef ci
#ifdef __cplusplus
}
#endif
