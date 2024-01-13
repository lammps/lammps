#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static doublereal c_b46 = .5;
int dlatrs_(char *uplo, char *trans, char *diag, char *normin, integer *n, doublereal *a,
            integer *lda, doublereal *x, doublereal *scale, doublereal *cnorm, integer *info,
            ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len, ftnlen normin_len)
{
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;
    integer i__, j;
    doublereal xj, rec, tjj;
    integer jinc;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, integer *);
    doublereal xbnd;
    integer imax;
    doublereal tmax, tjjs, xmax, grow, sumj, work[1];
    extern int dscal_(integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    doublereal tscal, uscal;
    extern doublereal dasum_(integer *, doublereal *, integer *);
    integer jlast;
    extern int daxpy_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *);
    logical upper;
    extern int dtrsv_(char *, char *, char *, integer *, doublereal *, integer *, doublereal *,
                      integer *, ftnlen, ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen),
        dlange_(char *, integer *, integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern int xerbla_(char *, integer *, ftnlen);
    doublereal bignum;
    logical notran;
    integer jfirst;
    doublereal smlnum;
    logical nounit;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    --cnorm;
    *info = 0;
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    notran = lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1);
    nounit = lsame_(diag, (char *)"N", (ftnlen)1, (ftnlen)1);
    if (!upper && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (!notran && !lsame_(trans, (char *)"T", (ftnlen)1, (ftnlen)1) &&
               !lsame_(trans, (char *)"C", (ftnlen)1, (ftnlen)1)) {
        *info = -2;
    } else if (!nounit && !lsame_(diag, (char *)"U", (ftnlen)1, (ftnlen)1)) {
        *info = -3;
    } else if (!lsame_(normin, (char *)"Y", (ftnlen)1, (ftnlen)1) &&
               !lsame_(normin, (char *)"N", (ftnlen)1, (ftnlen)1)) {
        *info = -4;
    } else if (*n < 0) {
        *info = -5;
    } else if (*lda < max(1, *n)) {
        *info = -7;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLATRS", &i__1, (ftnlen)6);
        return 0;
    }
    *scale = 1.;
    if (*n == 0) {
        return 0;
    }
    smlnum = dlamch_((char *)"Safe minimum", (ftnlen)12) / dlamch_((char *)"Precision", (ftnlen)9);
    bignum = 1. / smlnum;
    if (lsame_(normin, (char *)"N", (ftnlen)1, (ftnlen)1)) {
        if (upper) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = j - 1;
                cnorm[j] = dasum_(&i__2, &a[j * a_dim1 + 1], &c__1);
            }
        } else {
            i__1 = *n - 1;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *n - j;
                cnorm[j] = dasum_(&i__2, &a[j + 1 + j * a_dim1], &c__1);
            }
            cnorm[*n] = 0.;
        }
    }
    imax = idamax_(n, &cnorm[1], &c__1);
    tmax = cnorm[imax];
    if (tmax <= bignum) {
        tscal = 1.;
    } else {
        if (tmax <= dlamch_((char *)"Overflow", (ftnlen)8)) {
            tscal = 1. / (smlnum * tmax);
            dscal_(n, &tscal, &cnorm[1], &c__1);
        } else {
            tmax = 0.;
            if (upper) {
                i__1 = *n;
                for (j = 2; j <= i__1; ++j) {
                    i__2 = j - 1;
                    d__1 = dlange_((char *)"M", &i__2, &c__1, &a[j * a_dim1 + 1], &c__1, work, (ftnlen)1);
                    tmax = max(d__1, tmax);
                }
            } else {
                i__1 = *n - 1;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *n - j;
                    d__1 =
                        dlange_((char *)"M", &i__2, &c__1, &a[j + 1 + j * a_dim1], &c__1, work, (ftnlen)1);
                    tmax = max(d__1, tmax);
                }
            }
            if (tmax <= dlamch_((char *)"Overflow", (ftnlen)8)) {
                tscal = 1. / (smlnum * tmax);
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    if (cnorm[j] <= dlamch_((char *)"Overflow", (ftnlen)8)) {
                        cnorm[j] *= tscal;
                    } else {
                        cnorm[j] = 0.;
                        if (upper) {
                            i__2 = j - 1;
                            for (i__ = 1; i__ <= i__2; ++i__) {
                                cnorm[j] += tscal * (d__1 = a[i__ + j * a_dim1], abs(d__1));
                            }
                        } else {
                            i__2 = *n;
                            for (i__ = j + 1; i__ <= i__2; ++i__) {
                                cnorm[j] += tscal * (d__1 = a[i__ + j * a_dim1], abs(d__1));
                            }
                        }
                    }
                }
            } else {
                dtrsv_(uplo, trans, diag, n, &a[a_offset], lda, &x[1], &c__1, (ftnlen)1, (ftnlen)1,
                       (ftnlen)1);
                return 0;
            }
        }
    }
    j = idamax_(n, &x[1], &c__1);
    xmax = (d__1 = x[j], abs(d__1));
    xbnd = xmax;
    if (notran) {
        if (upper) {
            jfirst = *n;
            jlast = 1;
            jinc = -1;
        } else {
            jfirst = 1;
            jlast = *n;
            jinc = 1;
        }
        if (tscal != 1.) {
            grow = 0.;
            goto L50;
        }
        if (nounit) {
            grow = 1. / max(xbnd, smlnum);
            xbnd = grow;
            i__1 = jlast;
            i__2 = jinc;
            for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
                if (grow <= smlnum) {
                    goto L50;
                }
                tjj = (d__1 = a[j + j * a_dim1], abs(d__1));
                d__1 = xbnd, d__2 = min(1., tjj) * grow;
                xbnd = min(d__1, d__2);
                if (tjj + cnorm[j] >= smlnum) {
                    grow *= tjj / (tjj + cnorm[j]);
                } else {
                    grow = 0.;
                }
            }
            grow = xbnd;
        } else {
            d__1 = 1., d__2 = 1. / max(xbnd, smlnum);
            grow = min(d__1, d__2);
            i__2 = jlast;
            i__1 = jinc;
            for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {
                if (grow <= smlnum) {
                    goto L50;
                }
                grow *= 1. / (cnorm[j] + 1.);
            }
        }
    L50:;
    } else {
        if (upper) {
            jfirst = 1;
            jlast = *n;
            jinc = 1;
        } else {
            jfirst = *n;
            jlast = 1;
            jinc = -1;
        }
        if (tscal != 1.) {
            grow = 0.;
            goto L80;
        }
        if (nounit) {
            grow = 1. / max(xbnd, smlnum);
            xbnd = grow;
            i__1 = jlast;
            i__2 = jinc;
            for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
                if (grow <= smlnum) {
                    goto L80;
                }
                xj = cnorm[j] + 1.;
                d__1 = grow, d__2 = xbnd / xj;
                grow = min(d__1, d__2);
                tjj = (d__1 = a[j + j * a_dim1], abs(d__1));
                if (xj > tjj) {
                    xbnd *= tjj / xj;
                }
            }
            grow = min(grow, xbnd);
        } else {
            d__1 = 1., d__2 = 1. / max(xbnd, smlnum);
            grow = min(d__1, d__2);
            i__2 = jlast;
            i__1 = jinc;
            for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {
                if (grow <= smlnum) {
                    goto L80;
                }
                xj = cnorm[j] + 1.;
                grow /= xj;
            }
        }
    L80:;
    }
    if (grow * tscal > smlnum) {
        dtrsv_(uplo, trans, diag, n, &a[a_offset], lda, &x[1], &c__1, (ftnlen)1, (ftnlen)1,
               (ftnlen)1);
    } else {
        if (xmax > bignum) {
            *scale = bignum / xmax;
            dscal_(n, scale, &x[1], &c__1);
            xmax = bignum;
        }
        if (notran) {
            i__1 = jlast;
            i__2 = jinc;
            for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
                xj = (d__1 = x[j], abs(d__1));
                if (nounit) {
                    tjjs = a[j + j * a_dim1] * tscal;
                } else {
                    tjjs = tscal;
                    if (tscal == 1.) {
                        goto L100;
                    }
                }
                tjj = abs(tjjs);
                if (tjj > smlnum) {
                    if (tjj < 1.) {
                        if (xj > tjj * bignum) {
                            rec = 1. / xj;
                            dscal_(n, &rec, &x[1], &c__1);
                            *scale *= rec;
                            xmax *= rec;
                        }
                    }
                    x[j] /= tjjs;
                    xj = (d__1 = x[j], abs(d__1));
                } else if (tjj > 0.) {
                    if (xj > tjj * bignum) {
                        rec = tjj * bignum / xj;
                        if (cnorm[j] > 1.) {
                            rec /= cnorm[j];
                        }
                        dscal_(n, &rec, &x[1], &c__1);
                        *scale *= rec;
                        xmax *= rec;
                    }
                    x[j] /= tjjs;
                    xj = (d__1 = x[j], abs(d__1));
                } else {
                    i__3 = *n;
                    for (i__ = 1; i__ <= i__3; ++i__) {
                        x[i__] = 0.;
                    }
                    x[j] = 1.;
                    xj = 1.;
                    *scale = 0.;
                    xmax = 0.;
                }
            L100:
                if (xj > 1.) {
                    rec = 1. / xj;
                    if (cnorm[j] > (bignum - xmax) * rec) {
                        rec *= .5;
                        dscal_(n, &rec, &x[1], &c__1);
                        *scale *= rec;
                    }
                } else if (xj * cnorm[j] > bignum - xmax) {
                    dscal_(n, &c_b46, &x[1], &c__1);
                    *scale *= .5;
                }
                if (upper) {
                    if (j > 1) {
                        i__3 = j - 1;
                        d__1 = -x[j] * tscal;
                        daxpy_(&i__3, &d__1, &a[j * a_dim1 + 1], &c__1, &x[1], &c__1);
                        i__3 = j - 1;
                        i__ = idamax_(&i__3, &x[1], &c__1);
                        xmax = (d__1 = x[i__], abs(d__1));
                    }
                } else {
                    if (j < *n) {
                        i__3 = *n - j;
                        d__1 = -x[j] * tscal;
                        daxpy_(&i__3, &d__1, &a[j + 1 + j * a_dim1], &c__1, &x[j + 1], &c__1);
                        i__3 = *n - j;
                        i__ = j + idamax_(&i__3, &x[j + 1], &c__1);
                        xmax = (d__1 = x[i__], abs(d__1));
                    }
                }
            }
        } else {
            i__2 = jlast;
            i__1 = jinc;
            for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {
                xj = (d__1 = x[j], abs(d__1));
                uscal = tscal;
                rec = 1. / max(xmax, 1.);
                if (cnorm[j] > (bignum - xj) * rec) {
                    rec *= .5;
                    if (nounit) {
                        tjjs = a[j + j * a_dim1] * tscal;
                    } else {
                        tjjs = tscal;
                    }
                    tjj = abs(tjjs);
                    if (tjj > 1.) {
                        d__1 = 1., d__2 = rec * tjj;
                        rec = min(d__1, d__2);
                        uscal /= tjjs;
                    }
                    if (rec < 1.) {
                        dscal_(n, &rec, &x[1], &c__1);
                        *scale *= rec;
                        xmax *= rec;
                    }
                }
                sumj = 0.;
                if (uscal == 1.) {
                    if (upper) {
                        i__3 = j - 1;
                        sumj = ddot_(&i__3, &a[j * a_dim1 + 1], &c__1, &x[1], &c__1);
                    } else if (j < *n) {
                        i__3 = *n - j;
                        sumj = ddot_(&i__3, &a[j + 1 + j * a_dim1], &c__1, &x[j + 1], &c__1);
                    }
                } else {
                    if (upper) {
                        i__3 = j - 1;
                        for (i__ = 1; i__ <= i__3; ++i__) {
                            sumj += a[i__ + j * a_dim1] * uscal * x[i__];
                        }
                    } else if (j < *n) {
                        i__3 = *n;
                        for (i__ = j + 1; i__ <= i__3; ++i__) {
                            sumj += a[i__ + j * a_dim1] * uscal * x[i__];
                        }
                    }
                }
                if (uscal == tscal) {
                    x[j] -= sumj;
                    xj = (d__1 = x[j], abs(d__1));
                    if (nounit) {
                        tjjs = a[j + j * a_dim1] * tscal;
                    } else {
                        tjjs = tscal;
                        if (tscal == 1.) {
                            goto L150;
                        }
                    }
                    tjj = abs(tjjs);
                    if (tjj > smlnum) {
                        if (tjj < 1.) {
                            if (xj > tjj * bignum) {
                                rec = 1. / xj;
                                dscal_(n, &rec, &x[1], &c__1);
                                *scale *= rec;
                                xmax *= rec;
                            }
                        }
                        x[j] /= tjjs;
                    } else if (tjj > 0.) {
                        if (xj > tjj * bignum) {
                            rec = tjj * bignum / xj;
                            dscal_(n, &rec, &x[1], &c__1);
                            *scale *= rec;
                            xmax *= rec;
                        }
                        x[j] /= tjjs;
                    } else {
                        i__3 = *n;
                        for (i__ = 1; i__ <= i__3; ++i__) {
                            x[i__] = 0.;
                        }
                        x[j] = 1.;
                        *scale = 0.;
                        xmax = 0.;
                    }
                L150:;
                } else {
                    x[j] = x[j] / tjjs - sumj;
                }
                d__2 = xmax, d__3 = (d__1 = x[j], abs(d__1));
                xmax = max(d__2, d__3);
            }
        }
        *scale /= tscal;
    }
    if (tscal != 1.) {
        d__1 = 1. / tscal;
        dscal_(n, &d__1, &cnorm[1], &c__1);
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
