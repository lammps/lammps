#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static doublereal c_b17 = 0.;
static logical c_false = FALSE_;
static doublereal c_b29 = 1.;
static logical c_true = TRUE_;
int dtrevc3_(char *side, char *howmny, logical *select, integer *n, doublereal *t, integer *ldt,
             doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, integer *mm, integer *m,
             doublereal *work, integer *lwork, integer *info, ftnlen side_len, ftnlen howmny_len)
{
    address a__1[2];
    integer t_dim1, t_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1[2], i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;
    char ch__1[2];
    int s_lmp_cat(char *, char **, integer *, integer *, ftnlen);
    double sqrt(doublereal);
    integer i__, j, k;
    doublereal x[4];
    integer j1, j2, iscomplex[128], nb, ii, ki, ip, is, iv;
    doublereal wi, wr;
    integer ki2;
    doublereal rec, ulp, beta, emax;
    logical pair;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, integer *);
    logical allv;
    integer ierr;
    doublereal unfl, ovfl, smin;
    logical over;
    doublereal vmax;
    integer jnxt;
    extern int dscal_(integer *, doublereal *, doublereal *, integer *);
    doublereal scale;
    extern int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
                      integer *, doublereal *, integer *, doublereal *, doublereal *, integer *,
                      ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *,
                      doublereal *, integer *, doublereal *, doublereal *, integer *, ftnlen);
    doublereal remax;
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    logical leftv, bothv;
    extern int daxpy_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *);
    doublereal vcrit;
    logical somev;
    doublereal xnorm;
    extern int dlaln2_(logical *, integer *, integer *, doublereal *, doublereal *, doublereal *,
                       integer *, doublereal *, doublereal *, doublereal *, integer *, doublereal *,
                       doublereal *, doublereal *, integer *, doublereal *, doublereal *,
                       integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern int dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *,
                       integer *, ftnlen),
        xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    extern int dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *,
                       integer *, ftnlen);
    doublereal bignum;
    logical rightv;
    integer maxwrk;
    doublereal smlnum;
    logical lquery;
    --select;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    vl_dim1 = *ldvl;
    vl_offset = 1 + vl_dim1;
    vl -= vl_offset;
    vr_dim1 = *ldvr;
    vr_offset = 1 + vr_dim1;
    vr -= vr_offset;
    --work;
    bothv = lsame_(side, (char *)"B", (ftnlen)1, (ftnlen)1);
    rightv = lsame_(side, (char *)"R", (ftnlen)1, (ftnlen)1) || bothv;
    leftv = lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1) || bothv;
    allv = lsame_(howmny, (char *)"A", (ftnlen)1, (ftnlen)1);
    over = lsame_(howmny, (char *)"B", (ftnlen)1, (ftnlen)1);
    somev = lsame_(howmny, (char *)"S", (ftnlen)1, (ftnlen)1);
    *info = 0;
    i__1[0] = 1, a__1[0] = side;
    i__1[1] = 1, a__1[1] = howmny;
    s_lmp_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
    nb = ilaenv_(&c__1, (char *)"DTREVC", ch__1, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)2);
    i__2 = 1, i__3 = *n + (*n << 1) * nb;
    maxwrk = max(i__2, i__3);
    work[1] = (doublereal)maxwrk;
    lquery = *lwork == -1;
    if (!rightv && !leftv) {
        *info = -1;
    } else if (!allv && !over && !somev) {
        *info = -2;
    } else if (*n < 0) {
        *info = -4;
    } else if (*ldt < max(1, *n)) {
        *info = -6;
    } else if (*ldvl < 1 || leftv && *ldvl < *n) {
        *info = -8;
    } else if (*ldvr < 1 || rightv && *ldvr < *n) {
        *info = -10;
    } else {
        i__2 = 1, i__3 = *n * 3;
        if (*lwork < max(i__2, i__3) && !lquery) {
            *info = -14;
        } else {
            if (somev) {
                *m = 0;
                pair = FALSE_;
                i__2 = *n;
                for (j = 1; j <= i__2; ++j) {
                    if (pair) {
                        pair = FALSE_;
                        select[j] = FALSE_;
                    } else {
                        if (j < *n) {
                            if (t[j + 1 + j * t_dim1] == 0.) {
                                if (select[j]) {
                                    ++(*m);
                                }
                            } else {
                                pair = TRUE_;
                                if (select[j] || select[j + 1]) {
                                    select[j] = TRUE_;
                                    *m += 2;
                                }
                            }
                        } else {
                            if (select[*n]) {
                                ++(*m);
                            }
                        }
                    }
                }
            } else {
                *m = *n;
            }
            if (*mm < *m) {
                *info = -11;
            }
        }
    }
    if (*info != 0) {
        i__2 = -(*info);
        xerbla_((char *)"DTREVC3", &i__2, (ftnlen)7);
        return 0;
    } else if (lquery) {
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    if (over && *lwork >= *n + (*n << 4)) {
        nb = (*lwork - *n) / (*n << 1);
        nb = min(nb, 128);
        i__2 = (nb << 1) + 1;
        dlaset_((char *)"F", n, &i__2, &c_b17, &c_b17, &work[1], n, (ftnlen)1);
    } else {
        nb = 1;
    }
    unfl = dlamch_((char *)"Safe minimum", (ftnlen)12);
    ovfl = 1. / unfl;
    ulp = dlamch_((char *)"Precision", (ftnlen)9);
    smlnum = unfl * (*n / ulp);
    bignum = (1. - ulp) / smlnum;
    work[1] = 0.;
    i__2 = *n;
    for (j = 2; j <= i__2; ++j) {
        work[j] = 0.;
        i__3 = j - 1;
        for (i__ = 1; i__ <= i__3; ++i__) {
            work[j] += (d__1 = t[i__ + j * t_dim1], abs(d__1));
        }
    }
    if (rightv) {
        iv = 2;
        if (nb > 2) {
            iv = nb;
        }
        ip = 0;
        is = *m;
        for (ki = *n; ki >= 1; --ki) {
            if (ip == -1) {
                ip = 1;
                goto L140;
            } else if (ki == 1) {
                ip = 0;
            } else if (t[ki + (ki - 1) * t_dim1] == 0.) {
                ip = 0;
            } else {
                ip = -1;
            }
            if (somev) {
                if (ip == 0) {
                    if (!select[ki]) {
                        goto L140;
                    }
                } else {
                    if (!select[ki - 1]) {
                        goto L140;
                    }
                }
            }
            wr = t[ki + ki * t_dim1];
            wi = 0.;
            if (ip != 0) {
                wi = sqrt((d__1 = t[ki + (ki - 1) * t_dim1], abs(d__1))) *
                     sqrt((d__2 = t[ki - 1 + ki * t_dim1], abs(d__2)));
            }
            d__1 = ulp * (abs(wr) + abs(wi));
            smin = max(d__1, smlnum);
            if (ip == 0) {
                work[ki + iv * *n] = 1.;
                i__2 = ki - 1;
                for (k = 1; k <= i__2; ++k) {
                    work[k + iv * *n] = -t[k + ki * t_dim1];
                }
                jnxt = ki - 1;
                for (j = ki - 1; j >= 1; --j) {
                    if (j > jnxt) {
                        goto L60;
                    }
                    j1 = j;
                    j2 = j;
                    jnxt = j - 1;
                    if (j > 1) {
                        if (t[j + (j - 1) * t_dim1] != 0.) {
                            j1 = j - 1;
                            jnxt = j - 2;
                        }
                    }
                    if (j1 == j2) {
                        dlaln2_(&c_false, &c__1, &c__1, &smin, &c_b29, &t[j + j * t_dim1], ldt,
                                &c_b29, &c_b29, &work[j + iv * *n], n, &wr, &c_b17, x, &c__2,
                                &scale, &xnorm, &ierr);
                        if (xnorm > 1.) {
                            if (work[j] > bignum / xnorm) {
                                x[0] /= xnorm;
                                scale /= xnorm;
                            }
                        }
                        if (scale != 1.) {
                            dscal_(&ki, &scale, &work[iv * *n + 1], &c__1);
                        }
                        work[j + iv * *n] = x[0];
                        i__2 = j - 1;
                        d__1 = -x[0];
                        daxpy_(&i__2, &d__1, &t[j * t_dim1 + 1], &c__1, &work[iv * *n + 1], &c__1);
                    } else {
                        dlaln2_(&c_false, &c__2, &c__1, &smin, &c_b29, &t[j - 1 + (j - 1) * t_dim1],
                                ldt, &c_b29, &c_b29, &work[j - 1 + iv * *n], n, &wr, &c_b17, x,
                                &c__2, &scale, &xnorm, &ierr);
                        if (xnorm > 1.) {
                            d__1 = work[j - 1], d__2 = work[j];
                            beta = max(d__1, d__2);
                            if (beta > bignum / xnorm) {
                                x[0] /= xnorm;
                                x[1] /= xnorm;
                                scale /= xnorm;
                            }
                        }
                        if (scale != 1.) {
                            dscal_(&ki, &scale, &work[iv * *n + 1], &c__1);
                        }
                        work[j - 1 + iv * *n] = x[0];
                        work[j + iv * *n] = x[1];
                        i__2 = j - 2;
                        d__1 = -x[0];
                        daxpy_(&i__2, &d__1, &t[(j - 1) * t_dim1 + 1], &c__1, &work[iv * *n + 1],
                               &c__1);
                        i__2 = j - 2;
                        d__1 = -x[1];
                        daxpy_(&i__2, &d__1, &t[j * t_dim1 + 1], &c__1, &work[iv * *n + 1], &c__1);
                    }
                L60:;
                }
                if (!over) {
                    dcopy_(&ki, &work[iv * *n + 1], &c__1, &vr[is * vr_dim1 + 1], &c__1);
                    ii = idamax_(&ki, &vr[is * vr_dim1 + 1], &c__1);
                    remax = 1. / (d__1 = vr[ii + is * vr_dim1], abs(d__1));
                    dscal_(&ki, &remax, &vr[is * vr_dim1 + 1], &c__1);
                    i__2 = *n;
                    for (k = ki + 1; k <= i__2; ++k) {
                        vr[k + is * vr_dim1] = 0.;
                    }
                } else if (nb == 1) {
                    if (ki > 1) {
                        i__2 = ki - 1;
                        dgemv_((char *)"N", n, &i__2, &c_b29, &vr[vr_offset], ldvr, &work[iv * *n + 1],
                               &c__1, &work[ki + iv * *n], &vr[ki * vr_dim1 + 1], &c__1, (ftnlen)1);
                    }
                    ii = idamax_(n, &vr[ki * vr_dim1 + 1], &c__1);
                    remax = 1. / (d__1 = vr[ii + ki * vr_dim1], abs(d__1));
                    dscal_(n, &remax, &vr[ki * vr_dim1 + 1], &c__1);
                } else {
                    i__2 = *n;
                    for (k = ki + 1; k <= i__2; ++k) {
                        work[k + iv * *n] = 0.;
                    }
                    iscomplex[iv - 1] = ip;
                }
            } else {
                if ((d__1 = t[ki - 1 + ki * t_dim1], abs(d__1)) >=
                    (d__2 = t[ki + (ki - 1) * t_dim1], abs(d__2))) {
                    work[ki - 1 + (iv - 1) * *n] = 1.;
                    work[ki + iv * *n] = wi / t[ki - 1 + ki * t_dim1];
                } else {
                    work[ki - 1 + (iv - 1) * *n] = -wi / t[ki + (ki - 1) * t_dim1];
                    work[ki + iv * *n] = 1.;
                }
                work[ki + (iv - 1) * *n] = 0.;
                work[ki - 1 + iv * *n] = 0.;
                i__2 = ki - 2;
                for (k = 1; k <= i__2; ++k) {
                    work[k + (iv - 1) * *n] =
                        -work[ki - 1 + (iv - 1) * *n] * t[k + (ki - 1) * t_dim1];
                    work[k + iv * *n] = -work[ki + iv * *n] * t[k + ki * t_dim1];
                }
                jnxt = ki - 2;
                for (j = ki - 2; j >= 1; --j) {
                    if (j > jnxt) {
                        goto L90;
                    }
                    j1 = j;
                    j2 = j;
                    jnxt = j - 1;
                    if (j > 1) {
                        if (t[j + (j - 1) * t_dim1] != 0.) {
                            j1 = j - 1;
                            jnxt = j - 2;
                        }
                    }
                    if (j1 == j2) {
                        dlaln2_(&c_false, &c__1, &c__2, &smin, &c_b29, &t[j + j * t_dim1], ldt,
                                &c_b29, &c_b29, &work[j + (iv - 1) * *n], n, &wr, &wi, x, &c__2,
                                &scale, &xnorm, &ierr);
                        if (xnorm > 1.) {
                            if (work[j] > bignum / xnorm) {
                                x[0] /= xnorm;
                                x[2] /= xnorm;
                                scale /= xnorm;
                            }
                        }
                        if (scale != 1.) {
                            dscal_(&ki, &scale, &work[(iv - 1) * *n + 1], &c__1);
                            dscal_(&ki, &scale, &work[iv * *n + 1], &c__1);
                        }
                        work[j + (iv - 1) * *n] = x[0];
                        work[j + iv * *n] = x[2];
                        i__2 = j - 1;
                        d__1 = -x[0];
                        daxpy_(&i__2, &d__1, &t[j * t_dim1 + 1], &c__1, &work[(iv - 1) * *n + 1],
                               &c__1);
                        i__2 = j - 1;
                        d__1 = -x[2];
                        daxpy_(&i__2, &d__1, &t[j * t_dim1 + 1], &c__1, &work[iv * *n + 1], &c__1);
                    } else {
                        dlaln2_(&c_false, &c__2, &c__2, &smin, &c_b29, &t[j - 1 + (j - 1) * t_dim1],
                                ldt, &c_b29, &c_b29, &work[j - 1 + (iv - 1) * *n], n, &wr, &wi, x,
                                &c__2, &scale, &xnorm, &ierr);
                        if (xnorm > 1.) {
                            d__1 = work[j - 1], d__2 = work[j];
                            beta = max(d__1, d__2);
                            if (beta > bignum / xnorm) {
                                rec = 1. / xnorm;
                                x[0] *= rec;
                                x[2] *= rec;
                                x[1] *= rec;
                                x[3] *= rec;
                                scale *= rec;
                            }
                        }
                        if (scale != 1.) {
                            dscal_(&ki, &scale, &work[(iv - 1) * *n + 1], &c__1);
                            dscal_(&ki, &scale, &work[iv * *n + 1], &c__1);
                        }
                        work[j - 1 + (iv - 1) * *n] = x[0];
                        work[j + (iv - 1) * *n] = x[1];
                        work[j - 1 + iv * *n] = x[2];
                        work[j + iv * *n] = x[3];
                        i__2 = j - 2;
                        d__1 = -x[0];
                        daxpy_(&i__2, &d__1, &t[(j - 1) * t_dim1 + 1], &c__1,
                               &work[(iv - 1) * *n + 1], &c__1);
                        i__2 = j - 2;
                        d__1 = -x[1];
                        daxpy_(&i__2, &d__1, &t[j * t_dim1 + 1], &c__1, &work[(iv - 1) * *n + 1],
                               &c__1);
                        i__2 = j - 2;
                        d__1 = -x[2];
                        daxpy_(&i__2, &d__1, &t[(j - 1) * t_dim1 + 1], &c__1, &work[iv * *n + 1],
                               &c__1);
                        i__2 = j - 2;
                        d__1 = -x[3];
                        daxpy_(&i__2, &d__1, &t[j * t_dim1 + 1], &c__1, &work[iv * *n + 1], &c__1);
                    }
                L90:;
                }
                if (!over) {
                    dcopy_(&ki, &work[(iv - 1) * *n + 1], &c__1, &vr[(is - 1) * vr_dim1 + 1],
                           &c__1);
                    dcopy_(&ki, &work[iv * *n + 1], &c__1, &vr[is * vr_dim1 + 1], &c__1);
                    emax = 0.;
                    i__2 = ki;
                    for (k = 1; k <= i__2; ++k) {
                        d__3 = emax, d__4 = (d__1 = vr[k + (is - 1) * vr_dim1], abs(d__1)) +
                                            (d__2 = vr[k + is * vr_dim1], abs(d__2));
                        emax = max(d__3, d__4);
                    }
                    remax = 1. / emax;
                    dscal_(&ki, &remax, &vr[(is - 1) * vr_dim1 + 1], &c__1);
                    dscal_(&ki, &remax, &vr[is * vr_dim1 + 1], &c__1);
                    i__2 = *n;
                    for (k = ki + 1; k <= i__2; ++k) {
                        vr[k + (is - 1) * vr_dim1] = 0.;
                        vr[k + is * vr_dim1] = 0.;
                    }
                } else if (nb == 1) {
                    if (ki > 2) {
                        i__2 = ki - 2;
                        dgemv_((char *)"N", n, &i__2, &c_b29, &vr[vr_offset], ldvr,
                               &work[(iv - 1) * *n + 1], &c__1, &work[ki - 1 + (iv - 1) * *n],
                               &vr[(ki - 1) * vr_dim1 + 1], &c__1, (ftnlen)1);
                        i__2 = ki - 2;
                        dgemv_((char *)"N", n, &i__2, &c_b29, &vr[vr_offset], ldvr, &work[iv * *n + 1],
                               &c__1, &work[ki + iv * *n], &vr[ki * vr_dim1 + 1], &c__1, (ftnlen)1);
                    } else {
                        dscal_(n, &work[ki - 1 + (iv - 1) * *n], &vr[(ki - 1) * vr_dim1 + 1],
                               &c__1);
                        dscal_(n, &work[ki + iv * *n], &vr[ki * vr_dim1 + 1], &c__1);
                    }
                    emax = 0.;
                    i__2 = *n;
                    for (k = 1; k <= i__2; ++k) {
                        d__3 = emax, d__4 = (d__1 = vr[k + (ki - 1) * vr_dim1], abs(d__1)) +
                                            (d__2 = vr[k + ki * vr_dim1], abs(d__2));
                        emax = max(d__3, d__4);
                    }
                    remax = 1. / emax;
                    dscal_(n, &remax, &vr[(ki - 1) * vr_dim1 + 1], &c__1);
                    dscal_(n, &remax, &vr[ki * vr_dim1 + 1], &c__1);
                } else {
                    i__2 = *n;
                    for (k = ki + 1; k <= i__2; ++k) {
                        work[k + (iv - 1) * *n] = 0.;
                        work[k + iv * *n] = 0.;
                    }
                    iscomplex[iv - 2] = -ip;
                    iscomplex[iv - 1] = ip;
                    --iv;
                }
            }
            if (nb > 1) {
                if (ip == 0) {
                    ki2 = ki;
                } else {
                    ki2 = ki - 1;
                }
                if (iv <= 2 || ki2 == 1) {
                    i__2 = nb - iv + 1;
                    i__3 = ki2 + nb - iv;
                    dgemm_((char *)"N", (char *)"N", n, &i__2, &i__3, &c_b29, &vr[vr_offset], ldvr,
                           &work[iv * *n + 1], n, &c_b17, &work[(nb + iv) * *n + 1], n, (ftnlen)1,
                           (ftnlen)1);
                    i__2 = nb;
                    for (k = iv; k <= i__2; ++k) {
                        if (iscomplex[k - 1] == 0) {
                            ii = idamax_(n, &work[(nb + k) * *n + 1], &c__1);
                            remax = 1. / (d__1 = work[ii + (nb + k) * *n], abs(d__1));
                        } else if (iscomplex[k - 1] == 1) {
                            emax = 0.;
                            i__3 = *n;
                            for (ii = 1; ii <= i__3; ++ii) {
                                d__3 = emax,
                                d__4 = (d__1 = work[ii + (nb + k) * *n], abs(d__1)) +
                                       (d__2 = work[ii + (nb + k + 1) * *n], abs(d__2));
                                emax = max(d__3, d__4);
                            }
                            remax = 1. / emax;
                        }
                        dscal_(n, &remax, &work[(nb + k) * *n + 1], &c__1);
                    }
                    i__2 = nb - iv + 1;
                    dlacpy_((char *)"F", n, &i__2, &work[(nb + iv) * *n + 1], n, &vr[ki2 * vr_dim1 + 1],
                            ldvr, (ftnlen)1);
                    iv = nb;
                } else {
                    --iv;
                }
            }
            --is;
            if (ip != 0) {
                --is;
            }
        L140:;
        }
    }
    if (leftv) {
        iv = 1;
        ip = 0;
        is = 1;
        i__2 = *n;
        for (ki = 1; ki <= i__2; ++ki) {
            if (ip == 1) {
                ip = -1;
                goto L260;
            } else if (ki == *n) {
                ip = 0;
            } else if (t[ki + 1 + ki * t_dim1] == 0.) {
                ip = 0;
            } else {
                ip = 1;
            }
            if (somev) {
                if (!select[ki]) {
                    goto L260;
                }
            }
            wr = t[ki + ki * t_dim1];
            wi = 0.;
            if (ip != 0) {
                wi = sqrt((d__1 = t[ki + (ki + 1) * t_dim1], abs(d__1))) *
                     sqrt((d__2 = t[ki + 1 + ki * t_dim1], abs(d__2)));
            }
            d__1 = ulp * (abs(wr) + abs(wi));
            smin = max(d__1, smlnum);
            if (ip == 0) {
                work[ki + iv * *n] = 1.;
                i__3 = *n;
                for (k = ki + 1; k <= i__3; ++k) {
                    work[k + iv * *n] = -t[ki + k * t_dim1];
                }
                vmax = 1.;
                vcrit = bignum;
                jnxt = ki + 1;
                i__3 = *n;
                for (j = ki + 1; j <= i__3; ++j) {
                    if (j < jnxt) {
                        goto L170;
                    }
                    j1 = j;
                    j2 = j;
                    jnxt = j + 1;
                    if (j < *n) {
                        if (t[j + 1 + j * t_dim1] != 0.) {
                            j2 = j + 1;
                            jnxt = j + 2;
                        }
                    }
                    if (j1 == j2) {
                        if (work[j] > vcrit) {
                            rec = 1. / vmax;
                            i__4 = *n - ki + 1;
                            dscal_(&i__4, &rec, &work[ki + iv * *n], &c__1);
                            vmax = 1.;
                            vcrit = bignum;
                        }
                        i__4 = j - ki - 1;
                        work[j + iv * *n] -= ddot_(&i__4, &t[ki + 1 + j * t_dim1], &c__1,
                                                   &work[ki + 1 + iv * *n], &c__1);
                        dlaln2_(&c_false, &c__1, &c__1, &smin, &c_b29, &t[j + j * t_dim1], ldt,
                                &c_b29, &c_b29, &work[j + iv * *n], n, &wr, &c_b17, x, &c__2,
                                &scale, &xnorm, &ierr);
                        if (scale != 1.) {
                            i__4 = *n - ki + 1;
                            dscal_(&i__4, &scale, &work[ki + iv * *n], &c__1);
                        }
                        work[j + iv * *n] = x[0];
                        d__2 = (d__1 = work[j + iv * *n], abs(d__1));
                        vmax = max(d__2, vmax);
                        vcrit = bignum / vmax;
                    } else {
                        d__1 = work[j], d__2 = work[j + 1];
                        beta = max(d__1, d__2);
                        if (beta > vcrit) {
                            rec = 1. / vmax;
                            i__4 = *n - ki + 1;
                            dscal_(&i__4, &rec, &work[ki + iv * *n], &c__1);
                            vmax = 1.;
                            vcrit = bignum;
                        }
                        i__4 = j - ki - 1;
                        work[j + iv * *n] -= ddot_(&i__4, &t[ki + 1 + j * t_dim1], &c__1,
                                                   &work[ki + 1 + iv * *n], &c__1);
                        i__4 = j - ki - 1;
                        work[j + 1 + iv * *n] -= ddot_(&i__4, &t[ki + 1 + (j + 1) * t_dim1], &c__1,
                                                       &work[ki + 1 + iv * *n], &c__1);
                        dlaln2_(&c_true, &c__2, &c__1, &smin, &c_b29, &t[j + j * t_dim1], ldt,
                                &c_b29, &c_b29, &work[j + iv * *n], n, &wr, &c_b17, x, &c__2,
                                &scale, &xnorm, &ierr);
                        if (scale != 1.) {
                            i__4 = *n - ki + 1;
                            dscal_(&i__4, &scale, &work[ki + iv * *n], &c__1);
                        }
                        work[j + iv * *n] = x[0];
                        work[j + 1 + iv * *n] = x[1];
                        d__3 = (d__1 = work[j + iv * *n], abs(d__1)),
                        d__4 = (d__2 = work[j + 1 + iv * *n], abs(d__2)), d__3 = max(d__3, d__4);
                        vmax = max(d__3, vmax);
                        vcrit = bignum / vmax;
                    }
                L170:;
                }
                if (!over) {
                    i__3 = *n - ki + 1;
                    dcopy_(&i__3, &work[ki + iv * *n], &c__1, &vl[ki + is * vl_dim1], &c__1);
                    i__3 = *n - ki + 1;
                    ii = idamax_(&i__3, &vl[ki + is * vl_dim1], &c__1) + ki - 1;
                    remax = 1. / (d__1 = vl[ii + is * vl_dim1], abs(d__1));
                    i__3 = *n - ki + 1;
                    dscal_(&i__3, &remax, &vl[ki + is * vl_dim1], &c__1);
                    i__3 = ki - 1;
                    for (k = 1; k <= i__3; ++k) {
                        vl[k + is * vl_dim1] = 0.;
                    }
                } else if (nb == 1) {
                    if (ki < *n) {
                        i__3 = *n - ki;
                        dgemv_((char *)"N", n, &i__3, &c_b29, &vl[(ki + 1) * vl_dim1 + 1], ldvl,
                               &work[ki + 1 + iv * *n], &c__1, &work[ki + iv * *n],
                               &vl[ki * vl_dim1 + 1], &c__1, (ftnlen)1);
                    }
                    ii = idamax_(n, &vl[ki * vl_dim1 + 1], &c__1);
                    remax = 1. / (d__1 = vl[ii + ki * vl_dim1], abs(d__1));
                    dscal_(n, &remax, &vl[ki * vl_dim1 + 1], &c__1);
                } else {
                    i__3 = ki - 1;
                    for (k = 1; k <= i__3; ++k) {
                        work[k + iv * *n] = 0.;
                    }
                    iscomplex[iv - 1] = ip;
                }
            } else {
                if ((d__1 = t[ki + (ki + 1) * t_dim1], abs(d__1)) >=
                    (d__2 = t[ki + 1 + ki * t_dim1], abs(d__2))) {
                    work[ki + iv * *n] = wi / t[ki + (ki + 1) * t_dim1];
                    work[ki + 1 + (iv + 1) * *n] = 1.;
                } else {
                    work[ki + iv * *n] = 1.;
                    work[ki + 1 + (iv + 1) * *n] = -wi / t[ki + 1 + ki * t_dim1];
                }
                work[ki + 1 + iv * *n] = 0.;
                work[ki + (iv + 1) * *n] = 0.;
                i__3 = *n;
                for (k = ki + 2; k <= i__3; ++k) {
                    work[k + iv * *n] = -work[ki + iv * *n] * t[ki + k * t_dim1];
                    work[k + (iv + 1) * *n] =
                        -work[ki + 1 + (iv + 1) * *n] * t[ki + 1 + k * t_dim1];
                }
                vmax = 1.;
                vcrit = bignum;
                jnxt = ki + 2;
                i__3 = *n;
                for (j = ki + 2; j <= i__3; ++j) {
                    if (j < jnxt) {
                        goto L200;
                    }
                    j1 = j;
                    j2 = j;
                    jnxt = j + 1;
                    if (j < *n) {
                        if (t[j + 1 + j * t_dim1] != 0.) {
                            j2 = j + 1;
                            jnxt = j + 2;
                        }
                    }
                    if (j1 == j2) {
                        if (work[j] > vcrit) {
                            rec = 1. / vmax;
                            i__4 = *n - ki + 1;
                            dscal_(&i__4, &rec, &work[ki + iv * *n], &c__1);
                            i__4 = *n - ki + 1;
                            dscal_(&i__4, &rec, &work[ki + (iv + 1) * *n], &c__1);
                            vmax = 1.;
                            vcrit = bignum;
                        }
                        i__4 = j - ki - 2;
                        work[j + iv * *n] -= ddot_(&i__4, &t[ki + 2 + j * t_dim1], &c__1,
                                                   &work[ki + 2 + iv * *n], &c__1);
                        i__4 = j - ki - 2;
                        work[j + (iv + 1) * *n] -= ddot_(&i__4, &t[ki + 2 + j * t_dim1], &c__1,
                                                         &work[ki + 2 + (iv + 1) * *n], &c__1);
                        d__1 = -wi;
                        dlaln2_(&c_false, &c__1, &c__2, &smin, &c_b29, &t[j + j * t_dim1], ldt,
                                &c_b29, &c_b29, &work[j + iv * *n], n, &wr, &d__1, x, &c__2, &scale,
                                &xnorm, &ierr);
                        if (scale != 1.) {
                            i__4 = *n - ki + 1;
                            dscal_(&i__4, &scale, &work[ki + iv * *n], &c__1);
                            i__4 = *n - ki + 1;
                            dscal_(&i__4, &scale, &work[ki + (iv + 1) * *n], &c__1);
                        }
                        work[j + iv * *n] = x[0];
                        work[j + (iv + 1) * *n] = x[2];
                        d__3 = (d__1 = work[j + iv * *n], abs(d__1)),
                        d__4 = (d__2 = work[j + (iv + 1) * *n], abs(d__2)), d__3 = max(d__3, d__4);
                        vmax = max(d__3, vmax);
                        vcrit = bignum / vmax;
                    } else {
                        d__1 = work[j], d__2 = work[j + 1];
                        beta = max(d__1, d__2);
                        if (beta > vcrit) {
                            rec = 1. / vmax;
                            i__4 = *n - ki + 1;
                            dscal_(&i__4, &rec, &work[ki + iv * *n], &c__1);
                            i__4 = *n - ki + 1;
                            dscal_(&i__4, &rec, &work[ki + (iv + 1) * *n], &c__1);
                            vmax = 1.;
                            vcrit = bignum;
                        }
                        i__4 = j - ki - 2;
                        work[j + iv * *n] -= ddot_(&i__4, &t[ki + 2 + j * t_dim1], &c__1,
                                                   &work[ki + 2 + iv * *n], &c__1);
                        i__4 = j - ki - 2;
                        work[j + (iv + 1) * *n] -= ddot_(&i__4, &t[ki + 2 + j * t_dim1], &c__1,
                                                         &work[ki + 2 + (iv + 1) * *n], &c__1);
                        i__4 = j - ki - 2;
                        work[j + 1 + iv * *n] -= ddot_(&i__4, &t[ki + 2 + (j + 1) * t_dim1], &c__1,
                                                       &work[ki + 2 + iv * *n], &c__1);
                        i__4 = j - ki - 2;
                        work[j + 1 + (iv + 1) * *n] -=
                            ddot_(&i__4, &t[ki + 2 + (j + 1) * t_dim1], &c__1,
                                  &work[ki + 2 + (iv + 1) * *n], &c__1);
                        d__1 = -wi;
                        dlaln2_(&c_true, &c__2, &c__2, &smin, &c_b29, &t[j + j * t_dim1], ldt,
                                &c_b29, &c_b29, &work[j + iv * *n], n, &wr, &d__1, x, &c__2, &scale,
                                &xnorm, &ierr);
                        if (scale != 1.) {
                            i__4 = *n - ki + 1;
                            dscal_(&i__4, &scale, &work[ki + iv * *n], &c__1);
                            i__4 = *n - ki + 1;
                            dscal_(&i__4, &scale, &work[ki + (iv + 1) * *n], &c__1);
                        }
                        work[j + iv * *n] = x[0];
                        work[j + (iv + 1) * *n] = x[2];
                        work[j + 1 + iv * *n] = x[1];
                        work[j + 1 + (iv + 1) * *n] = x[3];
                        d__1 = abs(x[0]), d__2 = abs(x[2]), d__1 = max(d__1, d__2),
                        d__2 = abs(x[1]), d__1 = max(d__1, d__2), d__2 = abs(x[3]),
                        d__1 = max(d__1, d__2);
                        vmax = max(d__1, vmax);
                        vcrit = bignum / vmax;
                    }
                L200:;
                }
                if (!over) {
                    i__3 = *n - ki + 1;
                    dcopy_(&i__3, &work[ki + iv * *n], &c__1, &vl[ki + is * vl_dim1], &c__1);
                    i__3 = *n - ki + 1;
                    dcopy_(&i__3, &work[ki + (iv + 1) * *n], &c__1, &vl[ki + (is + 1) * vl_dim1],
                           &c__1);
                    emax = 0.;
                    i__3 = *n;
                    for (k = ki; k <= i__3; ++k) {
                        d__3 = emax, d__4 = (d__1 = vl[k + is * vl_dim1], abs(d__1)) +
                                            (d__2 = vl[k + (is + 1) * vl_dim1], abs(d__2));
                        emax = max(d__3, d__4);
                    }
                    remax = 1. / emax;
                    i__3 = *n - ki + 1;
                    dscal_(&i__3, &remax, &vl[ki + is * vl_dim1], &c__1);
                    i__3 = *n - ki + 1;
                    dscal_(&i__3, &remax, &vl[ki + (is + 1) * vl_dim1], &c__1);
                    i__3 = ki - 1;
                    for (k = 1; k <= i__3; ++k) {
                        vl[k + is * vl_dim1] = 0.;
                        vl[k + (is + 1) * vl_dim1] = 0.;
                    }
                } else if (nb == 1) {
                    if (ki < *n - 1) {
                        i__3 = *n - ki - 1;
                        dgemv_((char *)"N", n, &i__3, &c_b29, &vl[(ki + 2) * vl_dim1 + 1], ldvl,
                               &work[ki + 2 + iv * *n], &c__1, &work[ki + iv * *n],
                               &vl[ki * vl_dim1 + 1], &c__1, (ftnlen)1);
                        i__3 = *n - ki - 1;
                        dgemv_((char *)"N", n, &i__3, &c_b29, &vl[(ki + 2) * vl_dim1 + 1], ldvl,
                               &work[ki + 2 + (iv + 1) * *n], &c__1, &work[ki + 1 + (iv + 1) * *n],
                               &vl[(ki + 1) * vl_dim1 + 1], &c__1, (ftnlen)1);
                    } else {
                        dscal_(n, &work[ki + iv * *n], &vl[ki * vl_dim1 + 1], &c__1);
                        dscal_(n, &work[ki + 1 + (iv + 1) * *n], &vl[(ki + 1) * vl_dim1 + 1],
                               &c__1);
                    }
                    emax = 0.;
                    i__3 = *n;
                    for (k = 1; k <= i__3; ++k) {
                        d__3 = emax, d__4 = (d__1 = vl[k + ki * vl_dim1], abs(d__1)) +
                                            (d__2 = vl[k + (ki + 1) * vl_dim1], abs(d__2));
                        emax = max(d__3, d__4);
                    }
                    remax = 1. / emax;
                    dscal_(n, &remax, &vl[ki * vl_dim1 + 1], &c__1);
                    dscal_(n, &remax, &vl[(ki + 1) * vl_dim1 + 1], &c__1);
                } else {
                    i__3 = ki - 1;
                    for (k = 1; k <= i__3; ++k) {
                        work[k + iv * *n] = 0.;
                        work[k + (iv + 1) * *n] = 0.;
                    }
                    iscomplex[iv - 1] = ip;
                    iscomplex[iv] = -ip;
                    ++iv;
                }
            }
            if (nb > 1) {
                if (ip == 0) {
                    ki2 = ki;
                } else {
                    ki2 = ki + 1;
                }
                if (iv >= nb - 1 || ki2 == *n) {
                    i__3 = *n - ki2 + iv;
                    dgemm_((char *)"N", (char *)"N", n, &iv, &i__3, &c_b29, &vl[(ki2 - iv + 1) * vl_dim1 + 1], ldvl,
                           &work[ki2 - iv + 1 + *n], n, &c_b17, &work[(nb + 1) * *n + 1], n,
                           (ftnlen)1, (ftnlen)1);
                    i__3 = iv;
                    for (k = 1; k <= i__3; ++k) {
                        if (iscomplex[k - 1] == 0) {
                            ii = idamax_(n, &work[(nb + k) * *n + 1], &c__1);
                            remax = 1. / (d__1 = work[ii + (nb + k) * *n], abs(d__1));
                        } else if (iscomplex[k - 1] == 1) {
                            emax = 0.;
                            i__4 = *n;
                            for (ii = 1; ii <= i__4; ++ii) {
                                d__3 = emax,
                                d__4 = (d__1 = work[ii + (nb + k) * *n], abs(d__1)) +
                                       (d__2 = work[ii + (nb + k + 1) * *n], abs(d__2));
                                emax = max(d__3, d__4);
                            }
                            remax = 1. / emax;
                        }
                        dscal_(n, &remax, &work[(nb + k) * *n + 1], &c__1);
                    }
                    dlacpy_((char *)"F", n, &iv, &work[(nb + 1) * *n + 1], n,
                            &vl[(ki2 - iv + 1) * vl_dim1 + 1], ldvl, (ftnlen)1);
                    iv = 1;
                } else {
                    ++iv;
                }
            }
            ++is;
            if (ip != 0) {
                ++is;
            }
        L260:;
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
