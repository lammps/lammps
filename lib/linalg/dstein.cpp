#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__2 = 2;
static integer c__1 = 1;
static integer c_n1 = -1;
int dstein_(integer *n, doublereal *d__, doublereal *e, integer *m, doublereal *w, integer *iblock,
            integer *isplit, doublereal *z__, integer *ldz, doublereal *work, integer *iwork,
            integer *ifail, integer *info)
{
    integer z_dim1, z_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3, d__4, d__5;
    double sqrt(doublereal);
    integer i__, j, b1, j1, bn;
    doublereal xj, scl, eps, sep, nrm, tol;
    integer its;
    doublereal xjm, ztr, eps1;
    integer jblk, nblk;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer jmax;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern int dscal_(integer *, doublereal *, doublereal *, integer *);
    integer iseed[4], gpind, iinfo;
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *),
        daxpy_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *);
    doublereal ortol;
    integer indrv1, indrv2, indrv3, indrv4, indrv5;
    extern doublereal dlamch_(char *, ftnlen);
    extern int dlagtf_(integer *, doublereal *, doublereal *, doublereal *, doublereal *,
                       doublereal *, doublereal *, integer *, integer *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern int xerbla_(char *, integer *, ftnlen),
        dlagts_(integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *,
                integer *, doublereal *, doublereal *, integer *);
    integer nrmchk;
    extern int dlarnv_(integer *, integer *, integer *, doublereal *);
    integer blksiz;
    doublereal onenrm, dtpcrt, pertol;
    --d__;
    --e;
    --w;
    --iblock;
    --isplit;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --iwork;
    --ifail;
    *info = 0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
        ifail[i__] = 0;
    }
    if (*n < 0) {
        *info = -1;
    } else if (*m < 0 || *m > *n) {
        *info = -4;
    } else if (*ldz < max(1, *n)) {
        *info = -9;
    } else {
        i__1 = *m;
        for (j = 2; j <= i__1; ++j) {
            if (iblock[j] < iblock[j - 1]) {
                *info = -6;
                goto L30;
            }
            if (iblock[j] == iblock[j - 1] && w[j] < w[j - 1]) {
                *info = -5;
                goto L30;
            }
        }
    L30:;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DSTEIN", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0 || *m == 0) {
        return 0;
    } else if (*n == 1) {
        z__[z_dim1 + 1] = 1.;
        return 0;
    }
    eps = dlamch_((char *)"Precision", (ftnlen)9);
    for (i__ = 1; i__ <= 4; ++i__) {
        iseed[i__ - 1] = 1;
    }
    indrv1 = 0;
    indrv2 = indrv1 + *n;
    indrv3 = indrv2 + *n;
    indrv4 = indrv3 + *n;
    indrv5 = indrv4 + *n;
    j1 = 1;
    i__1 = iblock[*m];
    for (nblk = 1; nblk <= i__1; ++nblk) {
        if (nblk == 1) {
            b1 = 1;
        } else {
            b1 = isplit[nblk - 1] + 1;
        }
        bn = isplit[nblk];
        blksiz = bn - b1 + 1;
        if (blksiz == 1) {
            goto L60;
        }
        gpind = j1;
        onenrm = (d__1 = d__[b1], abs(d__1)) + (d__2 = e[b1], abs(d__2));
        d__3 = onenrm, d__4 = (d__1 = d__[bn], abs(d__1)) + (d__2 = e[bn - 1], abs(d__2));
        onenrm = max(d__3, d__4);
        i__2 = bn - 1;
        for (i__ = b1 + 1; i__ <= i__2; ++i__) {
            d__4 = onenrm, d__5 = (d__1 = d__[i__], abs(d__1)) + (d__2 = e[i__ - 1], abs(d__2)) +
                                  (d__3 = e[i__], abs(d__3));
            onenrm = max(d__4, d__5);
        }
        ortol = onenrm * .001;
        dtpcrt = sqrt(.1 / blksiz);
    L60:
        jblk = 0;
        i__2 = *m;
        for (j = j1; j <= i__2; ++j) {
            if (iblock[j] != nblk) {
                j1 = j;
                goto L160;
            }
            ++jblk;
            xj = w[j];
            if (blksiz == 1) {
                work[indrv1 + 1] = 1.;
                goto L120;
            }
            if (jblk > 1) {
                eps1 = (d__1 = eps * xj, abs(d__1));
                pertol = eps1 * 10.;
                sep = xj - xjm;
                if (sep < pertol) {
                    xj = xjm + pertol;
                }
            }
            its = 0;
            nrmchk = 0;
            dlarnv_(&c__2, iseed, &blksiz, &work[indrv1 + 1]);
            dcopy_(&blksiz, &d__[b1], &c__1, &work[indrv4 + 1], &c__1);
            i__3 = blksiz - 1;
            dcopy_(&i__3, &e[b1], &c__1, &work[indrv2 + 2], &c__1);
            i__3 = blksiz - 1;
            dcopy_(&i__3, &e[b1], &c__1, &work[indrv3 + 1], &c__1);
            tol = 0.;
            dlagtf_(&blksiz, &work[indrv4 + 1], &xj, &work[indrv2 + 2], &work[indrv3 + 1], &tol,
                    &work[indrv5 + 1], &iwork[1], &iinfo);
        L70:
            ++its;
            if (its > 5) {
                goto L100;
            }
            jmax = idamax_(&blksiz, &work[indrv1 + 1], &c__1);
            d__3 = eps, d__4 = (d__1 = work[indrv4 + blksiz], abs(d__1));
            scl = blksiz * onenrm * max(d__3, d__4) / (d__2 = work[indrv1 + jmax], abs(d__2));
            dscal_(&blksiz, &scl, &work[indrv1 + 1], &c__1);
            dlagts_(&c_n1, &blksiz, &work[indrv4 + 1], &work[indrv2 + 2], &work[indrv3 + 1],
                    &work[indrv5 + 1], &iwork[1], &work[indrv1 + 1], &tol, &iinfo);
            if (jblk == 1) {
                goto L90;
            }
            if ((d__1 = xj - xjm, abs(d__1)) > ortol) {
                gpind = j;
            }
            if (gpind != j) {
                i__3 = j - 1;
                for (i__ = gpind; i__ <= i__3; ++i__) {
                    ztr = -ddot_(&blksiz, &work[indrv1 + 1], &c__1, &z__[b1 + i__ * z_dim1], &c__1);
                    daxpy_(&blksiz, &ztr, &z__[b1 + i__ * z_dim1], &c__1, &work[indrv1 + 1], &c__1);
                }
            }
        L90:
            jmax = idamax_(&blksiz, &work[indrv1 + 1], &c__1);
            nrm = (d__1 = work[indrv1 + jmax], abs(d__1));
            if (nrm < dtpcrt) {
                goto L70;
            }
            ++nrmchk;
            if (nrmchk < 3) {
                goto L70;
            }
            goto L110;
        L100:
            ++(*info);
            ifail[*info] = j;
        L110:
            scl = 1. / dnrm2_(&blksiz, &work[indrv1 + 1], &c__1);
            jmax = idamax_(&blksiz, &work[indrv1 + 1], &c__1);
            if (work[indrv1 + jmax] < 0.) {
                scl = -scl;
            }
            dscal_(&blksiz, &scl, &work[indrv1 + 1], &c__1);
        L120:
            i__3 = *n;
            for (i__ = 1; i__ <= i__3; ++i__) {
                z__[i__ + j * z_dim1] = 0.;
            }
            i__3 = blksiz;
            for (i__ = 1; i__ <= i__3; ++i__) {
                z__[b1 + i__ - 1 + j * z_dim1] = work[indrv1 + i__];
            }
            xjm = xj;
        }
    L160:;
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
