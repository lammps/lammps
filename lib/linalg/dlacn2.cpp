#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
int dlacn2_(integer *n, doublereal *v, doublereal *x, integer *isgn, doublereal *est, integer *kase,
            integer *isave)
{
    integer i__1;
    doublereal d__1;
    integer i_lmp_dnnt(doublereal *);
    integer i__;
    doublereal xs, temp;
    extern doublereal dasum_(integer *, doublereal *, integer *);
    integer jlast;
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    extern integer idamax_(integer *, doublereal *, integer *);
    doublereal altsgn, estold;
    --isave;
    --isgn;
    --x;
    --v;
    if (*kase == 0) {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            x[i__] = 1. / (doublereal)(*n);
        }
        *kase = 1;
        isave[1] = 1;
        return 0;
    }
    switch (isave[1]) {
        case 1:
            goto L20;
        case 2:
            goto L40;
        case 3:
            goto L70;
        case 4:
            goto L110;
        case 5:
            goto L140;
    }
L20:
    if (*n == 1) {
        v[1] = x[1];
        *est = abs(v[1]);
        goto L150;
    }
    *est = dasum_(n, &x[1], &c__1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        if (x[i__] >= 0.) {
            x[i__] = 1.;
        } else {
            x[i__] = -1.;
        }
        isgn[i__] = i_lmp_dnnt(&x[i__]);
    }
    *kase = 2;
    isave[1] = 2;
    return 0;
L40:
    isave[2] = idamax_(n, &x[1], &c__1);
    isave[3] = 2;
L50:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        x[i__] = 0.;
    }
    x[isave[2]] = 1.;
    *kase = 1;
    isave[1] = 3;
    return 0;
L70:
    dcopy_(n, &x[1], &c__1, &v[1], &c__1);
    estold = *est;
    *est = dasum_(n, &v[1], &c__1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        if (x[i__] >= 0.) {
            xs = 1.;
        } else {
            xs = -1.;
        }
        if (i_lmp_dnnt(&xs) != isgn[i__]) {
            goto L90;
        }
    }
    goto L120;
L90:
    if (*est <= estold) {
        goto L120;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        if (x[i__] >= 0.) {
            x[i__] = 1.;
        } else {
            x[i__] = -1.;
        }
        isgn[i__] = i_lmp_dnnt(&x[i__]);
    }
    *kase = 2;
    isave[1] = 4;
    return 0;
L110:
    jlast = isave[2];
    isave[2] = idamax_(n, &x[1], &c__1);
    if (x[jlast] != (d__1 = x[isave[2]], abs(d__1)) && isave[3] < 5) {
        ++isave[3];
        goto L50;
    }
L120:
    altsgn = 1.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        x[i__] = altsgn * ((doublereal)(i__ - 1) / (doublereal)(*n - 1) + 1.);
        altsgn = -altsgn;
    }
    *kase = 1;
    isave[1] = 5;
    return 0;
L140:
    temp = dasum_(n, &x[1], &c__1) / (doublereal)(*n * 3) * 2.;
    if (temp > *est) {
        dcopy_(n, &x[1], &c__1, &v[1], &c__1);
        *est = temp;
    }
L150:
    *kase = 0;
    return 0;
}
#ifdef __cplusplus
}
#endif
