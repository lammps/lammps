#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
int dlaed9_(integer *k, integer *kstart, integer *kstop, integer *n, doublereal *d__, doublereal *q,
            integer *ldq, doublereal *rho, doublereal *dlamda, doublereal *w, doublereal *s,
            integer *lds, integer *info)
{
    integer q_dim1, q_offset, s_dim1, s_offset, i__1, i__2;
    doublereal d__1;
    double sqrt(doublereal), d_lmp_sign(doublereal *, doublereal *);
    integer i__, j;
    doublereal temp;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *),
        dlaed4_(integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *,
                doublereal *, integer *);
    extern doublereal dlamc3_(doublereal *, doublereal *);
    extern int xerbla_(char *, integer *, ftnlen);
    --d__;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --dlamda;
    --w;
    s_dim1 = *lds;
    s_offset = 1 + s_dim1;
    s -= s_offset;
    *info = 0;
    if (*k < 0) {
        *info = -1;
    } else if (*kstart < 1 || *kstart > max(1, *k)) {
        *info = -2;
    } else if (max(1, *kstop) < *kstart || *kstop > max(1, *k)) {
        *info = -3;
    } else if (*n < *k) {
        *info = -4;
    } else if (*ldq < max(1, *k)) {
        *info = -7;
    } else if (*lds < max(1, *k)) {
        *info = -12;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLAED9", &i__1, (ftnlen)6);
        return 0;
    }
    if (*k == 0) {
        return 0;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        dlamda[i__] = dlamc3_(&dlamda[i__], &dlamda[i__]) - dlamda[i__];
    }
    i__1 = *kstop;
    for (j = *kstart; j <= i__1; ++j) {
        dlaed4_(k, &j, &dlamda[1], &w[1], &q[j * q_dim1 + 1], rho, &d__[j], info);
        if (*info != 0) {
            goto L120;
        }
    }
    if (*k == 1 || *k == 2) {
        i__1 = *k;
        for (i__ = 1; i__ <= i__1; ++i__) {
            i__2 = *k;
            for (j = 1; j <= i__2; ++j) {
                s[j + i__ * s_dim1] = q[j + i__ * q_dim1];
            }
        }
        goto L120;
    }
    dcopy_(k, &w[1], &c__1, &s[s_offset], &c__1);
    i__1 = *ldq + 1;
    dcopy_(k, &q[q_offset], &i__1, &w[1], &c__1);
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
        i__2 = j - 1;
        for (i__ = 1; i__ <= i__2; ++i__) {
            w[i__] *= q[i__ + j * q_dim1] / (dlamda[i__] - dlamda[j]);
        }
        i__2 = *k;
        for (i__ = j + 1; i__ <= i__2; ++i__) {
            w[i__] *= q[i__ + j * q_dim1] / (dlamda[i__] - dlamda[j]);
        }
    }
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
        d__1 = sqrt(-w[i__]);
        w[i__] = d_lmp_sign(&d__1, &s[i__ + s_dim1]);
    }
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
        i__2 = *k;
        for (i__ = 1; i__ <= i__2; ++i__) {
            q[i__ + j * q_dim1] = w[i__] / q[i__ + j * q_dim1];
        }
        temp = dnrm2_(k, &q[j * q_dim1 + 1], &c__1);
        i__2 = *k;
        for (i__ = 1; i__ <= i__2; ++i__) {
            s[i__ + j * s_dim1] = q[i__ + j * q_dim1] / temp;
        }
    }
L120:
    return 0;
}
#ifdef __cplusplus
}
#endif
