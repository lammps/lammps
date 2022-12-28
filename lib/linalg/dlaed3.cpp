#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static doublereal c_b22 = 1.;
static doublereal c_b23 = 0.;
int dlaed3_(integer *k, integer *n, integer *n1, doublereal *d__, doublereal *q, integer *ldq,
            doublereal *rho, doublereal *dlamda, doublereal *q2, integer *indx, integer *ctot,
            doublereal *w, doublereal *s, integer *info)
{
    integer q_dim1, q_offset, i__1, i__2;
    doublereal d__1;
    double sqrt(doublereal), d_lmp_sign(doublereal *, doublereal *);
    integer i__, j, n2, n12, ii, n23, iq2;
    doublereal temp;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
                      integer *, doublereal *, integer *, doublereal *, doublereal *, integer *,
                      ftnlen, ftnlen),
        dcopy_(integer *, doublereal *, integer *, doublereal *, integer *),
        dlaed4_(integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *,
                doublereal *, integer *);
    extern doublereal dlamc3_(doublereal *, doublereal *);
    extern int dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *,
                       integer *, ftnlen),
        dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                ftnlen),
        xerbla_(char *, integer *, ftnlen);
    --d__;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --dlamda;
    --q2;
    --indx;
    --ctot;
    --w;
    --s;
    *info = 0;
    if (*k < 0) {
        *info = -1;
    } else if (*n < *k) {
        *info = -2;
    } else if (*ldq < max(1, *n)) {
        *info = -6;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLAED3", &i__1, (ftnlen)6);
        return 0;
    }
    if (*k == 0) {
        return 0;
    }
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
        dlamda[i__] = dlamc3_(&dlamda[i__], &dlamda[i__]) - dlamda[i__];
    }
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
        dlaed4_(k, &j, &dlamda[1], &w[1], &q[j * q_dim1 + 1], rho, &d__[j], info);
        if (*info != 0) {
            goto L120;
        }
    }
    if (*k == 1) {
        goto L110;
    }
    if (*k == 2) {
        i__1 = *k;
        for (j = 1; j <= i__1; ++j) {
            w[1] = q[j * q_dim1 + 1];
            w[2] = q[j * q_dim1 + 2];
            ii = indx[1];
            q[j * q_dim1 + 1] = w[ii];
            ii = indx[2];
            q[j * q_dim1 + 2] = w[ii];
        }
        goto L110;
    }
    dcopy_(k, &w[1], &c__1, &s[1], &c__1);
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
        w[i__] = d_lmp_sign(&d__1, &s[i__]);
    }
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
        i__2 = *k;
        for (i__ = 1; i__ <= i__2; ++i__) {
            s[i__] = w[i__] / q[i__ + j * q_dim1];
        }
        temp = dnrm2_(k, &s[1], &c__1);
        i__2 = *k;
        for (i__ = 1; i__ <= i__2; ++i__) {
            ii = indx[i__];
            q[i__ + j * q_dim1] = s[ii] / temp;
        }
    }
L110:
    n2 = *n - *n1;
    n12 = ctot[1] + ctot[2];
    n23 = ctot[2] + ctot[3];
    dlacpy_((char *)"A", &n23, k, &q[ctot[1] + 1 + q_dim1], ldq, &s[1], &n23, (ftnlen)1);
    iq2 = *n1 * n12 + 1;
    if (n23 != 0) {
        dgemm_((char *)"N", (char *)"N", &n2, k, &n23, &c_b22, &q2[iq2], &n2, &s[1], &n23, &c_b23,
               &q[*n1 + 1 + q_dim1], ldq, (ftnlen)1, (ftnlen)1);
    } else {
        dlaset_((char *)"A", &n2, k, &c_b23, &c_b23, &q[*n1 + 1 + q_dim1], ldq, (ftnlen)1);
    }
    dlacpy_((char *)"A", &n12, k, &q[q_offset], ldq, &s[1], &n12, (ftnlen)1);
    if (n12 != 0) {
        dgemm_((char *)"N", (char *)"N", n1, k, &n12, &c_b22, &q2[1], n1, &s[1], &n12, &c_b23, &q[q_offset], ldq,
               (ftnlen)1, (ftnlen)1);
    } else {
        dlaset_((char *)"A", n1, k, &c_b23, &c_b23, &q[q_dim1 + 1], ldq, (ftnlen)1);
    }
L120:
    return 0;
}
#ifdef __cplusplus
}
#endif
