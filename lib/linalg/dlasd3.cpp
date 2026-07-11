#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b12 = 1.;
static doublereal c_b25 = 0.;
int dlasd3_(integer *nl, integer *nr, integer *sqre, integer *k, doublereal *d__, doublereal *q,
            integer *ldq, doublereal *dsigma, doublereal *u, integer *ldu, doublereal *u2,
            integer *ldu2, doublereal *vt, integer *ldvt, doublereal *vt2, integer *ldvt2,
            integer *idxc, integer *ctot, doublereal *z__, integer *info)
{
    integer q_dim1, q_offset, u_dim1, u_offset, u2_dim1, u2_offset, vt_dim1, vt_offset, vt2_dim1,
        vt2_offset, i__1, i__2;
    doublereal d__1, d__2;
    double sqrt(doublereal), d_lmp_sign(doublereal *, doublereal *);
    integer i__, j, m, n, jc;
    doublereal rho;
    integer nlp1, nlp2, nrp1;
    doublereal temp;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
                      integer *, doublereal *, integer *, doublereal *, doublereal *, integer *,
                      ftnlen, ftnlen);
    integer ctemp;
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer ktemp;
    extern int dlasd4_(integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *,
                       doublereal *, doublereal *, integer *),
        dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *,
                doublereal *, integer *, integer *, ftnlen),
        dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *,
                ftnlen),
        xerbla_(char *, integer *, ftnlen);
    --d__;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --dsigma;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    u2_dim1 = *ldu2;
    u2_offset = 1 + u2_dim1;
    u2 -= u2_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    vt2_dim1 = *ldvt2;
    vt2_offset = 1 + vt2_dim1;
    vt2 -= vt2_offset;
    --idxc;
    --ctot;
    --z__;
    *info = 0;
    if (*nl < 1) {
        *info = -1;
    } else if (*nr < 1) {
        *info = -2;
    } else if (*sqre != 1 && *sqre != 0) {
        *info = -3;
    }
    n = *nl + *nr + 1;
    m = n + *sqre;
    nlp1 = *nl + 1;
    nlp2 = *nl + 2;
    if (*k < 1 || *k > n) {
        *info = -4;
    } else if (*ldq < *k) {
        *info = -7;
    } else if (*ldu < n) {
        *info = -10;
    } else if (*ldu2 < n) {
        *info = -12;
    } else if (*ldvt < m) {
        *info = -14;
    } else if (*ldvt2 < m) {
        *info = -16;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLASD3", &i__1, (ftnlen)6);
        return 0;
    }
    if (*k == 1) {
        d__[1] = abs(z__[1]);
        dcopy_(&m, &vt2[vt2_dim1 + 1], ldvt2, &vt[vt_dim1 + 1], ldvt);
        if (z__[1] > 0.) {
            dcopy_(&n, &u2[u2_dim1 + 1], &c__1, &u[u_dim1 + 1], &c__1);
        } else {
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                u[i__ + u_dim1] = -u2[i__ + u2_dim1];
            }
        }
        return 0;
    }
    dcopy_(k, &z__[1], &c__1, &q[q_offset], &c__1);
    rho = dnrm2_(k, &z__[1], &c__1);
    dlascl_((char *)"G", &c__0, &c__0, &rho, &c_b12, k, &c__1, &z__[1], k, info, (ftnlen)1);
    rho *= rho;
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
        dlasd4_(k, &j, &dsigma[1], &z__[1], &u[j * u_dim1 + 1], &rho, &d__[j], &vt[j * vt_dim1 + 1],
                info);
        if (*info != 0) {
            return 0;
        }
    }
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
        z__[i__] = u[i__ + *k * u_dim1] * vt[i__ + *k * vt_dim1];
        i__2 = i__ - 1;
        for (j = 1; j <= i__2; ++j) {
            z__[i__] *= u[i__ + j * u_dim1] * vt[i__ + j * vt_dim1] / (dsigma[i__] - dsigma[j]) /
                        (dsigma[i__] + dsigma[j]);
        }
        i__2 = *k - 1;
        for (j = i__; j <= i__2; ++j) {
            z__[i__] *= u[i__ + j * u_dim1] * vt[i__ + j * vt_dim1] /
                        (dsigma[i__] - dsigma[j + 1]) / (dsigma[i__] + dsigma[j + 1]);
        }
        d__2 = sqrt((d__1 = z__[i__], abs(d__1)));
        z__[i__] = d_lmp_sign(&d__2, &q[i__ + q_dim1]);
    }
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
        vt[i__ * vt_dim1 + 1] = z__[1] / u[i__ * u_dim1 + 1] / vt[i__ * vt_dim1 + 1];
        u[i__ * u_dim1 + 1] = -1.;
        i__2 = *k;
        for (j = 2; j <= i__2; ++j) {
            vt[j + i__ * vt_dim1] = z__[j] / u[j + i__ * u_dim1] / vt[j + i__ * vt_dim1];
            u[j + i__ * u_dim1] = dsigma[j] * vt[j + i__ * vt_dim1];
        }
        temp = dnrm2_(k, &u[i__ * u_dim1 + 1], &c__1);
        q[i__ * q_dim1 + 1] = u[i__ * u_dim1 + 1] / temp;
        i__2 = *k;
        for (j = 2; j <= i__2; ++j) {
            jc = idxc[j];
            q[j + i__ * q_dim1] = u[jc + i__ * u_dim1] / temp;
        }
    }
    if (*k == 2) {
        dgemm_((char *)"N", (char *)"N", &n, k, k, &c_b12, &u2[u2_offset], ldu2, &q[q_offset], ldq, &c_b25,
               &u[u_offset], ldu, (ftnlen)1, (ftnlen)1);
        goto L100;
    }
    if (ctot[1] > 0) {
        dgemm_((char *)"N", (char *)"N", nl, k, &ctot[1], &c_b12, &u2[(u2_dim1 << 1) + 1], ldu2, &q[q_dim1 + 2],
               ldq, &c_b25, &u[u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1);
        if (ctot[3] > 0) {
            ktemp = ctot[1] + 2 + ctot[2];
            dgemm_((char *)"N", (char *)"N", nl, k, &ctot[3], &c_b12, &u2[ktemp * u2_dim1 + 1], ldu2,
                   &q[ktemp + q_dim1], ldq, &c_b12, &u[u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1);
        }
    } else if (ctot[3] > 0) {
        ktemp = ctot[1] + 2 + ctot[2];
        dgemm_((char *)"N", (char *)"N", nl, k, &ctot[3], &c_b12, &u2[ktemp * u2_dim1 + 1], ldu2,
               &q[ktemp + q_dim1], ldq, &c_b25, &u[u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1);
    } else {
        dlacpy_((char *)"F", nl, k, &u2[u2_offset], ldu2, &u[u_offset], ldu, (ftnlen)1);
    }
    dcopy_(k, &q[q_dim1 + 1], ldq, &u[nlp1 + u_dim1], ldu);
    ktemp = ctot[1] + 2;
    ctemp = ctot[2] + ctot[3];
    dgemm_((char *)"N", (char *)"N", nr, k, &ctemp, &c_b12, &u2[nlp2 + ktemp * u2_dim1], ldu2, &q[ktemp + q_dim1],
           ldq, &c_b25, &u[nlp2 + u_dim1], ldu, (ftnlen)1, (ftnlen)1);
L100:
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
        temp = dnrm2_(k, &vt[i__ * vt_dim1 + 1], &c__1);
        q[i__ + q_dim1] = vt[i__ * vt_dim1 + 1] / temp;
        i__2 = *k;
        for (j = 2; j <= i__2; ++j) {
            jc = idxc[j];
            q[i__ + j * q_dim1] = vt[jc + i__ * vt_dim1] / temp;
        }
    }
    if (*k == 2) {
        dgemm_((char *)"N", (char *)"N", k, &m, k, &c_b12, &q[q_offset], ldq, &vt2[vt2_offset], ldvt2, &c_b25,
               &vt[vt_offset], ldvt, (ftnlen)1, (ftnlen)1);
        return 0;
    }
    ktemp = ctot[1] + 1;
    dgemm_((char *)"N", (char *)"N", k, &nlp1, &ktemp, &c_b12, &q[q_dim1 + 1], ldq, &vt2[vt2_dim1 + 1], ldvt2,
           &c_b25, &vt[vt_dim1 + 1], ldvt, (ftnlen)1, (ftnlen)1);
    ktemp = ctot[1] + 2 + ctot[2];
    if (ktemp <= *ldvt2) {
        dgemm_((char *)"N", (char *)"N", k, &nlp1, &ctot[3], &c_b12, &q[ktemp * q_dim1 + 1], ldq,
               &vt2[ktemp + vt2_dim1], ldvt2, &c_b12, &vt[vt_dim1 + 1], ldvt, (ftnlen)1, (ftnlen)1);
    }
    ktemp = ctot[1] + 1;
    nrp1 = *nr + *sqre;
    if (ktemp > 1) {
        i__1 = *k;
        for (i__ = 1; i__ <= i__1; ++i__) {
            q[i__ + ktemp * q_dim1] = q[i__ + q_dim1];
        }
        i__1 = m;
        for (i__ = nlp2; i__ <= i__1; ++i__) {
            vt2[ktemp + i__ * vt2_dim1] = vt2[i__ * vt2_dim1 + 1];
        }
    }
    ctemp = ctot[2] + 1 + ctot[3];
    dgemm_((char *)"N", (char *)"N", k, &nrp1, &ctemp, &c_b12, &q[ktemp * q_dim1 + 1], ldq,
           &vt2[ktemp + nlp2 * vt2_dim1], ldvt2, &c_b25, &vt[nlp2 * vt_dim1 + 1], ldvt, (ftnlen)1,
           (ftnlen)1);
    return 0;
}
#ifdef __cplusplus
}
#endif
