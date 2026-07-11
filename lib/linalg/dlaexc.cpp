#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c__4 = 4;
static logical c_false = FALSE_;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__3 = 3;
int dlaexc_(logical *wantq, integer *n, doublereal *t, integer *ldt, doublereal *q, integer *ldq,
            integer *j1, integer *n1, integer *n2, doublereal *work, integer *info)
{
    integer q_dim1, q_offset, t_dim1, t_offset, i__1;
    doublereal d__1, d__2, d__3;
    doublereal d__[16];
    integer k;
    doublereal u[3], x[4];
    integer j2, j3, j4;
    doublereal u1[3], u2[3];
    integer nd;
    doublereal cs, t11, t22, t33, sn, wi1, wi2, wr1, wr2, eps, tau, tau1, tau2;
    integer ierr;
    doublereal temp;
    extern int drot_(integer *, doublereal *, integer *, doublereal *, integer *, doublereal *,
                     doublereal *);
    doublereal scale, dnorm, xnorm;
    extern int dlanv2_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
                       doublereal *, doublereal *, doublereal *, doublereal *, doublereal *),
        dlasy2_(logical *, logical *, integer *, integer *, integer *, doublereal *, integer *,
                doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                integer *, doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen),
        dlange_(char *, integer *, integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern int dlarfg_(integer *, doublereal *, doublereal *, integer *, doublereal *),
        dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *,
                ftnlen),
        dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *),
        dlarfx_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                doublereal *, ftnlen);
    doublereal thresh, smlnum;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --work;
    *info = 0;
    if (*n == 0 || *n1 == 0 || *n2 == 0) {
        return 0;
    }
    if (*j1 + *n1 > *n) {
        return 0;
    }
    j2 = *j1 + 1;
    j3 = *j1 + 2;
    j4 = *j1 + 3;
    if (*n1 == 1 && *n2 == 1) {
        t11 = t[*j1 + *j1 * t_dim1];
        t22 = t[j2 + j2 * t_dim1];
        d__1 = t22 - t11;
        dlartg_(&t[*j1 + j2 * t_dim1], &d__1, &cs, &sn, &temp);
        if (j3 <= *n) {
            i__1 = *n - *j1 - 1;
            drot_(&i__1, &t[*j1 + j3 * t_dim1], ldt, &t[j2 + j3 * t_dim1], ldt, &cs, &sn);
        }
        i__1 = *j1 - 1;
        drot_(&i__1, &t[*j1 * t_dim1 + 1], &c__1, &t[j2 * t_dim1 + 1], &c__1, &cs, &sn);
        t[*j1 + *j1 * t_dim1] = t22;
        t[j2 + j2 * t_dim1] = t11;
        if (*wantq) {
            drot_(n, &q[*j1 * q_dim1 + 1], &c__1, &q[j2 * q_dim1 + 1], &c__1, &cs, &sn);
        }
    } else {
        nd = *n1 + *n2;
        dlacpy_((char *)"F", &nd, &nd, &t[*j1 + *j1 * t_dim1], ldt, d__, &c__4, (ftnlen)1);
        dnorm = dlange_((char *)"Max", &nd, &nd, d__, &c__4, &work[1], (ftnlen)3);
        eps = dlamch_((char *)"P", (ftnlen)1);
        smlnum = dlamch_((char *)"S", (ftnlen)1) / eps;
        d__1 = eps * 10. * dnorm;
        thresh = max(d__1, smlnum);
        dlasy2_(&c_false, &c_false, &c_n1, n1, n2, d__, &c__4, &d__[*n1 + 1 + (*n1 + 1 << 2) - 5],
                &c__4, &d__[(*n1 + 1 << 2) - 4], &c__4, &scale, x, &c__2, &xnorm, &ierr);
        k = *n1 + *n1 + *n2 - 3;
        switch (k) {
            case 1:
                goto L10;
            case 2:
                goto L20;
            case 3:
                goto L30;
        }
    L10:
        u[0] = scale;
        u[1] = x[0];
        u[2] = x[2];
        dlarfg_(&c__3, &u[2], u, &c__1, &tau);
        u[2] = 1.;
        t11 = t[*j1 + *j1 * t_dim1];
        dlarfx_((char *)"L", &c__3, &c__3, u, &tau, d__, &c__4, &work[1], (ftnlen)1);
        dlarfx_((char *)"R", &c__3, &c__3, u, &tau, d__, &c__4, &work[1], (ftnlen)1);
        d__2 = abs(d__[2]), d__3 = abs(d__[6]), d__2 = max(d__2, d__3),
        d__3 = (d__1 = d__[10] - t11, abs(d__1));
        if (max(d__2, d__3) > thresh) {
            goto L50;
        }
        i__1 = *n - *j1 + 1;
        dlarfx_((char *)"L", &c__3, &i__1, u, &tau, &t[*j1 + *j1 * t_dim1], ldt, &work[1], (ftnlen)1);
        dlarfx_((char *)"R", &j2, &c__3, u, &tau, &t[*j1 * t_dim1 + 1], ldt, &work[1], (ftnlen)1);
        t[j3 + *j1 * t_dim1] = 0.;
        t[j3 + j2 * t_dim1] = 0.;
        t[j3 + j3 * t_dim1] = t11;
        if (*wantq) {
            dlarfx_((char *)"R", n, &c__3, u, &tau, &q[*j1 * q_dim1 + 1], ldq, &work[1], (ftnlen)1);
        }
        goto L40;
    L20:
        u[0] = -x[0];
        u[1] = -x[1];
        u[2] = scale;
        dlarfg_(&c__3, u, &u[1], &c__1, &tau);
        u[0] = 1.;
        t33 = t[j3 + j3 * t_dim1];
        dlarfx_((char *)"L", &c__3, &c__3, u, &tau, d__, &c__4, &work[1], (ftnlen)1);
        dlarfx_((char *)"R", &c__3, &c__3, u, &tau, d__, &c__4, &work[1], (ftnlen)1);
        d__2 = abs(d__[1]), d__3 = abs(d__[2]), d__2 = max(d__2, d__3),
        d__3 = (d__1 = d__[0] - t33, abs(d__1));
        if (max(d__2, d__3) > thresh) {
            goto L50;
        }
        dlarfx_((char *)"R", &j3, &c__3, u, &tau, &t[*j1 * t_dim1 + 1], ldt, &work[1], (ftnlen)1);
        i__1 = *n - *j1;
        dlarfx_((char *)"L", &c__3, &i__1, u, &tau, &t[*j1 + j2 * t_dim1], ldt, &work[1], (ftnlen)1);
        t[*j1 + *j1 * t_dim1] = t33;
        t[j2 + *j1 * t_dim1] = 0.;
        t[j3 + *j1 * t_dim1] = 0.;
        if (*wantq) {
            dlarfx_((char *)"R", n, &c__3, u, &tau, &q[*j1 * q_dim1 + 1], ldq, &work[1], (ftnlen)1);
        }
        goto L40;
    L30:
        u1[0] = -x[0];
        u1[1] = -x[1];
        u1[2] = scale;
        dlarfg_(&c__3, u1, &u1[1], &c__1, &tau1);
        u1[0] = 1.;
        temp = -tau1 * (x[2] + u1[1] * x[3]);
        u2[0] = -temp * u1[1] - x[3];
        u2[1] = -temp * u1[2];
        u2[2] = scale;
        dlarfg_(&c__3, u2, &u2[1], &c__1, &tau2);
        u2[0] = 1.;
        dlarfx_((char *)"L", &c__3, &c__4, u1, &tau1, d__, &c__4, &work[1], (ftnlen)1);
        dlarfx_((char *)"R", &c__4, &c__3, u1, &tau1, d__, &c__4, &work[1], (ftnlen)1);
        dlarfx_((char *)"L", &c__3, &c__4, u2, &tau2, &d__[1], &c__4, &work[1], (ftnlen)1);
        dlarfx_((char *)"R", &c__4, &c__3, u2, &tau2, &d__[4], &c__4, &work[1], (ftnlen)1);
        d__1 = abs(d__[2]), d__2 = abs(d__[6]), d__1 = max(d__1, d__2), d__2 = abs(d__[3]),
        d__1 = max(d__1, d__2), d__2 = abs(d__[7]);
        if (max(d__1, d__2) > thresh) {
            goto L50;
        }
        i__1 = *n - *j1 + 1;
        dlarfx_((char *)"L", &c__3, &i__1, u1, &tau1, &t[*j1 + *j1 * t_dim1], ldt, &work[1], (ftnlen)1);
        dlarfx_((char *)"R", &j4, &c__3, u1, &tau1, &t[*j1 * t_dim1 + 1], ldt, &work[1], (ftnlen)1);
        i__1 = *n - *j1 + 1;
        dlarfx_((char *)"L", &c__3, &i__1, u2, &tau2, &t[j2 + *j1 * t_dim1], ldt, &work[1], (ftnlen)1);
        dlarfx_((char *)"R", &j4, &c__3, u2, &tau2, &t[j2 * t_dim1 + 1], ldt, &work[1], (ftnlen)1);
        t[j3 + *j1 * t_dim1] = 0.;
        t[j3 + j2 * t_dim1] = 0.;
        t[j4 + *j1 * t_dim1] = 0.;
        t[j4 + j2 * t_dim1] = 0.;
        if (*wantq) {
            dlarfx_((char *)"R", n, &c__3, u1, &tau1, &q[*j1 * q_dim1 + 1], ldq, &work[1], (ftnlen)1);
            dlarfx_((char *)"R", n, &c__3, u2, &tau2, &q[j2 * q_dim1 + 1], ldq, &work[1], (ftnlen)1);
        }
    L40:
        if (*n2 == 2) {
            dlanv2_(&t[*j1 + *j1 * t_dim1], &t[*j1 + j2 * t_dim1], &t[j2 + *j1 * t_dim1],
                    &t[j2 + j2 * t_dim1], &wr1, &wi1, &wr2, &wi2, &cs, &sn);
            i__1 = *n - *j1 - 1;
            drot_(&i__1, &t[*j1 + (*j1 + 2) * t_dim1], ldt, &t[j2 + (*j1 + 2) * t_dim1], ldt, &cs,
                  &sn);
            i__1 = *j1 - 1;
            drot_(&i__1, &t[*j1 * t_dim1 + 1], &c__1, &t[j2 * t_dim1 + 1], &c__1, &cs, &sn);
            if (*wantq) {
                drot_(n, &q[*j1 * q_dim1 + 1], &c__1, &q[j2 * q_dim1 + 1], &c__1, &cs, &sn);
            }
        }
        if (*n1 == 2) {
            j3 = *j1 + *n2;
            j4 = j3 + 1;
            dlanv2_(&t[j3 + j3 * t_dim1], &t[j3 + j4 * t_dim1], &t[j4 + j3 * t_dim1],
                    &t[j4 + j4 * t_dim1], &wr1, &wi1, &wr2, &wi2, &cs, &sn);
            if (j3 + 2 <= *n) {
                i__1 = *n - j3 - 1;
                drot_(&i__1, &t[j3 + (j3 + 2) * t_dim1], ldt, &t[j4 + (j3 + 2) * t_dim1], ldt, &cs,
                      &sn);
            }
            i__1 = j3 - 1;
            drot_(&i__1, &t[j3 * t_dim1 + 1], &c__1, &t[j4 * t_dim1 + 1], &c__1, &cs, &sn);
            if (*wantq) {
                drot_(n, &q[j3 * q_dim1 + 1], &c__1, &q[j4 * q_dim1 + 1], &c__1, &cs, &sn);
            }
        }
    }
    return 0;
L50:
    *info = 1;
    return 0;
}
#ifdef __cplusplus
}
#endif
