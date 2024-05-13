#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
int dlasdq_(char *uplo, integer *sqre, integer *n, integer *ncvt, integer *nru, integer *ncc,
            doublereal *d__, doublereal *e, doublereal *vt, integer *ldvt, doublereal *u,
            integer *ldu, doublereal *c__, integer *ldc, doublereal *work, integer *info,
            ftnlen uplo_len)
{
    integer c_dim1, c_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2;
    integer i__, j;
    doublereal r__, cs, sn;
    integer np1, isub;
    doublereal smin;
    integer sqre1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int dlasr_(char *, char *, char *, integer *, integer *, doublereal *, doublereal *,
                      doublereal *, integer *, ftnlen, ftnlen, ftnlen),
        dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer iuplo;
    extern int dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *),
        xerbla_(char *, integer *, ftnlen),
        dbdsqr_(char *, integer *, integer *, integer *, integer *, doublereal *, doublereal *,
                doublereal *, integer *, doublereal *, integer *, doublereal *, integer *,
                doublereal *, integer *, ftnlen);
    logical rotate;
    --d__;
    --e;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    *info = 0;
    iuplo = 0;
    if (lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1)) {
        iuplo = 1;
    }
    if (lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        iuplo = 2;
    }
    if (iuplo == 0) {
        *info = -1;
    } else if (*sqre < 0 || *sqre > 1) {
        *info = -2;
    } else if (*n < 0) {
        *info = -3;
    } else if (*ncvt < 0) {
        *info = -4;
    } else if (*nru < 0) {
        *info = -5;
    } else if (*ncc < 0) {
        *info = -6;
    } else if (*ncvt == 0 && *ldvt < 1 || *ncvt > 0 && *ldvt < max(1, *n)) {
        *info = -10;
    } else if (*ldu < max(1, *nru)) {
        *info = -12;
    } else if (*ncc == 0 && *ldc < 1 || *ncc > 0 && *ldc < max(1, *n)) {
        *info = -14;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLASDQ", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    rotate = *ncvt > 0 || *nru > 0 || *ncc > 0;
    np1 = *n + 1;
    sqre1 = *sqre;
    if (iuplo == 1 && sqre1 == 1) {
        i__1 = *n - 1;
        for (i__ = 1; i__ <= i__1; ++i__) {
            dlartg_(&d__[i__], &e[i__], &cs, &sn, &r__);
            d__[i__] = r__;
            e[i__] = sn * d__[i__ + 1];
            d__[i__ + 1] = cs * d__[i__ + 1];
            if (rotate) {
                work[i__] = cs;
                work[*n + i__] = sn;
            }
        }
        dlartg_(&d__[*n], &e[*n], &cs, &sn, &r__);
        d__[*n] = r__;
        e[*n] = 0.;
        if (rotate) {
            work[*n] = cs;
            work[*n + *n] = sn;
        }
        iuplo = 2;
        sqre1 = 0;
        if (*ncvt > 0) {
            dlasr_((char *)"L", (char *)"V", (char *)"F", &np1, ncvt, &work[1], &work[np1], &vt[vt_offset], ldvt, (ftnlen)1,
                   (ftnlen)1, (ftnlen)1);
        }
    }
    if (iuplo == 2) {
        i__1 = *n - 1;
        for (i__ = 1; i__ <= i__1; ++i__) {
            dlartg_(&d__[i__], &e[i__], &cs, &sn, &r__);
            d__[i__] = r__;
            e[i__] = sn * d__[i__ + 1];
            d__[i__ + 1] = cs * d__[i__ + 1];
            if (rotate) {
                work[i__] = cs;
                work[*n + i__] = sn;
            }
        }
        if (sqre1 == 1) {
            dlartg_(&d__[*n], &e[*n], &cs, &sn, &r__);
            d__[*n] = r__;
            if (rotate) {
                work[*n] = cs;
                work[*n + *n] = sn;
            }
        }
        if (*nru > 0) {
            if (sqre1 == 0) {
                dlasr_((char *)"R", (char *)"V", (char *)"F", nru, n, &work[1], &work[np1], &u[u_offset], ldu, (ftnlen)1,
                       (ftnlen)1, (ftnlen)1);
            } else {
                dlasr_((char *)"R", (char *)"V", (char *)"F", nru, &np1, &work[1], &work[np1], &u[u_offset], ldu, (ftnlen)1,
                       (ftnlen)1, (ftnlen)1);
            }
        }
        if (*ncc > 0) {
            if (sqre1 == 0) {
                dlasr_((char *)"L", (char *)"V", (char *)"F", n, ncc, &work[1], &work[np1], &c__[c_offset], ldc, (ftnlen)1,
                       (ftnlen)1, (ftnlen)1);
            } else {
                dlasr_((char *)"L", (char *)"V", (char *)"F", &np1, ncc, &work[1], &work[np1], &c__[c_offset], ldc,
                       (ftnlen)1, (ftnlen)1, (ftnlen)1);
            }
        }
    }
    dbdsqr_((char *)"U", n, ncvt, nru, ncc, &d__[1], &e[1], &vt[vt_offset], ldvt, &u[u_offset], ldu,
            &c__[c_offset], ldc, &work[1], info, (ftnlen)1);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        isub = i__;
        smin = d__[i__];
        i__2 = *n;
        for (j = i__ + 1; j <= i__2; ++j) {
            if (d__[j] < smin) {
                isub = j;
                smin = d__[j];
            }
        }
        if (isub != i__) {
            d__[isub] = d__[i__];
            d__[i__] = smin;
            if (*ncvt > 0) {
                dswap_(ncvt, &vt[isub + vt_dim1], ldvt, &vt[i__ + vt_dim1], ldvt);
            }
            if (*nru > 0) {
                dswap_(nru, &u[isub * u_dim1 + 1], &c__1, &u[i__ * u_dim1 + 1], &c__1);
            }
            if (*ncc > 0) {
                dswap_(ncc, &c__[isub + c_dim1], ldc, &c__[i__ + c_dim1], ldc);
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
