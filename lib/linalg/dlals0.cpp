#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b5 = -1.;
static integer c__1 = 1;
static doublereal c_b11 = 1.;
static doublereal c_b13 = 0.;
static integer c__0 = 0;
int dlals0_(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *nrhs, doublereal *b,
            integer *ldb, doublereal *bx, integer *ldbx, integer *perm, integer *givptr,
            integer *givcol, integer *ldgcol, doublereal *givnum, integer *ldgnum,
            doublereal *poles, doublereal *difl, doublereal *difr, doublereal *z__, integer *k,
            doublereal *c__, doublereal *s, doublereal *work, integer *info)
{
    integer givcol_dim1, givcol_offset, b_dim1, b_offset, bx_dim1, bx_offset, difr_dim1,
        difr_offset, givnum_dim1, givnum_offset, poles_dim1, poles_offset, i__1, i__2;
    doublereal d__1;
    integer i__, j, m, n;
    doublereal dj;
    integer nlp1;
    doublereal temp;
    extern int drot_(integer *, doublereal *, integer *, doublereal *, integer *, doublereal *,
                     doublereal *);
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern int dscal_(integer *, doublereal *, doublereal *, integer *);
    doublereal diflj, difrj, dsigj;
    extern int dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *,
                      doublereal *, integer *, doublereal *, doublereal *, integer *, ftnlen),
        dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    extern doublereal dlamc3_(doublereal *, doublereal *);
    extern int dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *,
                       integer *, doublereal *, integer *, integer *, ftnlen),
        dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *,
                ftnlen),
        xerbla_(char *, integer *, ftnlen);
    doublereal dsigjp;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    bx_dim1 = *ldbx;
    bx_offset = 1 + bx_dim1;
    bx -= bx_offset;
    --perm;
    givcol_dim1 = *ldgcol;
    givcol_offset = 1 + givcol_dim1;
    givcol -= givcol_offset;
    difr_dim1 = *ldgnum;
    difr_offset = 1 + difr_dim1;
    difr -= difr_offset;
    poles_dim1 = *ldgnum;
    poles_offset = 1 + poles_dim1;
    poles -= poles_offset;
    givnum_dim1 = *ldgnum;
    givnum_offset = 1 + givnum_dim1;
    givnum -= givnum_offset;
    --difl;
    --z__;
    --work;
    *info = 0;
    n = *nl + *nr + 1;
    if (*icompq < 0 || *icompq > 1) {
        *info = -1;
    } else if (*nl < 1) {
        *info = -2;
    } else if (*nr < 1) {
        *info = -3;
    } else if (*sqre < 0 || *sqre > 1) {
        *info = -4;
    } else if (*nrhs < 1) {
        *info = -5;
    } else if (*ldb < n) {
        *info = -7;
    } else if (*ldbx < n) {
        *info = -9;
    } else if (*givptr < 0) {
        *info = -11;
    } else if (*ldgcol < n) {
        *info = -13;
    } else if (*ldgnum < n) {
        *info = -15;
    } else if (*k < 1) {
        *info = -20;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLALS0", &i__1, (ftnlen)6);
        return 0;
    }
    m = n + *sqre;
    nlp1 = *nl + 1;
    if (*icompq == 0) {
        i__1 = *givptr;
        for (i__ = 1; i__ <= i__1; ++i__) {
            drot_(nrhs, &b[givcol[i__ + (givcol_dim1 << 1)] + b_dim1], ldb,
                  &b[givcol[i__ + givcol_dim1] + b_dim1], ldb, &givnum[i__ + (givnum_dim1 << 1)],
                  &givnum[i__ + givnum_dim1]);
        }
        dcopy_(nrhs, &b[nlp1 + b_dim1], ldb, &bx[bx_dim1 + 1], ldbx);
        i__1 = n;
        for (i__ = 2; i__ <= i__1; ++i__) {
            dcopy_(nrhs, &b[perm[i__] + b_dim1], ldb, &bx[i__ + bx_dim1], ldbx);
        }
        if (*k == 1) {
            dcopy_(nrhs, &bx[bx_offset], ldbx, &b[b_offset], ldb);
            if (z__[1] < 0.) {
                dscal_(nrhs, &c_b5, &b[b_offset], ldb);
            }
        } else {
            i__1 = *k;
            for (j = 1; j <= i__1; ++j) {
                diflj = difl[j];
                dj = poles[j + poles_dim1];
                dsigj = -poles[j + (poles_dim1 << 1)];
                if (j < *k) {
                    difrj = -difr[j + difr_dim1];
                    dsigjp = -poles[j + 1 + (poles_dim1 << 1)];
                }
                if (z__[j] == 0. || poles[j + (poles_dim1 << 1)] == 0.) {
                    work[j] = 0.;
                } else {
                    work[j] = -poles[j + (poles_dim1 << 1)] * z__[j] / diflj /
                              (poles[j + (poles_dim1 << 1)] + dj);
                }
                i__2 = j - 1;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    if (z__[i__] == 0. || poles[i__ + (poles_dim1 << 1)] == 0.) {
                        work[i__] = 0.;
                    } else {
                        work[i__] = poles[i__ + (poles_dim1 << 1)] * z__[i__] /
                                    (dlamc3_(&poles[i__ + (poles_dim1 << 1)], &dsigj) - diflj) /
                                    (poles[i__ + (poles_dim1 << 1)] + dj);
                    }
                }
                i__2 = *k;
                for (i__ = j + 1; i__ <= i__2; ++i__) {
                    if (z__[i__] == 0. || poles[i__ + (poles_dim1 << 1)] == 0.) {
                        work[i__] = 0.;
                    } else {
                        work[i__] = poles[i__ + (poles_dim1 << 1)] * z__[i__] /
                                    (dlamc3_(&poles[i__ + (poles_dim1 << 1)], &dsigjp) + difrj) /
                                    (poles[i__ + (poles_dim1 << 1)] + dj);
                    }
                }
                work[1] = -1.;
                temp = dnrm2_(k, &work[1], &c__1);
                dgemv_((char *)"T", k, nrhs, &c_b11, &bx[bx_offset], ldbx, &work[1], &c__1, &c_b13,
                       &b[j + b_dim1], ldb, (ftnlen)1);
                dlascl_((char *)"G", &c__0, &c__0, &temp, &c_b11, &c__1, nrhs, &b[j + b_dim1], ldb, info,
                        (ftnlen)1);
            }
        }
        if (*k < max(m, n)) {
            i__1 = n - *k;
            dlacpy_((char *)"A", &i__1, nrhs, &bx[*k + 1 + bx_dim1], ldbx, &b[*k + 1 + b_dim1], ldb,
                    (ftnlen)1);
        }
    } else {
        if (*k == 1) {
            dcopy_(nrhs, &b[b_offset], ldb, &bx[bx_offset], ldbx);
        } else {
            i__1 = *k;
            for (j = 1; j <= i__1; ++j) {
                dsigj = poles[j + (poles_dim1 << 1)];
                if (z__[j] == 0.) {
                    work[j] = 0.;
                } else {
                    work[j] = -z__[j] / difl[j] / (dsigj + poles[j + poles_dim1]) /
                              difr[j + (difr_dim1 << 1)];
                }
                i__2 = j - 1;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    if (z__[j] == 0.) {
                        work[i__] = 0.;
                    } else {
                        d__1 = -poles[i__ + 1 + (poles_dim1 << 1)];
                        work[i__] = z__[j] / (dlamc3_(&dsigj, &d__1) - difr[i__ + difr_dim1]) /
                                    (dsigj + poles[i__ + poles_dim1]) /
                                    difr[i__ + (difr_dim1 << 1)];
                    }
                }
                i__2 = *k;
                for (i__ = j + 1; i__ <= i__2; ++i__) {
                    if (z__[j] == 0.) {
                        work[i__] = 0.;
                    } else {
                        d__1 = -poles[i__ + (poles_dim1 << 1)];
                        work[i__] = z__[j] / (dlamc3_(&dsigj, &d__1) - difl[i__]) /
                                    (dsigj + poles[i__ + poles_dim1]) /
                                    difr[i__ + (difr_dim1 << 1)];
                    }
                }
                dgemv_((char *)"T", k, nrhs, &c_b11, &b[b_offset], ldb, &work[1], &c__1, &c_b13,
                       &bx[j + bx_dim1], ldbx, (ftnlen)1);
            }
        }
        if (*sqre == 1) {
            dcopy_(nrhs, &b[m + b_dim1], ldb, &bx[m + bx_dim1], ldbx);
            drot_(nrhs, &bx[bx_dim1 + 1], ldbx, &bx[m + bx_dim1], ldbx, c__, s);
        }
        if (*k < max(m, n)) {
            i__1 = n - *k;
            dlacpy_((char *)"A", &i__1, nrhs, &b[*k + 1 + b_dim1], ldb, &bx[*k + 1 + bx_dim1], ldbx,
                    (ftnlen)1);
        }
        dcopy_(nrhs, &bx[bx_dim1 + 1], ldbx, &b[nlp1 + b_dim1], ldb);
        if (*sqre == 1) {
            dcopy_(nrhs, &bx[m + bx_dim1], ldbx, &b[m + b_dim1], ldb);
        }
        i__1 = n;
        for (i__ = 2; i__ <= i__1; ++i__) {
            dcopy_(nrhs, &bx[i__ + bx_dim1], ldbx, &b[perm[i__] + b_dim1], ldb);
        }
        for (i__ = *givptr; i__ >= 1; --i__) {
            d__1 = -givnum[i__ + givnum_dim1];
            drot_(nrhs, &b[givcol[i__ + (givcol_dim1 << 1)] + b_dim1], ldb,
                  &b[givcol[i__ + givcol_dim1] + b_dim1], ldb, &givnum[i__ + (givnum_dim1 << 1)],
                  &d__1);
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
