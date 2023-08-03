#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b8 = 1.;
int dlasd8_(integer *icompq, integer *k, doublereal *d__, doublereal *z__, doublereal *vf,
            doublereal *vl, doublereal *difl, doublereal *difr, integer *lddifr, doublereal *dsigma,
            doublereal *work, integer *info)
{
    integer difr_dim1, difr_offset, i__1, i__2;
    doublereal d__1, d__2;
    double sqrt(doublereal), d_lmp_sign(doublereal *, doublereal *);
    integer i__, j;
    doublereal dj, rho;
    integer iwk1, iwk2, iwk3;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, integer *);
    doublereal temp;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    integer iwk2i, iwk3i;
    doublereal diflj, difrj, dsigj;
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    extern doublereal dlamc3_(doublereal *, doublereal *);
    extern int dlasd4_(integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *,
                       doublereal *, doublereal *, integer *),
        dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *,
                doublereal *, integer *, integer *, ftnlen),
        dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                ftnlen),
        xerbla_(char *, integer *, ftnlen);
    doublereal dsigjp;
    --d__;
    --z__;
    --vf;
    --vl;
    --difl;
    difr_dim1 = *lddifr;
    difr_offset = 1 + difr_dim1;
    difr -= difr_offset;
    --dsigma;
    --work;
    *info = 0;
    if (*icompq < 0 || *icompq > 1) {
        *info = -1;
    } else if (*k < 1) {
        *info = -2;
    } else if (*lddifr < *k) {
        *info = -9;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLASD8", &i__1, (ftnlen)6);
        return 0;
    }
    if (*k == 1) {
        d__[1] = abs(z__[1]);
        difl[1] = d__[1];
        if (*icompq == 1) {
            difl[2] = 1.;
            difr[(difr_dim1 << 1) + 1] = 1.;
        }
        return 0;
    }
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
        dsigma[i__] = dlamc3_(&dsigma[i__], &dsigma[i__]) - dsigma[i__];
    }
    iwk1 = 1;
    iwk2 = iwk1 + *k;
    iwk3 = iwk2 + *k;
    iwk2i = iwk2 - 1;
    iwk3i = iwk3 - 1;
    rho = dnrm2_(k, &z__[1], &c__1);
    dlascl_((char *)"G", &c__0, &c__0, &rho, &c_b8, k, &c__1, &z__[1], k, info, (ftnlen)1);
    rho *= rho;
    dlaset_((char *)"A", k, &c__1, &c_b8, &c_b8, &work[iwk3], k, (ftnlen)1);
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
        dlasd4_(k, &j, &dsigma[1], &z__[1], &work[iwk1], &rho, &d__[j], &work[iwk2], info);
        if (*info != 0) {
            return 0;
        }
        work[iwk3i + j] = work[iwk3i + j] * work[j] * work[iwk2i + j];
        difl[j] = -work[j];
        difr[j + difr_dim1] = -work[j + 1];
        i__2 = j - 1;
        for (i__ = 1; i__ <= i__2; ++i__) {
            work[iwk3i + i__] = work[iwk3i + i__] * work[i__] * work[iwk2i + i__] /
                                (dsigma[i__] - dsigma[j]) / (dsigma[i__] + dsigma[j]);
        }
        i__2 = *k;
        for (i__ = j + 1; i__ <= i__2; ++i__) {
            work[iwk3i + i__] = work[iwk3i + i__] * work[i__] * work[iwk2i + i__] /
                                (dsigma[i__] - dsigma[j]) / (dsigma[i__] + dsigma[j]);
        }
    }
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
        d__2 = sqrt((d__1 = work[iwk3i + i__], abs(d__1)));
        z__[i__] = d_lmp_sign(&d__2, &z__[i__]);
    }
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
        diflj = difl[j];
        dj = d__[j];
        dsigj = -dsigma[j];
        if (j < *k) {
            difrj = -difr[j + difr_dim1];
            dsigjp = -dsigma[j + 1];
        }
        work[j] = -z__[j] / diflj / (dsigma[j] + dj);
        i__2 = j - 1;
        for (i__ = 1; i__ <= i__2; ++i__) {
            work[i__] = z__[i__] / (dlamc3_(&dsigma[i__], &dsigj) - diflj) / (dsigma[i__] + dj);
        }
        i__2 = *k;
        for (i__ = j + 1; i__ <= i__2; ++i__) {
            work[i__] = z__[i__] / (dlamc3_(&dsigma[i__], &dsigjp) + difrj) / (dsigma[i__] + dj);
        }
        temp = dnrm2_(k, &work[1], &c__1);
        work[iwk2i + j] = ddot_(k, &work[1], &c__1, &vf[1], &c__1) / temp;
        work[iwk3i + j] = ddot_(k, &work[1], &c__1, &vl[1], &c__1) / temp;
        if (*icompq == 1) {
            difr[j + (difr_dim1 << 1)] = temp;
        }
    }
    dcopy_(k, &work[iwk2], &c__1, &vf[1], &c__1);
    dcopy_(k, &work[iwk3], &c__1, &vl[1], &c__1);
    return 0;
}
#ifdef __cplusplus
}
#endif
