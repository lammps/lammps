#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__2 = 2;
static integer c__1 = 1;
static doublereal c_b24 = 1.;
static doublereal c_b26 = 0.;
int dlaeda_(integer *n, integer *tlvls, integer *curlvl, integer *curpbm, integer *prmptr,
            integer *perm, integer *givptr, integer *givcol, doublereal *givnum, doublereal *q,
            integer *qptr, doublereal *z__, doublereal *ztemp, integer *info)
{
    integer i__1, i__2, i__3;
    integer pow_lmp_ii(integer *, integer *);
    double sqrt(doublereal);
    integer i__, k, mid, ptr;
    extern int drot_(integer *, doublereal *, integer *, doublereal *, integer *, doublereal *,
                     doublereal *);
    integer curr, bsiz1, bsiz2, psiz1, psiz2, zptr1;
    extern int dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *,
                      doublereal *, integer *, doublereal *, doublereal *, integer *, ftnlen),
        dcopy_(integer *, doublereal *, integer *, doublereal *, integer *),
        xerbla_(char *, integer *, ftnlen);
    --ztemp;
    --z__;
    --qptr;
    --q;
    givnum -= 3;
    givcol -= 3;
    --givptr;
    --perm;
    --prmptr;
    *info = 0;
    if (*n < 0) {
        *info = -1;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLAEDA", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    mid = *n / 2 + 1;
    ptr = 1;
    i__1 = *curlvl - 1;
    curr = ptr + *curpbm * pow_lmp_ii(&c__2, curlvl) + pow_lmp_ii(&c__2, &i__1) - 1;
    bsiz1 = (integer)(sqrt((doublereal)(qptr[curr + 1] - qptr[curr])) + .5);
    bsiz2 = (integer)(sqrt((doublereal)(qptr[curr + 2] - qptr[curr + 1])) + .5);
    i__1 = mid - bsiz1 - 1;
    for (k = 1; k <= i__1; ++k) {
        z__[k] = 0.;
    }
    dcopy_(&bsiz1, &q[qptr[curr] + bsiz1 - 1], &bsiz1, &z__[mid - bsiz1], &c__1);
    dcopy_(&bsiz2, &q[qptr[curr + 1]], &bsiz2, &z__[mid], &c__1);
    i__1 = *n;
    for (k = mid + bsiz2; k <= i__1; ++k) {
        z__[k] = 0.;
    }
    ptr = pow_lmp_ii(&c__2, tlvls) + 1;
    i__1 = *curlvl - 1;
    for (k = 1; k <= i__1; ++k) {
        i__2 = *curlvl - k;
        i__3 = *curlvl - k - 1;
        curr = ptr + *curpbm * pow_lmp_ii(&c__2, &i__2) + pow_lmp_ii(&c__2, &i__3) - 1;
        psiz1 = prmptr[curr + 1] - prmptr[curr];
        psiz2 = prmptr[curr + 2] - prmptr[curr + 1];
        zptr1 = mid - psiz1;
        i__2 = givptr[curr + 1] - 1;
        for (i__ = givptr[curr]; i__ <= i__2; ++i__) {
            drot_(&c__1, &z__[zptr1 + givcol[(i__ << 1) + 1] - 1], &c__1,
                  &z__[zptr1 + givcol[(i__ << 1) + 2] - 1], &c__1, &givnum[(i__ << 1) + 1],
                  &givnum[(i__ << 1) + 2]);
        }
        i__2 = givptr[curr + 2] - 1;
        for (i__ = givptr[curr + 1]; i__ <= i__2; ++i__) {
            drot_(&c__1, &z__[mid - 1 + givcol[(i__ << 1) + 1]], &c__1,
                  &z__[mid - 1 + givcol[(i__ << 1) + 2]], &c__1, &givnum[(i__ << 1) + 1],
                  &givnum[(i__ << 1) + 2]);
        }
        psiz1 = prmptr[curr + 1] - prmptr[curr];
        psiz2 = prmptr[curr + 2] - prmptr[curr + 1];
        i__2 = psiz1 - 1;
        for (i__ = 0; i__ <= i__2; ++i__) {
            ztemp[i__ + 1] = z__[zptr1 + perm[prmptr[curr] + i__] - 1];
        }
        i__2 = psiz2 - 1;
        for (i__ = 0; i__ <= i__2; ++i__) {
            ztemp[psiz1 + i__ + 1] = z__[mid + perm[prmptr[curr + 1] + i__] - 1];
        }
        bsiz1 = (integer)(sqrt((doublereal)(qptr[curr + 1] - qptr[curr])) + .5);
        bsiz2 = (integer)(sqrt((doublereal)(qptr[curr + 2] - qptr[curr + 1])) + .5);
        if (bsiz1 > 0) {
            dgemv_((char *)"T", &bsiz1, &bsiz1, &c_b24, &q[qptr[curr]], &bsiz1, &ztemp[1], &c__1, &c_b26,
                   &z__[zptr1], &c__1, (ftnlen)1);
        }
        i__2 = psiz1 - bsiz1;
        dcopy_(&i__2, &ztemp[bsiz1 + 1], &c__1, &z__[zptr1 + bsiz1], &c__1);
        if (bsiz2 > 0) {
            dgemv_((char *)"T", &bsiz2, &bsiz2, &c_b24, &q[qptr[curr + 1]], &bsiz2, &ztemp[psiz1 + 1],
                   &c__1, &c_b26, &z__[mid], &c__1, (ftnlen)1);
        }
        i__2 = psiz2 - bsiz2;
        dcopy_(&i__2, &ztemp[psiz1 + bsiz2 + 1], &c__1, &z__[mid + bsiz2], &c__1);
        i__2 = *tlvls - k;
        ptr += pow_lmp_ii(&c__2, &i__2);
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
