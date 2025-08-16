#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
int dorghr_(integer *n, integer *ilo, integer *ihi, doublereal *a, integer *lda, doublereal *tau,
            doublereal *work, integer *lwork, integer *info)
{
    integer a_dim1, a_offset, i__1, i__2;
    integer i__, j, nb, nh, iinfo;
    extern int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    extern int dorgqr_(integer *, integer *, integer *, doublereal *, integer *, doublereal *,
                       doublereal *, integer *, integer *);
    integer lwkopt;
    logical lquery;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    *info = 0;
    nh = *ihi - *ilo;
    lquery = *lwork == -1;
    if (*n < 0) {
        *info = -1;
    } else if (*ilo < 1 || *ilo > max(1, *n)) {
        *info = -2;
    } else if (*ihi < min(*ilo, *n) || *ihi > *n) {
        *info = -3;
    } else if (*lda < max(1, *n)) {
        *info = -5;
    } else if (*lwork < max(1, nh) && !lquery) {
        *info = -8;
    }
    if (*info == 0) {
        nb = ilaenv_(&c__1, (char *)"DORGQR", (char *)" ", &nh, &nh, &nh, &c_n1, (ftnlen)6, (ftnlen)1);
        lwkopt = max(1, nh) * nb;
        work[1] = (doublereal)lwkopt;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DORGHR", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    if (*n == 0) {
        work[1] = 1.;
        return 0;
    }
    i__1 = *ilo + 1;
    for (j = *ihi; j >= i__1; --j) {
        i__2 = j - 1;
        for (i__ = 1; i__ <= i__2; ++i__) {
            a[i__ + j * a_dim1] = 0.;
        }
        i__2 = *ihi;
        for (i__ = j + 1; i__ <= i__2; ++i__) {
            a[i__ + j * a_dim1] = a[i__ + (j - 1) * a_dim1];
        }
        i__2 = *n;
        for (i__ = *ihi + 1; i__ <= i__2; ++i__) {
            a[i__ + j * a_dim1] = 0.;
        }
    }
    i__1 = *ilo;
    for (j = 1; j <= i__1; ++j) {
        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) {
            a[i__ + j * a_dim1] = 0.;
        }
        a[j + j * a_dim1] = 1.;
    }
    i__1 = *n;
    for (j = *ihi + 1; j <= i__1; ++j) {
        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) {
            a[i__ + j * a_dim1] = 0.;
        }
        a[j + j * a_dim1] = 1.;
    }
    if (nh > 0) {
        dorgqr_(&nh, &nh, &nh, &a[*ilo + 1 + (*ilo + 1) * a_dim1], lda, &tau[*ilo], &work[1], lwork,
                &iinfo);
    }
    work[1] = (doublereal)lwkopt;
    return 0;
}
#ifdef __cplusplus
}
#endif
