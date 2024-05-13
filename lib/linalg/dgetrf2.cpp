#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static doublereal c_b13 = 1.;
static doublereal c_b16 = -1.;
int dgetrf2_(integer *m, integer *n, doublereal *a, integer *lda, integer *ipiv, integer *info)
{
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;
    integer i__, n1, n2;
    doublereal temp;
    extern int dscal_(integer *, doublereal *, doublereal *, integer *),
        dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
               integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, ftnlen,
               ftnlen);
    integer iinfo;
    doublereal sfmin;
    extern int dtrsm_(char *, char *, char *, char *, integer *, integer *, doublereal *,
                      doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen,
                      ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern int xerbla_(char *, integer *, ftnlen),
        dlaswp_(integer *, doublereal *, integer *, integer *, integer *, integer *, integer *);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    *info = 0;
    if (*m < 0) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*lda < max(1, *m)) {
        *info = -4;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DGETRF2", &i__1, (ftnlen)7);
        return 0;
    }
    if (*m == 0 || *n == 0) {
        return 0;
    }
    if (*m == 1) {
        ipiv[1] = 1;
        if (a[a_dim1 + 1] == 0.) {
            *info = 1;
        }
    } else if (*n == 1) {
        sfmin = dlamch_((char *)"S", (ftnlen)1);
        i__ = idamax_(m, &a[a_dim1 + 1], &c__1);
        ipiv[1] = i__;
        if (a[i__ + a_dim1] != 0.) {
            if (i__ != 1) {
                temp = a[a_dim1 + 1];
                a[a_dim1 + 1] = a[i__ + a_dim1];
                a[i__ + a_dim1] = temp;
            }
            if ((d__1 = a[a_dim1 + 1], abs(d__1)) >= sfmin) {
                i__1 = *m - 1;
                d__1 = 1. / a[a_dim1 + 1];
                dscal_(&i__1, &d__1, &a[a_dim1 + 2], &c__1);
            } else {
                i__1 = *m - 1;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    a[i__ + 1 + a_dim1] /= a[a_dim1 + 1];
                }
            }
        } else {
            *info = 1;
        }
    } else {
        n1 = min(*m, *n) / 2;
        n2 = *n - n1;
        dgetrf2_(m, &n1, &a[a_offset], lda, &ipiv[1], &iinfo);
        if (*info == 0 && iinfo > 0) {
            *info = iinfo;
        }
        dlaswp_(&n2, &a[(n1 + 1) * a_dim1 + 1], lda, &c__1, &n1, &ipiv[1], &c__1);
        dtrsm_((char *)"L", (char *)"L", (char *)"N", (char *)"U", &n1, &n2, &c_b13, &a[a_offset], lda, &a[(n1 + 1) * a_dim1 + 1],
               lda, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
        i__1 = *m - n1;
        dgemm_((char *)"N", (char *)"N", &i__1, &n2, &n1, &c_b16, &a[n1 + 1 + a_dim1], lda,
               &a[(n1 + 1) * a_dim1 + 1], lda, &c_b13, &a[n1 + 1 + (n1 + 1) * a_dim1], lda,
               (ftnlen)1, (ftnlen)1);
        i__1 = *m - n1;
        dgetrf2_(&i__1, &n2, &a[n1 + 1 + (n1 + 1) * a_dim1], lda, &ipiv[n1 + 1], &iinfo);
        if (*info == 0 && iinfo > 0) {
            *info = iinfo + n1;
        }
        i__1 = min(*m, *n);
        for (i__ = n1 + 1; i__ <= i__1; ++i__) {
            ipiv[i__] += n1;
        }
        i__1 = n1 + 1;
        i__2 = min(*m, *n);
        dlaswp_(&n1, &a[a_dim1 + 1], lda, &i__1, &i__2, &ipiv[1], &c__1);
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
