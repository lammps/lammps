#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublecomplex c_b1 = {1., 0.};
static integer c__1 = 1;
int zgetrf2_(integer *m, integer *n, doublecomplex *a, integer *lda, integer *ipiv, integer *info)
{
    integer a_dim1, a_offset, i__1, i__2;
    doublecomplex z__1;
    double z_lmp_abs(doublecomplex *);
    void z_lmp_div(doublecomplex *, doublecomplex *, doublecomplex *);
    integer i__, n1, n2;
    doublecomplex temp;
    integer iinfo;
    doublereal sfmin;
    extern int zscal_(integer *, doublecomplex *, doublecomplex *, integer *),
        zgemm_(char *, char *, integer *, integer *, integer *, doublecomplex *, doublecomplex *,
               integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *,
               ftnlen, ftnlen),
        ztrsm_(char *, char *, char *, char *, integer *, integer *, doublecomplex *,
               doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, ftnlen, ftnlen,
               ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    extern int xerbla_(char *, integer *, ftnlen);
    extern integer izamax_(integer *, doublecomplex *, integer *);
    extern int zlaswp_(integer *, doublecomplex *, integer *, integer *, integer *, integer *,
                       integer *);
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
        xerbla_((char *)"ZGETRF2", &i__1, (ftnlen)7);
        return 0;
    }
    if (*m == 0 || *n == 0) {
        return 0;
    }
    if (*m == 1) {
        ipiv[1] = 1;
        i__1 = a_dim1 + 1;
        if (a[i__1].r == 0. && a[i__1].i == 0.) {
            *info = 1;
        }
    } else if (*n == 1) {
        sfmin = dlamch_((char *)"S", (ftnlen)1);
        i__ = izamax_(m, &a[a_dim1 + 1], &c__1);
        ipiv[1] = i__;
        i__1 = i__ + a_dim1;
        if (a[i__1].r != 0. || a[i__1].i != 0.) {
            if (i__ != 1) {
                i__1 = a_dim1 + 1;
                temp.r = a[i__1].r, temp.i = a[i__1].i;
                i__1 = a_dim1 + 1;
                i__2 = i__ + a_dim1;
                a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
                i__1 = i__ + a_dim1;
                a[i__1].r = temp.r, a[i__1].i = temp.i;
            }
            if (z_lmp_abs(&a[a_dim1 + 1]) >= sfmin) {
                i__1 = *m - 1;
                z_lmp_div(&z__1, &c_b1, &a[a_dim1 + 1]);
                zscal_(&i__1, &z__1, &a[a_dim1 + 2], &c__1);
            } else {
                i__1 = *m - 1;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    i__2 = i__ + 1 + a_dim1;
                    z_lmp_div(&z__1, &a[i__ + 1 + a_dim1], &a[a_dim1 + 1]);
                    a[i__2].r = z__1.r, a[i__2].i = z__1.i;
                }
            }
        } else {
            *info = 1;
        }
    } else {
        n1 = min(*m, *n) / 2;
        n2 = *n - n1;
        zgetrf2_(m, &n1, &a[a_offset], lda, &ipiv[1], &iinfo);
        if (*info == 0 && iinfo > 0) {
            *info = iinfo;
        }
        zlaswp_(&n2, &a[(n1 + 1) * a_dim1 + 1], lda, &c__1, &n1, &ipiv[1], &c__1);
        ztrsm_((char *)"L", (char *)"L", (char *)"N", (char *)"U", &n1, &n2, &c_b1, &a[a_offset], lda, &a[(n1 + 1) * a_dim1 + 1],
               lda, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
        i__1 = *m - n1;
        z__1.r = -1., z__1.i = -0.;
        zgemm_((char *)"N", (char *)"N", &i__1, &n2, &n1, &z__1, &a[n1 + 1 + a_dim1], lda,
               &a[(n1 + 1) * a_dim1 + 1], lda, &c_b1, &a[n1 + 1 + (n1 + 1) * a_dim1], lda,
               (ftnlen)1, (ftnlen)1);
        i__1 = *m - n1;
        zgetrf2_(&i__1, &n2, &a[n1 + 1 + (n1 + 1) * a_dim1], lda, &ipiv[n1 + 1], &iinfo);
        if (*info == 0 && iinfo > 0) {
            *info = iinfo + n1;
        }
        i__1 = min(*m, *n);
        for (i__ = n1 + 1; i__ <= i__1; ++i__) {
            ipiv[i__] += n1;
        }
        i__1 = n1 + 1;
        i__2 = min(*m, *n);
        zlaswp_(&n1, &a[a_dim1 + 1], lda, &i__1, &i__2, &ipiv[1], &c__1);
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
