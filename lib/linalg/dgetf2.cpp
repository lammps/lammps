#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static doublereal c_b8 = -1.;
int dgetf2_(integer *m, integer *n, doublereal *a, integer *lda, integer *ipiv, integer *info)
{
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;
    integer i__, j, jp;
    extern int dger_(integer *, integer *, doublereal *, doublereal *, integer *, doublereal *,
                     integer *, doublereal *, integer *),
        dscal_(integer *, doublereal *, doublereal *, integer *);
    doublereal sfmin;
    extern int dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern int xerbla_(char *, integer *, ftnlen);
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
        xerbla_((char *)"DGETF2", &i__1, (ftnlen)6);
        return 0;
    }
    if (*m == 0 || *n == 0) {
        return 0;
    }
    sfmin = dlamch_((char *)"S", (ftnlen)1);
    i__1 = min(*m, *n);
    for (j = 1; j <= i__1; ++j) {
        i__2 = *m - j + 1;
        jp = j - 1 + idamax_(&i__2, &a[j + j * a_dim1], &c__1);
        ipiv[j] = jp;
        if (a[jp + j * a_dim1] != 0.) {
            if (jp != j) {
                dswap_(n, &a[j + a_dim1], lda, &a[jp + a_dim1], lda);
            }
            if (j < *m) {
                if ((d__1 = a[j + j * a_dim1], abs(d__1)) >= sfmin) {
                    i__2 = *m - j;
                    d__1 = 1. / a[j + j * a_dim1];
                    dscal_(&i__2, &d__1, &a[j + 1 + j * a_dim1], &c__1);
                } else {
                    i__2 = *m - j;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        a[j + i__ + j * a_dim1] /= a[j + j * a_dim1];
                    }
                }
            }
        } else if (*info == 0) {
            *info = j;
        }
        if (j < min(*m, *n)) {
            i__2 = *m - j;
            i__3 = *n - j;
            dger_(&i__2, &i__3, &c_b8, &a[j + 1 + j * a_dim1], &c__1, &a[j + (j + 1) * a_dim1], lda,
                  &a[j + 1 + (j + 1) * a_dim1], lda);
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
