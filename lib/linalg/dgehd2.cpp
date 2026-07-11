#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
int dgehd2_(integer *n, integer *ilo, integer *ihi, doublereal *a, integer *lda, doublereal *tau,
            doublereal *work, integer *info)
{
    integer a_dim1, a_offset, i__1, i__2, i__3;
    integer i__;
    extern int dlarfg_(integer *, doublereal *, doublereal *, integer *, doublereal *),
        xerbla_(char *, integer *, ftnlen),
        dlarf1f_(char *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                 integer *, doublereal *, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    *info = 0;
    if (*n < 0) {
        *info = -1;
    } else if (*ilo < 1 || *ilo > max(1, *n)) {
        *info = -2;
    } else if (*ihi < min(*ilo, *n) || *ihi > *n) {
        *info = -3;
    } else if (*lda < max(1, *n)) {
        *info = -5;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DGEHD2", &i__1, (ftnlen)6);
        return 0;
    }
    i__1 = *ihi - 1;
    for (i__ = *ilo; i__ <= i__1; ++i__) {
        i__2 = *ihi - i__;
        i__3 = i__ + 2;
        dlarfg_(&i__2, &a[i__ + 1 + i__ * a_dim1], &a[min(i__3, *n) + i__ * a_dim1], &c__1,
                &tau[i__]);
        i__2 = *ihi - i__;
        dlarf1f_((char *)"R", ihi, &i__2, &a[i__ + 1 + i__ * a_dim1], &c__1, &tau[i__],
                 &a[(i__ + 1) * a_dim1 + 1], lda, &work[1], (ftnlen)1);
        i__2 = *ihi - i__;
        i__3 = *n - i__;
        dlarf1f_((char *)"L", &i__2, &i__3, &a[i__ + 1 + i__ * a_dim1], &c__1, &tau[i__],
                 &a[i__ + 1 + (i__ + 1) * a_dim1], lda, &work[1], (ftnlen)1);
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
