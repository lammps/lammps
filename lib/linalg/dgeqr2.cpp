#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
int dgeqr2_(integer *m, integer *n, doublereal *a, integer *lda, doublereal *tau, doublereal *work,
            integer *info)
{
    integer a_dim1, a_offset, i__1, i__2, i__3;
    integer i__, k;
    doublereal aii;
    extern int dlarf_(char *, integer *, integer *, doublereal *, integer *, doublereal *,
                      doublereal *, integer *, doublereal *, ftnlen),
        dlarfg_(integer *, doublereal *, doublereal *, integer *, doublereal *),
        xerbla_(char *, integer *, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
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
        xerbla_((char *)"DGEQR2", &i__1, (ftnlen)6);
        return 0;
    }
    k = min(*m, *n);
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
        i__2 = *m - i__ + 1;
        i__3 = i__ + 1;
        dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[min(i__3, *m) + i__ * a_dim1], &c__1, &tau[i__]);
        if (i__ < *n) {
            aii = a[i__ + i__ * a_dim1];
            a[i__ + i__ * a_dim1] = 1.;
            i__2 = *m - i__ + 1;
            i__3 = *n - i__;
            dlarf_((char *)"Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &tau[i__],
                   &a[i__ + (i__ + 1) * a_dim1], lda, &work[1], (ftnlen)4);
            a[i__ + i__ * a_dim1] = aii;
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
