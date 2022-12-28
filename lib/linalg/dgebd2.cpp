#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
int dgebd2_(integer *m, integer *n, doublereal *a, integer *lda, doublereal *d__, doublereal *e,
            doublereal *tauq, doublereal *taup, doublereal *work, integer *info)
{
    integer a_dim1, a_offset, i__1, i__2, i__3;
    integer i__;
    extern int dlarf_(char *, integer *, integer *, doublereal *, integer *, doublereal *,
                      doublereal *, integer *, doublereal *, ftnlen),
        dlarfg_(integer *, doublereal *, doublereal *, integer *, doublereal *),
        xerbla_(char *, integer *, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tauq;
    --taup;
    --work;
    *info = 0;
    if (*m < 0) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*lda < max(1, *m)) {
        *info = -4;
    }
    if (*info < 0) {
        i__1 = -(*info);
        xerbla_((char *)"DGEBD2", &i__1, (ftnlen)6);
        return 0;
    }
    if (*m >= *n) {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            i__2 = *m - i__ + 1;
            i__3 = i__ + 1;
            dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[min(i__3, *m) + i__ * a_dim1], &c__1,
                    &tauq[i__]);
            d__[i__] = a[i__ + i__ * a_dim1];
            a[i__ + i__ * a_dim1] = 1.;
            if (i__ < *n) {
                i__2 = *m - i__ + 1;
                i__3 = *n - i__;
                dlarf_((char *)"Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &tauq[i__],
                       &a[i__ + (i__ + 1) * a_dim1], lda, &work[1], (ftnlen)4);
            }
            a[i__ + i__ * a_dim1] = d__[i__];
            if (i__ < *n) {
                i__2 = *n - i__;
                i__3 = i__ + 2;
                dlarfg_(&i__2, &a[i__ + (i__ + 1) * a_dim1], &a[i__ + min(i__3, *n) * a_dim1], lda,
                        &taup[i__]);
                e[i__] = a[i__ + (i__ + 1) * a_dim1];
                a[i__ + (i__ + 1) * a_dim1] = 1.;
                i__2 = *m - i__;
                i__3 = *n - i__;
                dlarf_((char *)"Right", &i__2, &i__3, &a[i__ + (i__ + 1) * a_dim1], lda, &taup[i__],
                       &a[i__ + 1 + (i__ + 1) * a_dim1], lda, &work[1], (ftnlen)5);
                a[i__ + (i__ + 1) * a_dim1] = e[i__];
            } else {
                taup[i__] = 0.;
            }
        }
    } else {
        i__1 = *m;
        for (i__ = 1; i__ <= i__1; ++i__) {
            i__2 = *n - i__ + 1;
            i__3 = i__ + 1;
            dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[i__ + min(i__3, *n) * a_dim1], lda,
                    &taup[i__]);
            d__[i__] = a[i__ + i__ * a_dim1];
            a[i__ + i__ * a_dim1] = 1.;
            if (i__ < *m) {
                i__2 = *m - i__;
                i__3 = *n - i__ + 1;
                dlarf_((char *)"Right", &i__2, &i__3, &a[i__ + i__ * a_dim1], lda, &taup[i__],
                       &a[i__ + 1 + i__ * a_dim1], lda, &work[1], (ftnlen)5);
            }
            a[i__ + i__ * a_dim1] = d__[i__];
            if (i__ < *m) {
                i__2 = *m - i__;
                i__3 = i__ + 2;
                dlarfg_(&i__2, &a[i__ + 1 + i__ * a_dim1], &a[min(i__3, *m) + i__ * a_dim1], &c__1,
                        &tauq[i__]);
                e[i__] = a[i__ + 1 + i__ * a_dim1];
                a[i__ + 1 + i__ * a_dim1] = 1.;
                i__2 = *m - i__;
                i__3 = *n - i__;
                dlarf_((char *)"Left", &i__2, &i__3, &a[i__ + 1 + i__ * a_dim1], &c__1, &tauq[i__],
                       &a[i__ + 1 + (i__ + 1) * a_dim1], lda, &work[1], (ftnlen)4);
                a[i__ + 1 + i__ * a_dim1] = e[i__];
            } else {
                tauq[i__] = 0.;
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
