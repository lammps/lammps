#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static doublereal c_b16 = -1.;
int zpptrf_(char *uplo, integer *n, doublecomplex *ap, integer *info, ftnlen uplo_len)
{
    integer i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1;
    double sqrt(doublereal);
    integer j, jc, jj;
    doublereal ajj;
    extern int zhpr_(char *, integer *, doublereal *, doublecomplex *, integer *, doublecomplex *,
                     ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern VOID zdotc_(doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *,
                       integer *);
    logical upper;
    extern int ztpsv_(char *, char *, char *, integer *, doublecomplex *, doublecomplex *,
                      integer *, ftnlen, ftnlen, ftnlen),
        xerbla_(char *, integer *, ftnlen),
        zdscal_(integer *, doublereal *, doublecomplex *, integer *);
    --ap;
    *info = 0;
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    if (!upper && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"ZPPTRF", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    if (upper) {
        jj = 0;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            jc = jj + 1;
            jj += j;
            if (j > 1) {
                i__2 = j - 1;
                ztpsv_((char *)"Upper", (char *)"Conjugate transpose", (char *)"Non-unit", &i__2, &ap[1], &ap[jc], &c__1,
                       (ftnlen)5, (ftnlen)19, (ftnlen)8);
            }
            i__2 = jj;
            i__3 = j - 1;
            zdotc_(&z__1, &i__3, &ap[jc], &c__1, &ap[jc], &c__1);
            ajj = ap[i__2].r - z__1.r;
            if (ajj <= 0.) {
                i__2 = jj;
                ap[i__2].r = ajj, ap[i__2].i = 0.;
                goto L30;
            }
            i__2 = jj;
            d__1 = sqrt(ajj);
            ap[i__2].r = d__1, ap[i__2].i = 0.;
        }
    } else {
        jj = 1;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__2 = jj;
            ajj = ap[i__2].r;
            if (ajj <= 0.) {
                i__2 = jj;
                ap[i__2].r = ajj, ap[i__2].i = 0.;
                goto L30;
            }
            ajj = sqrt(ajj);
            i__2 = jj;
            ap[i__2].r = ajj, ap[i__2].i = 0.;
            if (j < *n) {
                i__2 = *n - j;
                d__1 = 1. / ajj;
                zdscal_(&i__2, &d__1, &ap[jj + 1], &c__1);
                i__2 = *n - j;
                zhpr_((char *)"Lower", &i__2, &c_b16, &ap[jj + 1], &c__1, &ap[jj + *n - j + 1], (ftnlen)5);
                jj = jj + *n - j + 1;
            }
        }
    }
    goto L40;
L30:
    *info = j;
L40:
    return 0;
}
#ifdef __cplusplus
}
#endif
