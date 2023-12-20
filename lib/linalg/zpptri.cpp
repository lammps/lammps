#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b8 = 1.;
static integer c__1 = 1;
int zpptri_(char *uplo, integer *n, doublecomplex *ap, integer *info, ftnlen uplo_len)
{
    integer i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1;
    integer j, jc, jj;
    doublereal ajj;
    integer jjn;
    extern int zhpr_(char *, integer *, doublereal *, doublecomplex *, integer *, doublecomplex *,
                     ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern VOID zdotc_(doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *,
                       integer *);
    logical upper;
    extern int ztpmv_(char *, char *, char *, integer *, doublecomplex *, doublecomplex *,
                      integer *, ftnlen, ftnlen, ftnlen),
        xerbla_(char *, integer *, ftnlen),
        zdscal_(integer *, doublereal *, doublecomplex *, integer *),
        ztptri_(char *, char *, integer *, doublecomplex *, integer *, ftnlen, ftnlen);
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
        xerbla_((char *)"ZPPTRI", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    ztptri_(uplo, (char *)"Non-unit", n, &ap[1], info, (ftnlen)1, (ftnlen)8);
    if (*info > 0) {
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
                zhpr_((char *)"Upper", &i__2, &c_b8, &ap[jc], &c__1, &ap[1], (ftnlen)5);
            }
            i__2 = jj;
            ajj = ap[i__2].r;
            zdscal_(&j, &ajj, &ap[jc], &c__1);
        }
    } else {
        jj = 1;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            jjn = jj + *n - j + 1;
            i__2 = jj;
            i__3 = *n - j + 1;
            zdotc_(&z__1, &i__3, &ap[jj], &c__1, &ap[jj], &c__1);
            d__1 = z__1.r;
            ap[i__2].r = d__1, ap[i__2].i = 0.;
            if (j < *n) {
                i__2 = *n - j;
                ztpmv_((char *)"Lower", (char *)"Conjugate transpose", (char *)"Non-unit", &i__2, &ap[jjn], &ap[jj + 1],
                       &c__1, (ftnlen)5, (ftnlen)19, (ftnlen)8);
            }
            jj = jjn;
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
