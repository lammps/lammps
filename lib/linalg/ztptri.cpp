#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublecomplex c_b1 = {1., 0.};
static integer c__1 = 1;
int ztptri_(char *uplo, char *diag, integer *n, doublecomplex *ap, integer *info, ftnlen uplo_len,
            ftnlen diag_len)
{
    integer i__1, i__2;
    doublecomplex z__1;
    void z_lmp_div(doublecomplex *, doublecomplex *, doublecomplex *);
    integer j, jc, jj;
    doublecomplex ajj;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int zscal_(integer *, doublecomplex *, doublecomplex *, integer *);
    logical upper;
    extern int ztpmv_(char *, char *, char *, integer *, doublecomplex *, doublecomplex *,
                      integer *, ftnlen, ftnlen, ftnlen),
        xerbla_(char *, integer *, ftnlen);
    integer jclast;
    logical nounit;
    --ap;
    *info = 0;
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    nounit = lsame_(diag, (char *)"N", (ftnlen)1, (ftnlen)1);
    if (!upper && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (!nounit && !lsame_(diag, (char *)"U", (ftnlen)1, (ftnlen)1)) {
        *info = -2;
    } else if (*n < 0) {
        *info = -3;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"ZTPTRI", &i__1, (ftnlen)6);
        return 0;
    }
    if (nounit) {
        if (upper) {
            jj = 0;
            i__1 = *n;
            for (*info = 1; *info <= i__1; ++(*info)) {
                jj += *info;
                i__2 = jj;
                if (ap[i__2].r == 0. && ap[i__2].i == 0.) {
                    return 0;
                }
            }
        } else {
            jj = 1;
            i__1 = *n;
            for (*info = 1; *info <= i__1; ++(*info)) {
                i__2 = jj;
                if (ap[i__2].r == 0. && ap[i__2].i == 0.) {
                    return 0;
                }
                jj = jj + *n - *info + 1;
            }
        }
        *info = 0;
    }
    if (upper) {
        jc = 1;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            if (nounit) {
                i__2 = jc + j - 1;
                z_lmp_div(&z__1, &c_b1, &ap[jc + j - 1]);
                ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
                i__2 = jc + j - 1;
                z__1.r = -ap[i__2].r, z__1.i = -ap[i__2].i;
                ajj.r = z__1.r, ajj.i = z__1.i;
            } else {
                z__1.r = -1., z__1.i = -0.;
                ajj.r = z__1.r, ajj.i = z__1.i;
            }
            i__2 = j - 1;
            ztpmv_((char *)"Upper", (char *)"No transpose", diag, &i__2, &ap[1], &ap[jc], &c__1, (ftnlen)5,
                   (ftnlen)12, (ftnlen)1);
            i__2 = j - 1;
            zscal_(&i__2, &ajj, &ap[jc], &c__1);
            jc += j;
        }
    } else {
        jc = *n * (*n + 1) / 2;
        for (j = *n; j >= 1; --j) {
            if (nounit) {
                i__1 = jc;
                z_lmp_div(&z__1, &c_b1, &ap[jc]);
                ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
                i__1 = jc;
                z__1.r = -ap[i__1].r, z__1.i = -ap[i__1].i;
                ajj.r = z__1.r, ajj.i = z__1.i;
            } else {
                z__1.r = -1., z__1.i = -0.;
                ajj.r = z__1.r, ajj.i = z__1.i;
            }
            if (j < *n) {
                i__1 = *n - j;
                ztpmv_((char *)"Lower", (char *)"No transpose", diag, &i__1, &ap[jclast], &ap[jc + 1], &c__1,
                       (ftnlen)5, (ftnlen)12, (ftnlen)1);
                i__1 = *n - j;
                zscal_(&i__1, &ajj, &ap[jc + 1], &c__1);
            }
            jclast = jc;
            jc = jc - *n + j - 2;
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
