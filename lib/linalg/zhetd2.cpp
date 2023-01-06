#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublecomplex c_b2 = {0., 0.};
static integer c__1 = 1;
int zhetd2_(char *uplo, integer *n, doublecomplex *a, integer *lda, doublereal *d__, doublereal *e,
            doublecomplex *tau, integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4;
    integer i__;
    doublecomplex taui;
    extern int zher2_(char *, integer *, doublecomplex *, doublecomplex *, integer *,
                      doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    doublecomplex alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern VOID zdotc_(doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *,
                       integer *);
    extern int zhemv_(char *, integer *, doublecomplex *, doublecomplex *, integer *,
                      doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *,
                      ftnlen);
    logical upper;
    extern int zaxpy_(integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *,
                      integer *),
        xerbla_(char *, integer *, ftnlen),
        zlarfg_(integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tau;
    *info = 0;
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    if (!upper && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*lda < max(1, *n)) {
        *info = -4;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"ZHETD2", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n <= 0) {
        return 0;
    }
    if (upper) {
        i__1 = *n + *n * a_dim1;
        i__2 = *n + *n * a_dim1;
        d__1 = a[i__2].r;
        a[i__1].r = d__1, a[i__1].i = 0.;
        for (i__ = *n - 1; i__ >= 1; --i__) {
            i__1 = i__ + (i__ + 1) * a_dim1;
            alpha.r = a[i__1].r, alpha.i = a[i__1].i;
            zlarfg_(&i__, &alpha, &a[(i__ + 1) * a_dim1 + 1], &c__1, &taui);
            e[i__] = alpha.r;
            if (taui.r != 0. || taui.i != 0.) {
                i__1 = i__ + (i__ + 1) * a_dim1;
                a[i__1].r = 1., a[i__1].i = 0.;
                zhemv_(uplo, &i__, &taui, &a[a_offset], lda, &a[(i__ + 1) * a_dim1 + 1], &c__1,
                       &c_b2, &tau[1], &c__1, (ftnlen)1);
                z__3.r = -.5, z__3.i = -0.;
                z__2.r = z__3.r * taui.r - z__3.i * taui.i,
                z__2.i = z__3.r * taui.i + z__3.i * taui.r;
                zdotc_(&z__4, &i__, &tau[1], &c__1, &a[(i__ + 1) * a_dim1 + 1], &c__1);
                z__1.r = z__2.r * z__4.r - z__2.i * z__4.i,
                z__1.i = z__2.r * z__4.i + z__2.i * z__4.r;
                alpha.r = z__1.r, alpha.i = z__1.i;
                zaxpy_(&i__, &alpha, &a[(i__ + 1) * a_dim1 + 1], &c__1, &tau[1], &c__1);
                z__1.r = -1., z__1.i = -0.;
                zher2_(uplo, &i__, &z__1, &a[(i__ + 1) * a_dim1 + 1], &c__1, &tau[1], &c__1,
                       &a[a_offset], lda, (ftnlen)1);
            } else {
                i__1 = i__ + i__ * a_dim1;
                i__2 = i__ + i__ * a_dim1;
                d__1 = a[i__2].r;
                a[i__1].r = d__1, a[i__1].i = 0.;
            }
            i__1 = i__ + (i__ + 1) * a_dim1;
            i__2 = i__;
            a[i__1].r = e[i__2], a[i__1].i = 0.;
            i__1 = i__ + 1 + (i__ + 1) * a_dim1;
            d__[i__ + 1] = a[i__1].r;
            i__1 = i__;
            tau[i__1].r = taui.r, tau[i__1].i = taui.i;
        }
        i__1 = a_dim1 + 1;
        d__[1] = a[i__1].r;
    } else {
        i__1 = a_dim1 + 1;
        i__2 = a_dim1 + 1;
        d__1 = a[i__2].r;
        a[i__1].r = d__1, a[i__1].i = 0.;
        i__1 = *n - 1;
        for (i__ = 1; i__ <= i__1; ++i__) {
            i__2 = i__ + 1 + i__ * a_dim1;
            alpha.r = a[i__2].r, alpha.i = a[i__2].i;
            i__2 = *n - i__;
            i__3 = i__ + 2;
            zlarfg_(&i__2, &alpha, &a[min(i__3, *n) + i__ * a_dim1], &c__1, &taui);
            e[i__] = alpha.r;
            if (taui.r != 0. || taui.i != 0.) {
                i__2 = i__ + 1 + i__ * a_dim1;
                a[i__2].r = 1., a[i__2].i = 0.;
                i__2 = *n - i__;
                zhemv_(uplo, &i__2, &taui, &a[i__ + 1 + (i__ + 1) * a_dim1], lda,
                       &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b2, &tau[i__], &c__1, (ftnlen)1);
                z__3.r = -.5, z__3.i = -0.;
                z__2.r = z__3.r * taui.r - z__3.i * taui.i,
                z__2.i = z__3.r * taui.i + z__3.i * taui.r;
                i__2 = *n - i__;
                zdotc_(&z__4, &i__2, &tau[i__], &c__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
                z__1.r = z__2.r * z__4.r - z__2.i * z__4.i,
                z__1.i = z__2.r * z__4.i + z__2.i * z__4.r;
                alpha.r = z__1.r, alpha.i = z__1.i;
                i__2 = *n - i__;
                zaxpy_(&i__2, &alpha, &a[i__ + 1 + i__ * a_dim1], &c__1, &tau[i__], &c__1);
                i__2 = *n - i__;
                z__1.r = -1., z__1.i = -0.;
                zher2_(uplo, &i__2, &z__1, &a[i__ + 1 + i__ * a_dim1], &c__1, &tau[i__], &c__1,
                       &a[i__ + 1 + (i__ + 1) * a_dim1], lda, (ftnlen)1);
            } else {
                i__2 = i__ + 1 + (i__ + 1) * a_dim1;
                i__3 = i__ + 1 + (i__ + 1) * a_dim1;
                d__1 = a[i__3].r;
                a[i__2].r = d__1, a[i__2].i = 0.;
            }
            i__2 = i__ + 1 + i__ * a_dim1;
            i__3 = i__;
            a[i__2].r = e[i__3], a[i__2].i = 0.;
            i__2 = i__ + i__ * a_dim1;
            d__[i__] = a[i__2].r;
            i__2 = i__;
            tau[i__2].r = taui.r, tau[i__2].i = taui.i;
        }
        i__1 = *n + *n * a_dim1;
        d__[*n] = a[i__1].r;
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
