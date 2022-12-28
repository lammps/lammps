#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
int zunm2r_(char *side, char *trans, integer *m, integer *n, integer *k, doublecomplex *a,
            integer *lda, doublecomplex *tau, doublecomplex *c__, integer *ldc, doublecomplex *work,
            integer *info, ftnlen side_len, ftnlen trans_len)
{
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3;
    doublecomplex z__1;
    void d_lmp_cnjg(doublecomplex *, doublecomplex *);
    integer i__, i1, i2, i3, ic, jc, mi, ni, nq;
    doublecomplex aii;
    logical left;
    doublecomplex taui;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int zlarf_(char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *,
                      doublecomplex *, integer *, doublecomplex *, ftnlen),
        xerbla_(char *, integer *, ftnlen);
    logical notran;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    *info = 0;
    left = lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1);
    notran = lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1);
    if (left) {
        nq = *m;
    } else {
        nq = *n;
    }
    if (!left && !lsame_(side, (char *)"R", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (!notran && !lsame_(trans, (char *)"C", (ftnlen)1, (ftnlen)1)) {
        *info = -2;
    } else if (*m < 0) {
        *info = -3;
    } else if (*n < 0) {
        *info = -4;
    } else if (*k < 0 || *k > nq) {
        *info = -5;
    } else if (*lda < max(1, nq)) {
        *info = -7;
    } else if (*ldc < max(1, *m)) {
        *info = -10;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"ZUNM2R", &i__1, (ftnlen)6);
        return 0;
    }
    if (*m == 0 || *n == 0 || *k == 0) {
        return 0;
    }
    if (left && !notran || !left && notran) {
        i1 = 1;
        i2 = *k;
        i3 = 1;
    } else {
        i1 = *k;
        i2 = 1;
        i3 = -1;
    }
    if (left) {
        ni = *n;
        jc = 1;
    } else {
        mi = *m;
        ic = 1;
    }
    i__1 = i2;
    i__2 = i3;
    for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
        if (left) {
            mi = *m - i__ + 1;
            ic = i__;
        } else {
            ni = *n - i__ + 1;
            jc = i__;
        }
        if (notran) {
            i__3 = i__;
            taui.r = tau[i__3].r, taui.i = tau[i__3].i;
        } else {
            d_lmp_cnjg(&z__1, &tau[i__]);
            taui.r = z__1.r, taui.i = z__1.i;
        }
        i__3 = i__ + i__ * a_dim1;
        aii.r = a[i__3].r, aii.i = a[i__3].i;
        i__3 = i__ + i__ * a_dim1;
        a[i__3].r = 1., a[i__3].i = 0.;
        zlarf_(side, &mi, &ni, &a[i__ + i__ * a_dim1], &c__1, &taui, &c__[ic + jc * c_dim1], ldc,
               &work[1], (ftnlen)1);
        i__3 = i__ + i__ * a_dim1;
        a[i__3].r = aii.r, a[i__3].i = aii.i;
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
