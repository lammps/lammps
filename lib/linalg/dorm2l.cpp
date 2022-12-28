#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
int dorm2l_(char *side, char *trans, integer *m, integer *n, integer *k, doublereal *a,
            integer *lda, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work,
            integer *info, ftnlen side_len, ftnlen trans_len)
{
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;
    integer i__, i1, i2, i3, mi, ni, nq;
    doublereal aii;
    logical left;
    extern int dlarf_(char *, integer *, integer *, doublereal *, integer *, doublereal *,
                      doublereal *, integer *, doublereal *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int xerbla_(char *, integer *, ftnlen);
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
    } else if (!notran && !lsame_(trans, (char *)"T", (ftnlen)1, (ftnlen)1)) {
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
        xerbla_((char *)"DORM2L", &i__1, (ftnlen)6);
        return 0;
    }
    if (*m == 0 || *n == 0 || *k == 0) {
        return 0;
    }
    if (left && notran || !left && !notran) {
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
    } else {
        mi = *m;
    }
    i__1 = i2;
    i__2 = i3;
    for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
        if (left) {
            mi = *m - *k + i__;
        } else {
            ni = *n - *k + i__;
        }
        aii = a[nq - *k + i__ + i__ * a_dim1];
        a[nq - *k + i__ + i__ * a_dim1] = 1.;
        dlarf_(side, &mi, &ni, &a[i__ * a_dim1 + 1], &c__1, &tau[i__], &c__[c_offset], ldc,
               &work[1], (ftnlen)1);
        a[nq - *k + i__ + i__ * a_dim1] = aii;
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
