#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
int zunmtr_(char *side, char *uplo, char *trans, integer *m, integer *n, doublecomplex *a,
            integer *lda, doublecomplex *tau, doublecomplex *c__, integer *ldc, doublecomplex *work,
            integer *lwork, integer *info, ftnlen side_len, ftnlen uplo_len, ftnlen trans_len)
{
    address a__1[2];
    integer a_dim1, a_offset, c_dim1, c_offset, i__1[2], i__2, i__3;
    char ch__1[2];
    int s_lmp_cat(char *, char **, integer *, integer *, ftnlen);
    integer i1, i2, nb, mi, ni, nq, nw;
    logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer iinfo;
    logical upper;
    extern int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    integer lwkopt;
    logical lquery;
    extern int zunmql_(char *, char *, integer *, integer *, integer *, doublecomplex *, integer *,
                       doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *,
                       integer *, ftnlen, ftnlen),
        zunmqr_(char *, char *, integer *, integer *, integer *, doublecomplex *, integer *,
                doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, integer *,
                ftnlen, ftnlen);
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
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    lquery = *lwork == -1;
    if (left) {
        nq = *m;
        nw = max(1, *n);
    } else {
        nq = *n;
        nw = max(1, *m);
    }
    if (!left && !lsame_(side, (char *)"R", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (!upper && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -2;
    } else if (!lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1) &&
               !lsame_(trans, (char *)"C", (ftnlen)1, (ftnlen)1)) {
        *info = -3;
    } else if (*m < 0) {
        *info = -4;
    } else if (*n < 0) {
        *info = -5;
    } else if (*lda < max(1, nq)) {
        *info = -7;
    } else if (*ldc < max(1, *m)) {
        *info = -10;
    } else if (*lwork < nw && !lquery) {
        *info = -12;
    }
    if (*info == 0) {
        if (upper) {
            if (left) {
                i__1[0] = 1, a__1[0] = side;
                i__1[1] = 1, a__1[1] = trans;
                s_lmp_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
                i__2 = *m - 1;
                i__3 = *m - 1;
                nb = ilaenv_(&c__1, (char *)"ZUNMQL", ch__1, &i__2, n, &i__3, &c_n1, (ftnlen)6, (ftnlen)2);
            } else {
                i__1[0] = 1, a__1[0] = side;
                i__1[1] = 1, a__1[1] = trans;
                s_lmp_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
                i__2 = *n - 1;
                i__3 = *n - 1;
                nb = ilaenv_(&c__1, (char *)"ZUNMQL", ch__1, m, &i__2, &i__3, &c_n1, (ftnlen)6, (ftnlen)2);
            }
        } else {
            if (left) {
                i__1[0] = 1, a__1[0] = side;
                i__1[1] = 1, a__1[1] = trans;
                s_lmp_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
                i__2 = *m - 1;
                i__3 = *m - 1;
                nb = ilaenv_(&c__1, (char *)"ZUNMQR", ch__1, &i__2, n, &i__3, &c_n1, (ftnlen)6, (ftnlen)2);
            } else {
                i__1[0] = 1, a__1[0] = side;
                i__1[1] = 1, a__1[1] = trans;
                s_lmp_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
                i__2 = *n - 1;
                i__3 = *n - 1;
                nb = ilaenv_(&c__1, (char *)"ZUNMQR", ch__1, m, &i__2, &i__3, &c_n1, (ftnlen)6, (ftnlen)2);
            }
        }
        lwkopt = nw * nb;
        work[1].r = (doublereal)lwkopt, work[1].i = 0.;
    }
    if (*info != 0) {
        i__2 = -(*info);
        xerbla_((char *)"ZUNMTR", &i__2, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    if (*m == 0 || *n == 0 || nq == 1) {
        work[1].r = 1., work[1].i = 0.;
        return 0;
    }
    if (left) {
        mi = *m - 1;
        ni = *n;
    } else {
        mi = *m;
        ni = *n - 1;
    }
    if (upper) {
        i__2 = nq - 1;
        zunmql_(side, trans, &mi, &ni, &i__2, &a[(a_dim1 << 1) + 1], lda, &tau[1], &c__[c_offset],
                ldc, &work[1], lwork, &iinfo, (ftnlen)1, (ftnlen)1);
    } else {
        if (left) {
            i1 = 2;
            i2 = 1;
        } else {
            i1 = 1;
            i2 = 2;
        }
        i__2 = nq - 1;
        zunmqr_(side, trans, &mi, &ni, &i__2, &a[a_dim1 + 2], lda, &tau[1], &c__[i1 + i2 * c_dim1],
                ldc, &work[1], lwork, &iinfo, (ftnlen)1, (ftnlen)1);
    }
    work[1].r = (doublereal)lwkopt, work[1].i = 0.;
    return 0;
}
#ifdef __cplusplus
}
#endif
