#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__65 = 65;
int dormqr_(char *side, char *trans, integer *m, integer *n, integer *k, doublereal *a,
            integer *lda, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work,
            integer *lwork, integer *info, ftnlen side_len, ftnlen trans_len)
{
    address a__1[2];
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3[2], i__4, i__5;
    char ch__1[2];
    int s_lmp_cat(char *, char **, integer *, integer *, ftnlen);
    integer i__, i1, i2, i3, ib, ic, jc, nb, mi, ni, nq, nw, iwt;
    logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer nbmin, iinfo;
    extern int dorm2r_(char *, char *, integer *, integer *, integer *, doublereal *, integer *,
                       doublereal *, doublereal *, integer *, doublereal *, integer *, ftnlen,
                       ftnlen),
        dlarfb_(char *, char *, char *, char *, integer *, integer *, integer *, doublereal *,
                integer *, doublereal *, integer *, doublereal *, integer *, doublereal *,
                integer *, ftnlen, ftnlen, ftnlen, ftnlen),
        dlarft_(char *, char *, integer *, integer *, doublereal *, integer *, doublereal *,
                doublereal *, integer *, ftnlen, ftnlen),
        xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    logical notran;
    integer ldwork, lwkopt;
    logical lquery;
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
    } else if (*lwork < nw && !lquery) {
        *info = -12;
    }
    if (*info == 0) {
        i__3[0] = 1, a__1[0] = side;
        i__3[1] = 1, a__1[1] = trans;
        s_lmp_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
        i__1 = 64, i__2 = ilaenv_(&c__1, (char *)"DORMQR", ch__1, m, n, k, &c_n1, (ftnlen)6, (ftnlen)2);
        nb = min(i__1, i__2);
        lwkopt = nw * nb + 4160;
        work[1] = (doublereal)lwkopt;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DORMQR", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    if (*m == 0 || *n == 0 || *k == 0) {
        work[1] = 1.;
        return 0;
    }
    nbmin = 2;
    ldwork = nw;
    if (nb > 1 && nb < *k) {
        if (*lwork < lwkopt) {
            nb = (*lwork - 4160) / ldwork;
            i__3[0] = 1, a__1[0] = side;
            i__3[1] = 1, a__1[1] = trans;
            s_lmp_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
            i__1 = 2, i__2 = ilaenv_(&c__2, (char *)"DORMQR", ch__1, m, n, k, &c_n1, (ftnlen)6, (ftnlen)2);
            nbmin = max(i__1, i__2);
        }
    }
    if (nb < nbmin || nb >= *k) {
        dorm2r_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[c_offset], ldc, &work[1],
                &iinfo, (ftnlen)1, (ftnlen)1);
    } else {
        iwt = nw * nb + 1;
        if (left && !notran || !left && notran) {
            i1 = 1;
            i2 = *k;
            i3 = nb;
        } else {
            i1 = (*k - 1) / nb * nb + 1;
            i2 = 1;
            i3 = -nb;
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
            i__4 = nb, i__5 = *k - i__ + 1;
            ib = min(i__4, i__5);
            i__4 = nq - i__ + 1;
            dlarft_((char *)"Forward", (char *)"Columnwise", &i__4, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__],
                    &work[iwt], &c__65, (ftnlen)7, (ftnlen)10);
            if (left) {
                mi = *m - i__ + 1;
                ic = i__;
            } else {
                ni = *n - i__ + 1;
                jc = i__;
            }
            dlarfb_(side, trans, (char *)"Forward", (char *)"Columnwise", &mi, &ni, &ib, &a[i__ + i__ * a_dim1],
                    lda, &work[iwt], &c__65, &c__[ic + jc * c_dim1], ldc, &work[1], &ldwork,
                    (ftnlen)1, (ftnlen)1, (ftnlen)7, (ftnlen)10);
        }
    }
    work[1] = (doublereal)lwkopt;
    return 0;
}
#ifdef __cplusplus
}
#endif
