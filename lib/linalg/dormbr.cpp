#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
int dormbr_(char *vect, char *side, char *trans, integer *m, integer *n, integer *k, doublereal *a,
            integer *lda, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work,
            integer *lwork, integer *info, ftnlen vect_len, ftnlen side_len, ftnlen trans_len)
{
    address a__1[2];
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3[2];
    char ch__1[2];
    int s_lmp_cat(char *, char **, integer *, integer *, ftnlen);
    integer i1, i2, nb, mi, ni, nq, nw;
    logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer iinfo;
    extern int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    extern int dormlq_(char *, char *, integer *, integer *, integer *, doublereal *, integer *,
                       doublereal *, doublereal *, integer *, doublereal *, integer *, integer *,
                       ftnlen, ftnlen);
    logical notran;
    extern int dormqr_(char *, char *, integer *, integer *, integer *, doublereal *, integer *,
                       doublereal *, doublereal *, integer *, doublereal *, integer *, integer *,
                       ftnlen, ftnlen);
    logical applyq;
    char transt[1];
    integer lwkopt;
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
    applyq = lsame_(vect, (char *)"Q", (ftnlen)1, (ftnlen)1);
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
    if (!applyq && !lsame_(vect, (char *)"P", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (!left && !lsame_(side, (char *)"R", (ftnlen)1, (ftnlen)1)) {
        *info = -2;
    } else if (!notran && !lsame_(trans, (char *)"T", (ftnlen)1, (ftnlen)1)) {
        *info = -3;
    } else if (*m < 0) {
        *info = -4;
    } else if (*n < 0) {
        *info = -5;
    } else if (*k < 0) {
        *info = -6;
    } else {
        i__1 = 1, i__2 = min(nq, *k);
        if (applyq && *lda < max(1, nq) || !applyq && *lda < max(i__1, i__2)) {
            *info = -8;
        } else if (*ldc < max(1, *m)) {
            *info = -11;
        } else if (*lwork < nw && !lquery) {
            *info = -13;
        }
    }
    if (*info == 0) {
        if (applyq) {
            if (left) {
                i__3[0] = 1, a__1[0] = side;
                i__3[1] = 1, a__1[1] = trans;
                s_lmp_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
                i__1 = *m - 1;
                i__2 = *m - 1;
                nb = ilaenv_(&c__1, (char *)"DORMQR", ch__1, &i__1, n, &i__2, &c_n1, (ftnlen)6, (ftnlen)2);
            } else {
                i__3[0] = 1, a__1[0] = side;
                i__3[1] = 1, a__1[1] = trans;
                s_lmp_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
                i__1 = *n - 1;
                i__2 = *n - 1;
                nb = ilaenv_(&c__1, (char *)"DORMQR", ch__1, m, &i__1, &i__2, &c_n1, (ftnlen)6, (ftnlen)2);
            }
        } else {
            if (left) {
                i__3[0] = 1, a__1[0] = side;
                i__3[1] = 1, a__1[1] = trans;
                s_lmp_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
                i__1 = *m - 1;
                i__2 = *m - 1;
                nb = ilaenv_(&c__1, (char *)"DORMLQ", ch__1, &i__1, n, &i__2, &c_n1, (ftnlen)6, (ftnlen)2);
            } else {
                i__3[0] = 1, a__1[0] = side;
                i__3[1] = 1, a__1[1] = trans;
                s_lmp_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
                i__1 = *n - 1;
                i__2 = *n - 1;
                nb = ilaenv_(&c__1, (char *)"DORMLQ", ch__1, m, &i__1, &i__2, &c_n1, (ftnlen)6, (ftnlen)2);
            }
        }
        lwkopt = nw * nb;
        work[1] = (doublereal)lwkopt;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DORMBR", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    work[1] = 1.;
    if (*m == 0 || *n == 0) {
        return 0;
    }
    if (applyq) {
        if (nq >= *k) {
            dormqr_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[c_offset], ldc, &work[1],
                    lwork, &iinfo, (ftnlen)1, (ftnlen)1);
        } else if (nq > 1) {
            if (left) {
                mi = *m - 1;
                ni = *n;
                i1 = 2;
                i2 = 1;
            } else {
                mi = *m;
                ni = *n - 1;
                i1 = 1;
                i2 = 2;
            }
            i__1 = nq - 1;
            dormqr_(side, trans, &mi, &ni, &i__1, &a[a_dim1 + 2], lda, &tau[1],
                    &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &iinfo, (ftnlen)1, (ftnlen)1);
        }
    } else {
        if (notran) {
            *(unsigned char *)transt = 'T';
        } else {
            *(unsigned char *)transt = 'N';
        }
        if (nq > *k) {
            dormlq_(side, transt, m, n, k, &a[a_offset], lda, &tau[1], &c__[c_offset], ldc,
                    &work[1], lwork, &iinfo, (ftnlen)1, (ftnlen)1);
        } else if (nq > 1) {
            if (left) {
                mi = *m - 1;
                ni = *n;
                i1 = 2;
                i2 = 1;
            } else {
                mi = *m;
                ni = *n - 1;
                i1 = 1;
                i2 = 2;
            }
            i__1 = nq - 1;
            dormlq_(side, transt, &mi, &ni, &i__1, &a[(a_dim1 << 1) + 1], lda, &tau[1],
                    &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &iinfo, (ftnlen)1, (ftnlen)1);
        }
    }
    work[1] = (doublereal)lwkopt;
    return 0;
}
#ifdef __cplusplus
}
#endif
