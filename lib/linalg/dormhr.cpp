#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
int dormhr_(char *side, char *trans, integer *m, integer *n, integer *ilo, integer *ihi,
            doublereal *a, integer *lda, doublereal *tau, doublereal *c__, integer *ldc,
            doublereal *work, integer *lwork, integer *info, ftnlen side_len, ftnlen trans_len)
{
    address a__1[2];
    integer a_dim1, a_offset, c_dim1, c_offset, i__1[2], i__2;
    char ch__1[2];
    int s_lmp_cat(char *, char **, integer *, integer *, ftnlen);
    integer i1, i2, nb, mi, nh, ni, nq, nw;
    logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer iinfo;
    extern int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    extern int dormqr_(char *, char *, integer *, integer *, integer *, doublereal *, integer *,
                       doublereal *, doublereal *, integer *, doublereal *, integer *, integer *,
                       ftnlen, ftnlen);
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
    nh = *ihi - *ilo;
    left = lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1);
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
    } else if (!lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1) &&
               !lsame_(trans, (char *)"T", (ftnlen)1, (ftnlen)1)) {
        *info = -2;
    } else if (*m < 0) {
        *info = -3;
    } else if (*n < 0) {
        *info = -4;
    } else if (*ilo < 1 || *ilo > max(1, nq)) {
        *info = -5;
    } else if (*ihi < min(*ilo, nq) || *ihi > nq) {
        *info = -6;
    } else if (*lda < max(1, nq)) {
        *info = -8;
    } else if (*ldc < max(1, *m)) {
        *info = -11;
    } else if (*lwork < nw && !lquery) {
        *info = -13;
    }
    if (*info == 0) {
        if (left) {
            i__1[0] = 1, a__1[0] = side;
            i__1[1] = 1, a__1[1] = trans;
            s_lmp_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
            nb = ilaenv_(&c__1, (char *)"DORMQR", ch__1, &nh, n, &nh, &c_n1, (ftnlen)6, (ftnlen)2);
        } else {
            i__1[0] = 1, a__1[0] = side;
            i__1[1] = 1, a__1[1] = trans;
            s_lmp_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
            nb = ilaenv_(&c__1, (char *)"DORMQR", ch__1, m, &nh, &nh, &c_n1, (ftnlen)6, (ftnlen)2);
        }
        lwkopt = nw * nb;
        work[1] = (doublereal)lwkopt;
    }
    if (*info != 0) {
        i__2 = -(*info);
        xerbla_((char *)"DORMHR", &i__2, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    if (*m == 0 || *n == 0 || nh == 0) {
        work[1] = 1.;
        return 0;
    }
    if (left) {
        mi = nh;
        ni = *n;
        i1 = *ilo + 1;
        i2 = 1;
    } else {
        mi = *m;
        ni = nh;
        i1 = 1;
        i2 = *ilo + 1;
    }
    dormqr_(side, trans, &mi, &ni, &nh, &a[*ilo + 1 + *ilo * a_dim1], lda, &tau[*ilo],
            &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &iinfo, (ftnlen)1, (ftnlen)1);
    work[1] = (doublereal)lwkopt;
    return 0;
}
#ifdef __cplusplus
}
#endif
