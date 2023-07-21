#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
int dgelqf_(integer *m, integer *n, doublereal *a, integer *lda, doublereal *tau, doublereal *work,
            integer *lwork, integer *info)
{
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    integer i__, k, ib, nb, nx, iws, nbmin, iinfo;
    extern int dgelq2_(integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                       integer *),
        dlarfb_(char *, char *, char *, char *, integer *, integer *, integer *, doublereal *,
                integer *, doublereal *, integer *, doublereal *, integer *, doublereal *,
                integer *, ftnlen, ftnlen, ftnlen, ftnlen),
        dlarft_(char *, char *, integer *, integer *, doublereal *, integer *, doublereal *,
                doublereal *, integer *, ftnlen, ftnlen),
        xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    integer ldwork, lwkopt;
    logical lquery;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    *info = 0;
    nb = ilaenv_(&c__1, (char *)"DGELQF", (char *)" ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
    lwkopt = *m * nb;
    work[1] = (doublereal)lwkopt;
    lquery = *lwork == -1;
    if (*m < 0) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*lda < max(1, *m)) {
        *info = -4;
    } else if (*lwork < max(1, *m) && !lquery) {
        *info = -7;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DGELQF", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    k = min(*m, *n);
    if (k == 0) {
        work[1] = 1.;
        return 0;
    }
    nbmin = 2;
    nx = 0;
    iws = *m;
    if (nb > 1 && nb < k) {
        i__1 = 0, i__2 = ilaenv_(&c__3, (char *)"DGELQF", (char *)" ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
        nx = max(i__1, i__2);
        if (nx < k) {
            ldwork = *m;
            iws = ldwork * nb;
            if (*lwork < iws) {
                nb = *lwork / ldwork;
                i__1 = 2,
                i__2 = ilaenv_(&c__2, (char *)"DGELQF", (char *)" ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
                nbmin = max(i__1, i__2);
            }
        }
    }
    if (nb >= nbmin && nb < k && nx < k) {
        i__1 = k - nx;
        i__2 = nb;
        for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
            i__3 = k - i__ + 1;
            ib = min(i__3, nb);
            i__3 = *n - i__ + 1;
            dgelq2_(&ib, &i__3, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1], &iinfo);
            if (i__ + ib <= *m) {
                i__3 = *n - i__ + 1;
                dlarft_((char *)"Forward", (char *)"Rowwise", &i__3, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__],
                        &work[1], &ldwork, (ftnlen)7, (ftnlen)7);
                i__3 = *m - i__ - ib + 1;
                i__4 = *n - i__ + 1;
                dlarfb_((char *)"Right", (char *)"No transpose", (char *)"Forward", (char *)"Rowwise", &i__3, &i__4, &ib,
                        &a[i__ + i__ * a_dim1], lda, &work[1], &ldwork, &a[i__ + ib + i__ * a_dim1],
                        lda, &work[ib + 1], &ldwork, (ftnlen)5, (ftnlen)12, (ftnlen)7, (ftnlen)7);
            }
        }
    } else {
        i__ = 1;
    }
    if (i__ <= k) {
        i__2 = *m - i__ + 1;
        i__1 = *n - i__ + 1;
        dgelq2_(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1], &iinfo);
    }
    work[1] = (doublereal)iws;
    return 0;
}
#ifdef __cplusplus
}
#endif
