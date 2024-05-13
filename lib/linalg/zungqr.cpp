#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
int zungqr_(integer *m, integer *n, integer *k, doublecomplex *a, integer *lda, doublecomplex *tau,
            doublecomplex *work, integer *lwork, integer *info)
{
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    integer i__, j, l, ib, nb, ki, kk, nx, iws, nbmin, iinfo;
    extern int zung2r_(integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *,
                       doublecomplex *, integer *),
        xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    extern int zlarfb_(char *, char *, char *, char *, integer *, integer *, integer *,
                       doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *,
                       integer *, doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    integer ldwork;
    extern int zlarft_(char *, char *, integer *, integer *, doublecomplex *, integer *,
                       doublecomplex *, doublecomplex *, integer *, ftnlen, ftnlen);
    integer lwkopt;
    logical lquery;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    *info = 0;
    nb = ilaenv_(&c__1, (char *)"ZUNGQR", (char *)" ", m, n, k, &c_n1, (ftnlen)6, (ftnlen)1);
    lwkopt = max(1, *n) * nb;
    work[1].r = (doublereal)lwkopt, work[1].i = 0.;
    lquery = *lwork == -1;
    if (*m < 0) {
        *info = -1;
    } else if (*n < 0 || *n > *m) {
        *info = -2;
    } else if (*k < 0 || *k > *n) {
        *info = -3;
    } else if (*lda < max(1, *m)) {
        *info = -5;
    } else if (*lwork < max(1, *n) && !lquery) {
        *info = -8;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"ZUNGQR", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    if (*n <= 0) {
        work[1].r = 1., work[1].i = 0.;
        return 0;
    }
    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < *k) {
        i__1 = 0, i__2 = ilaenv_(&c__3, (char *)"ZUNGQR", (char *)" ", m, n, k, &c_n1, (ftnlen)6, (ftnlen)1);
        nx = max(i__1, i__2);
        if (nx < *k) {
            ldwork = *n;
            iws = ldwork * nb;
            if (*lwork < iws) {
                nb = *lwork / ldwork;
                i__1 = 2,
                i__2 = ilaenv_(&c__2, (char *)"ZUNGQR", (char *)" ", m, n, k, &c_n1, (ftnlen)6, (ftnlen)1);
                nbmin = max(i__1, i__2);
            }
        }
    }
    if (nb >= nbmin && nb < *k && nx < *k) {
        ki = (*k - nx - 1) / nb * nb;
        i__1 = *k, i__2 = ki + nb;
        kk = min(i__1, i__2);
        i__1 = *n;
        for (j = kk + 1; j <= i__1; ++j) {
            i__2 = kk;
            for (i__ = 1; i__ <= i__2; ++i__) {
                i__3 = i__ + j * a_dim1;
                a[i__3].r = 0., a[i__3].i = 0.;
            }
        }
    } else {
        kk = 0;
    }
    if (kk < *n) {
        i__1 = *m - kk;
        i__2 = *n - kk;
        i__3 = *k - kk;
        zung2r_(&i__1, &i__2, &i__3, &a[kk + 1 + (kk + 1) * a_dim1], lda, &tau[kk + 1], &work[1],
                &iinfo);
    }
    if (kk > 0) {
        i__1 = -nb;
        for (i__ = ki + 1; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
            i__2 = nb, i__3 = *k - i__ + 1;
            ib = min(i__2, i__3);
            if (i__ + ib <= *n) {
                i__2 = *m - i__ + 1;
                zlarft_((char *)"Forward", (char *)"Columnwise", &i__2, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__],
                        &work[1], &ldwork, (ftnlen)7, (ftnlen)10);
                i__2 = *m - i__ + 1;
                i__3 = *n - i__ - ib + 1;
                zlarfb_((char *)"Left", (char *)"No transpose", (char *)"Forward", (char *)"Columnwise", &i__2, &i__3, &ib,
                        &a[i__ + i__ * a_dim1], lda, &work[1], &ldwork,
                        &a[i__ + (i__ + ib) * a_dim1], lda, &work[ib + 1], &ldwork, (ftnlen)4,
                        (ftnlen)12, (ftnlen)7, (ftnlen)10);
            }
            i__2 = *m - i__ + 1;
            zung2r_(&i__2, &ib, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1], &iinfo);
            i__2 = i__ + ib - 1;
            for (j = i__; j <= i__2; ++j) {
                i__3 = i__ - 1;
                for (l = 1; l <= i__3; ++l) {
                    i__4 = l + j * a_dim1;
                    a[i__4].r = 0., a[i__4].i = 0.;
                }
            }
        }
    }
    work[1].r = (doublereal)iws, work[1].i = 0.;
    return 0;
}
#ifdef __cplusplus
}
#endif
