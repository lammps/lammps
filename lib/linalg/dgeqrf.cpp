#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
int dgeqrf_(integer *m, integer *n, doublereal *a, integer *lda, doublereal *tau, doublereal *work,
            integer *lwork, integer *info)
{
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    integer i__, k, ib, nb, nx, iws, nbmin, iinfo;
    extern int dgeqr2_(integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
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
    k = min(*m, *n);
    *info = 0;
    nb = ilaenv_(&c__1, (char *)"DGEQRF", (char *)" ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
    lquery = *lwork == -1;
    if (*m < 0) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*lda < max(1, *m)) {
        *info = -4;
    } else if (!lquery) {
        if (*lwork <= 0 || *m > 0 && *lwork < max(1, *n)) {
            *info = -7;
        }
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DGEQRF", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        if (k == 0) {
            lwkopt = 1;
        } else {
            lwkopt = *n * nb;
        }
        work[1] = (doublereal)lwkopt;
        return 0;
    }
    if (k == 0) {
        work[1] = 1.;
        return 0;
    }
    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < k) {
        i__1 = 0, i__2 = ilaenv_(&c__3, (char *)"DGEQRF", (char *)" ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
        nx = max(i__1, i__2);
        if (nx < k) {
            ldwork = *n;
            iws = ldwork * nb;
            if (*lwork < iws) {
                nb = *lwork / ldwork;
                i__1 = 2,
                i__2 = ilaenv_(&c__2, (char *)"DGEQRF", (char *)" ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
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
            i__3 = *m - i__ + 1;
            dgeqr2_(&i__3, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1], &iinfo);
            if (i__ + ib <= *n) {
                i__3 = *m - i__ + 1;
                dlarft_((char *)"Forward", (char *)"Columnwise", &i__3, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__],
                        &work[1], &ldwork, (ftnlen)7, (ftnlen)10);
                i__3 = *m - i__ + 1;
                i__4 = *n - i__ - ib + 1;
                dlarfb_((char *)"Left", (char *)"Transpose", (char *)"Forward", (char *)"Columnwise", &i__3, &i__4, &ib,
                        &a[i__ + i__ * a_dim1], lda, &work[1], &ldwork,
                        &a[i__ + (i__ + ib) * a_dim1], lda, &work[ib + 1], &ldwork, (ftnlen)4,
                        (ftnlen)9, (ftnlen)7, (ftnlen)10);
            }
        }
    } else {
        i__ = 1;
    }
    if (i__ <= k) {
        i__2 = *m - i__ + 1;
        i__1 = *n - i__ + 1;
        dgeqr2_(&i__2, &i__1, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1], &iinfo);
    }
    work[1] = (doublereal)iws;
    return 0;
}
#ifdef __cplusplus
}
#endif
