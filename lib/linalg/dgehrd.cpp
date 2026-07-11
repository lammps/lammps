#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
static integer c__65 = 65;
static doublereal c_b25 = -1.;
static doublereal c_b26 = 1.;
int dgehrd_(integer *n, integer *ilo, integer *ihi, doublereal *a, integer *lda, doublereal *tau,
            doublereal *work, integer *lwork, integer *info)
{
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    integer i__, j, ib;
    doublereal ei;
    integer nb, nh, nx, iwt;
    extern int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
                      integer *, doublereal *, integer *, doublereal *, doublereal *, integer *,
                      ftnlen, ftnlen);
    integer nbmin, iinfo;
    extern int dtrmm_(char *, char *, char *, char *, integer *, integer *, doublereal *,
                      doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen,
                      ftnlen),
        daxpy_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *),
        dgehd2_(integer *, integer *, integer *, doublereal *, integer *, doublereal *,
                doublereal *, integer *),
        dlahr2_(integer *, integer *, integer *, doublereal *, integer *, doublereal *,
                doublereal *, integer *, doublereal *, integer *),
        dlarfb_(char *, char *, char *, char *, integer *, integer *, integer *, doublereal *,
                integer *, doublereal *, integer *, doublereal *, integer *, doublereal *,
                integer *, ftnlen, ftnlen, ftnlen, ftnlen),
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
    lquery = *lwork == -1;
    if (*n < 0) {
        *info = -1;
    } else if (*ilo < 1 || *ilo > max(1, *n)) {
        *info = -2;
    } else if (*ihi < min(*ilo, *n) || *ihi > *n) {
        *info = -3;
    } else if (*lda < max(1, *n)) {
        *info = -5;
    } else if (*lwork < max(1, *n) && !lquery) {
        *info = -8;
    }
    nh = *ihi - *ilo + 1;
    if (*info == 0) {
        if (nh <= 1) {
            lwkopt = 1;
        } else {
            i__1 = 64,
            i__2 = ilaenv_(&c__1, (char *)"DGEHRD", (char *)" ", n, ilo, ihi, &c_n1, (ftnlen)6, (ftnlen)1);
            nb = min(i__1, i__2);
            lwkopt = *n * nb + 4160;
        }
        work[1] = (doublereal)lwkopt;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DGEHRD", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    i__1 = *ilo - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
        tau[i__] = 0.;
    }
    i__1 = *n - 1;
    for (i__ = max(1, *ihi); i__ <= i__1; ++i__) {
        tau[i__] = 0.;
    }
    if (nh <= 1) {
        work[1] = 1.;
        return 0;
    }
    i__1 = 64, i__2 = ilaenv_(&c__1, (char *)"DGEHRD", (char *)" ", n, ilo, ihi, &c_n1, (ftnlen)6, (ftnlen)1);
    nb = min(i__1, i__2);
    nbmin = 2;
    if (nb > 1 && nb < nh) {
        i__1 = nb, i__2 = ilaenv_(&c__3, (char *)"DGEHRD", (char *)" ", n, ilo, ihi, &c_n1, (ftnlen)6, (ftnlen)1);
        nx = max(i__1, i__2);
        if (nx < nh) {
            if (*lwork < lwkopt) {
                i__1 = 2,
                i__2 = ilaenv_(&c__2, (char *)"DGEHRD", (char *)" ", n, ilo, ihi, &c_n1, (ftnlen)6, (ftnlen)1);
                nbmin = max(i__1, i__2);
                if (*lwork >= *n * nbmin + 4160) {
                    nb = (*lwork - 4160) / *n;
                } else {
                    nb = 1;
                }
            }
        }
    }
    ldwork = *n;
    if (nb < nbmin || nb >= nh) {
        i__ = *ilo;
    } else {
        iwt = *n * nb + 1;
        i__1 = *ihi - 1 - nx;
        i__2 = nb;
        for (i__ = *ilo; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
            i__3 = nb, i__4 = *ihi - i__;
            ib = min(i__3, i__4);
            dlahr2_(ihi, &i__, &ib, &a[i__ * a_dim1 + 1], lda, &tau[i__], &work[iwt], &c__65,
                    &work[1], &ldwork);
            ei = a[i__ + ib + (i__ + ib - 1) * a_dim1];
            a[i__ + ib + (i__ + ib - 1) * a_dim1] = 1.;
            i__3 = *ihi - i__ - ib + 1;
            dgemm_((char *)"N", (char *)"T", ihi, &i__3, &ib, &c_b25, &work[1], &ldwork,
                   &a[i__ + ib + i__ * a_dim1], lda, &c_b26, &a[(i__ + ib) * a_dim1 + 1], lda,
                   (ftnlen)1, (ftnlen)1);
            a[i__ + ib + (i__ + ib - 1) * a_dim1] = ei;
            i__3 = ib - 1;
            dtrmm_((char *)"R", (char *)"L", (char *)"T", (char *)"U", &i__, &i__3, &c_b26, &a[i__ + 1 + i__ * a_dim1], lda,
                   &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
            i__3 = ib - 2;
            for (j = 0; j <= i__3; ++j) {
                daxpy_(&i__, &c_b25, &work[ldwork * j + 1], &c__1, &a[(i__ + j + 1) * a_dim1 + 1],
                       &c__1);
            }
            i__3 = *ihi - i__;
            i__4 = *n - i__ - ib + 1;
            dlarfb_((char *)"L", (char *)"T", (char *)"F", (char *)"C", &i__3, &i__4, &ib, &a[i__ + 1 + i__ * a_dim1], lda,
                    &work[iwt], &c__65, &a[i__ + 1 + (i__ + ib) * a_dim1], lda, &work[1], &ldwork,
                    (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
        }
    }
    dgehd2_(n, &i__, ihi, &a[a_offset], lda, &tau[1], &work[1], &iinfo);
    work[1] = (doublereal)lwkopt;
    return 0;
}
#ifdef __cplusplus
}
#endif
