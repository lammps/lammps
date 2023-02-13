#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__3 = 3;
static integer c__2 = 2;
static doublereal c_b22 = -1.;
static doublereal c_b23 = 1.;
int dsytrd_(char *uplo, integer *n, doublereal *a, integer *lda, doublereal *d__, doublereal *e,
            doublereal *tau, doublereal *work, integer *lwork, integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, i__1, i__2, i__3;
    integer i__, j, nb, kk, nx, iws;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer nbmin, iinfo;
    logical upper;
    extern int dsytd2_(char *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                       doublereal *, integer *, ftnlen),
        dsyr2k_(char *, char *, integer *, integer *, doublereal *, doublereal *, integer *,
                doublereal *, integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
        dlatrd_(char *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                doublereal *, integer *, ftnlen),
        xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    integer ldwork, lwkopt;
    logical lquery;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tau;
    --work;
    *info = 0;
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    lquery = *lwork == -1;
    if (!upper && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*lda < max(1, *n)) {
        *info = -4;
    } else if (*lwork < 1 && !lquery) {
        *info = -9;
    }
    if (*info == 0) {
        nb = ilaenv_(&c__1, (char *)"DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
        lwkopt = *n * nb;
        work[1] = (doublereal)lwkopt;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DSYTRD", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    if (*n == 0) {
        work[1] = 1.;
        return 0;
    }
    nx = *n;
    iws = 1;
    if (nb > 1 && nb < *n) {
        i__1 = nb,
        i__2 = ilaenv_(&c__3, (char *)"DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
        nx = max(i__1, i__2);
        if (nx < *n) {
            ldwork = *n;
            iws = ldwork * nb;
            if (*lwork < iws) {
                i__1 = *lwork / ldwork;
                nb = max(i__1, 1);
                nbmin =
                    ilaenv_(&c__2, (char *)"DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
                if (nb < nbmin) {
                    nx = *n;
                }
            }
        } else {
            nx = *n;
        }
    } else {
        nb = 1;
    }
    if (upper) {
        kk = *n - (*n - nx + nb - 1) / nb * nb;
        i__1 = kk + 1;
        i__2 = -nb;
        for (i__ = *n - nb + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
            i__3 = i__ + nb - 1;
            dlatrd_(uplo, &i__3, &nb, &a[a_offset], lda, &e[1], &tau[1], &work[1], &ldwork,
                    (ftnlen)1);
            i__3 = i__ - 1;
            dsyr2k_(uplo, (char *)"No transpose", &i__3, &nb, &c_b22, &a[i__ * a_dim1 + 1], lda, &work[1],
                    &ldwork, &c_b23, &a[a_offset], lda, (ftnlen)1, (ftnlen)12);
            i__3 = i__ + nb - 1;
            for (j = i__; j <= i__3; ++j) {
                a[j - 1 + j * a_dim1] = e[j - 1];
                d__[j] = a[j + j * a_dim1];
            }
        }
        dsytd2_(uplo, &kk, &a[a_offset], lda, &d__[1], &e[1], &tau[1], &iinfo, (ftnlen)1);
    } else {
        i__2 = *n - nx;
        i__1 = nb;
        for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
            i__3 = *n - i__ + 1;
            dlatrd_(uplo, &i__3, &nb, &a[i__ + i__ * a_dim1], lda, &e[i__], &tau[i__], &work[1],
                    &ldwork, (ftnlen)1);
            i__3 = *n - i__ - nb + 1;
            dsyr2k_(uplo, (char *)"No transpose", &i__3, &nb, &c_b22, &a[i__ + nb + i__ * a_dim1], lda,
                    &work[nb + 1], &ldwork, &c_b23, &a[i__ + nb + (i__ + nb) * a_dim1], lda,
                    (ftnlen)1, (ftnlen)12);
            i__3 = i__ + nb - 1;
            for (j = i__; j <= i__3; ++j) {
                a[j + 1 + j * a_dim1] = e[j];
                d__[j] = a[j + j * a_dim1];
            }
        }
        i__1 = *n - i__ + 1;
        dsytd2_(uplo, &i__1, &a[i__ + i__ * a_dim1], lda, &d__[i__], &e[i__], &tau[i__], &iinfo,
                (ftnlen)1);
    }
    work[1] = (doublereal)lwkopt;
    return 0;
}
#ifdef __cplusplus
}
#endif
