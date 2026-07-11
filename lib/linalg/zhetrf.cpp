#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
int zhetrf_(char *uplo, integer *n, doublecomplex *a, integer *lda, integer *ipiv,
            doublecomplex *work, integer *lwork, integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, i__1, i__2;
    integer j, k, kb, nb, iws;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer nbmin, iinfo;
    logical upper;
    extern int zhetf2_(char *, integer *, doublecomplex *, integer *, integer *, integer *, ftnlen),
        zlahef_(char *, integer *, integer *, integer *, doublecomplex *, integer *, integer *,
                doublecomplex *, integer *, integer *, ftnlen),
        xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    integer ldwork, lwkopt;
    logical lquery;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
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
        *info = -7;
    }
    if (*info == 0) {
        nb = ilaenv_(&c__1, (char *)"ZHETRF", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
        i__1 = 1, i__2 = *n * nb;
        lwkopt = max(i__1, i__2);
        work[1].r = (doublereal)lwkopt, work[1].i = 0.;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"ZHETRF", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    nbmin = 2;
    ldwork = *n;
    if (nb > 1 && nb < *n) {
        iws = ldwork * nb;
        if (*lwork < iws) {
            i__1 = *lwork / ldwork;
            nb = max(i__1, 1);
            i__1 = 2,
            i__2 = ilaenv_(&c__2, (char *)"ZHETRF", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
            nbmin = max(i__1, i__2);
        }
    } else {
        iws = 1;
    }
    if (nb < nbmin) {
        nb = *n;
    }
    if (upper) {
        k = *n;
    L10:
        if (k < 1) {
            goto L40;
        }
        if (k > nb) {
            zlahef_(uplo, &k, &nb, &kb, &a[a_offset], lda, &ipiv[1], &work[1], n, &iinfo,
                    (ftnlen)1);
        } else {
            zhetf2_(uplo, &k, &a[a_offset], lda, &ipiv[1], &iinfo, (ftnlen)1);
            kb = k;
        }
        if (*info == 0 && iinfo > 0) {
            *info = iinfo;
        }
        k -= kb;
        goto L10;
    } else {
        k = 1;
    L20:
        if (k > *n) {
            goto L40;
        }
        if (k <= *n - nb) {
            i__1 = *n - k + 1;
            zlahef_(uplo, &i__1, &nb, &kb, &a[k + k * a_dim1], lda, &ipiv[k], &work[1], n, &iinfo,
                    (ftnlen)1);
        } else {
            i__1 = *n - k + 1;
            zhetf2_(uplo, &i__1, &a[k + k * a_dim1], lda, &ipiv[k], &iinfo, (ftnlen)1);
            kb = *n - k + 1;
        }
        if (*info == 0 && iinfo > 0) {
            *info = iinfo + k - 1;
        }
        i__1 = k + kb - 1;
        for (j = k; j <= i__1; ++j) {
            if (ipiv[j] > 0) {
                ipiv[j] = ipiv[j] + k - 1;
            } else {
                ipiv[j] = ipiv[j] - k + 1;
            }
        }
        k += kb;
        goto L20;
    }
L40:
    work[1].r = (doublereal)lwkopt, work[1].i = 0.;
    return 0;
}
#ifdef __cplusplus
}
#endif
