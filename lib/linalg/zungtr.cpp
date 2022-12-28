#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
int zungtr_(char *uplo, integer *n, doublecomplex *a, integer *lda, doublecomplex *tau,
            doublecomplex *work, integer *lwork, integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    integer i__, j, nb;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer iinfo;
    logical upper;
    extern int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    integer lwkopt;
    logical lquery;
    extern int zungql_(integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *,
                       doublecomplex *, integer *, integer *),
        zungqr_(integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *,
                doublecomplex *, integer *, integer *);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    *info = 0;
    lquery = *lwork == -1;
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    if (!upper && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*lda < max(1, *n)) {
        *info = -4;
    } else {
        i__1 = 1, i__2 = *n - 1;
        if (*lwork < max(i__1, i__2) && !lquery) {
            *info = -7;
        }
    }
    if (*info == 0) {
        if (upper) {
            i__1 = *n - 1;
            i__2 = *n - 1;
            i__3 = *n - 1;
            nb = ilaenv_(&c__1, (char *)"ZUNGQL", (char *)" ", &i__1, &i__2, &i__3, &c_n1, (ftnlen)6, (ftnlen)1);
        } else {
            i__1 = *n - 1;
            i__2 = *n - 1;
            i__3 = *n - 1;
            nb = ilaenv_(&c__1, (char *)"ZUNGQR", (char *)" ", &i__1, &i__2, &i__3, &c_n1, (ftnlen)6, (ftnlen)1);
        }
        i__1 = 1, i__2 = *n - 1;
        lwkopt = max(i__1, i__2) * nb;
        work[1].r = (doublereal)lwkopt, work[1].i = 0.;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"ZUNGTR", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    if (*n == 0) {
        work[1].r = 1., work[1].i = 0.;
        return 0;
    }
    if (upper) {
        i__1 = *n - 1;
        for (j = 1; j <= i__1; ++j) {
            i__2 = j - 1;
            for (i__ = 1; i__ <= i__2; ++i__) {
                i__3 = i__ + j * a_dim1;
                i__4 = i__ + (j + 1) * a_dim1;
                a[i__3].r = a[i__4].r, a[i__3].i = a[i__4].i;
            }
            i__2 = *n + j * a_dim1;
            a[i__2].r = 0., a[i__2].i = 0.;
        }
        i__1 = *n - 1;
        for (i__ = 1; i__ <= i__1; ++i__) {
            i__2 = i__ + *n * a_dim1;
            a[i__2].r = 0., a[i__2].i = 0.;
        }
        i__1 = *n + *n * a_dim1;
        a[i__1].r = 1., a[i__1].i = 0.;
        i__1 = *n - 1;
        i__2 = *n - 1;
        i__3 = *n - 1;
        zungql_(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], &work[1], lwork, &iinfo);
    } else {
        for (j = *n; j >= 2; --j) {
            i__1 = j * a_dim1 + 1;
            a[i__1].r = 0., a[i__1].i = 0.;
            i__1 = *n;
            for (i__ = j + 1; i__ <= i__1; ++i__) {
                i__2 = i__ + j * a_dim1;
                i__3 = i__ + (j - 1) * a_dim1;
                a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
            }
        }
        i__1 = a_dim1 + 1;
        a[i__1].r = 1., a[i__1].i = 0.;
        i__1 = *n;
        for (i__ = 2; i__ <= i__1; ++i__) {
            i__2 = i__ + a_dim1;
            a[i__2].r = 0., a[i__2].i = 0.;
        }
        if (*n > 1) {
            i__1 = *n - 1;
            i__2 = *n - 1;
            i__3 = *n - 1;
            zungqr_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[1], &work[1], lwork,
                    &iinfo);
        }
    }
    work[1].r = (doublereal)lwkopt, work[1].i = 0.;
    return 0;
}
#ifdef __cplusplus
}
#endif
