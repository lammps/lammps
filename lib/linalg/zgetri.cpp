#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublecomplex c_b2 = {1., 0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
int zgetri_(integer *n, doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *work,
            integer *lwork, integer *info)
{
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1;
    integer i__, j, jb, nb, jj, jp, nn, iws, nbmin;
    extern int zgemm_(char *, char *, integer *, integer *, integer *, doublecomplex *,
                      doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *,
                      doublecomplex *, integer *, ftnlen, ftnlen),
        zgemv_(char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *,
               doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *, ftnlen),
        zswap_(integer *, doublecomplex *, integer *, doublecomplex *, integer *),
        ztrsm_(char *, char *, char *, char *, integer *, integer *, doublecomplex *,
               doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, ftnlen, ftnlen,
               ftnlen),
        xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    integer ldwork, lwkopt;
    logical lquery;
    extern int ztrtri_(char *, char *, integer *, doublecomplex *, integer *, integer *, ftnlen,
                       ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    --work;
    *info = 0;
    nb = ilaenv_(&c__1, (char *)"ZGETRI", (char *)" ", n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
    i__1 = 1, i__2 = *n * nb;
    lwkopt = max(i__1, i__2);
    work[1].r = (doublereal)lwkopt, work[1].i = 0.;
    lquery = *lwork == -1;
    if (*n < 0) {
        *info = -1;
    } else if (*lda < max(1, *n)) {
        *info = -3;
    } else if (*lwork < max(1, *n) && !lquery) {
        *info = -6;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"ZGETRI", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    ztrtri_((char *)"U", (char *)"N", n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)1);
    if (*info > 0) {
        return 0;
    }
    nbmin = 2;
    ldwork = *n;
    if (nb > 1 && nb < *n) {
        i__1 = ldwork * nb;
        iws = max(i__1, 1);
        if (*lwork < iws) {
            nb = *lwork / ldwork;
            i__1 = 2,
            i__2 = ilaenv_(&c__2, (char *)"ZGETRI", (char *)" ", n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
            nbmin = max(i__1, i__2);
        }
    } else {
        iws = *n;
    }
    if (nb < nbmin || nb >= *n) {
        for (j = *n; j >= 1; --j) {
            i__1 = *n;
            for (i__ = j + 1; i__ <= i__1; ++i__) {
                i__2 = i__;
                i__3 = i__ + j * a_dim1;
                work[i__2].r = a[i__3].r, work[i__2].i = a[i__3].i;
                i__2 = i__ + j * a_dim1;
                a[i__2].r = 0., a[i__2].i = 0.;
            }
            if (j < *n) {
                i__1 = *n - j;
                z__1.r = -1., z__1.i = -0.;
                zgemv_((char *)"N", n, &i__1, &z__1, &a[(j + 1) * a_dim1 + 1], lda, &work[j + 1], &c__1,
                       &c_b2, &a[j * a_dim1 + 1], &c__1, (ftnlen)1);
            }
        }
    } else {
        nn = (*n - 1) / nb * nb + 1;
        i__1 = -nb;
        for (j = nn; i__1 < 0 ? j >= 1 : j <= 1; j += i__1) {
            i__2 = nb, i__3 = *n - j + 1;
            jb = min(i__2, i__3);
            i__2 = j + jb - 1;
            for (jj = j; jj <= i__2; ++jj) {
                i__3 = *n;
                for (i__ = jj + 1; i__ <= i__3; ++i__) {
                    i__4 = i__ + (jj - j) * ldwork;
                    i__5 = i__ + jj * a_dim1;
                    work[i__4].r = a[i__5].r, work[i__4].i = a[i__5].i;
                    i__4 = i__ + jj * a_dim1;
                    a[i__4].r = 0., a[i__4].i = 0.;
                }
            }
            if (j + jb <= *n) {
                i__2 = *n - j - jb + 1;
                z__1.r = -1., z__1.i = -0.;
                zgemm_((char *)"N", (char *)"N", n, &jb, &i__2, &z__1, &a[(j + jb) * a_dim1 + 1], lda,
                       &work[j + jb], &ldwork, &c_b2, &a[j * a_dim1 + 1], lda, (ftnlen)1,
                       (ftnlen)1);
            }
            ztrsm_((char *)"R", (char *)"L", (char *)"N", (char *)"U", n, &jb, &c_b2, &work[j], &ldwork, &a[j * a_dim1 + 1], lda,
                   (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
        }
    }
    for (j = *n - 1; j >= 1; --j) {
        jp = ipiv[j];
        if (jp != j) {
            zswap_(n, &a[j * a_dim1 + 1], &c__1, &a[jp * a_dim1 + 1], &c__1);
        }
    }
    work[1].r = (doublereal)iws, work[1].i = 0.;
    return 0;
}
#ifdef __cplusplus
}
#endif
