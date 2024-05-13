#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static doublereal c_b20 = -1.;
static doublereal c_b22 = 1.;
int dgetri_(integer *n, doublereal *a, integer *lda, integer *ipiv, doublereal *work,
            integer *lwork, integer *info)
{
    integer a_dim1, a_offset, i__1, i__2, i__3;
    integer i__, j, jb, nb, jj, jp, nn, iws;
    extern int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
                      integer *, doublereal *, integer *, doublereal *, doublereal *, integer *,
                      ftnlen, ftnlen),
        dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *, doublereal *,
               integer *, doublereal *, doublereal *, integer *, ftnlen);
    integer nbmin;
    extern int dswap_(integer *, doublereal *, integer *, doublereal *, integer *),
        dtrsm_(char *, char *, char *, char *, integer *, integer *, doublereal *, doublereal *,
               integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen),
        xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    integer ldwork;
    extern int dtrtri_(char *, char *, integer *, doublereal *, integer *, integer *, ftnlen,
                       ftnlen);
    integer lwkopt;
    logical lquery;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    --work;
    *info = 0;
    nb = ilaenv_(&c__1, (char *)"DGETRI", (char *)" ", n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
    lwkopt = *n * nb;
    work[1] = (doublereal)lwkopt;
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
        xerbla_((char *)"DGETRI", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    dtrtri_((char *)"Upper", (char *)"Non-unit", n, &a[a_offset], lda, info, (ftnlen)5, (ftnlen)8);
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
            i__2 = ilaenv_(&c__2, (char *)"DGETRI", (char *)" ", n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
            nbmin = max(i__1, i__2);
        }
    } else {
        iws = *n;
    }
    if (nb < nbmin || nb >= *n) {
        for (j = *n; j >= 1; --j) {
            i__1 = *n;
            for (i__ = j + 1; i__ <= i__1; ++i__) {
                work[i__] = a[i__ + j * a_dim1];
                a[i__ + j * a_dim1] = 0.;
            }
            if (j < *n) {
                i__1 = *n - j;
                dgemv_((char *)"No transpose", n, &i__1, &c_b20, &a[(j + 1) * a_dim1 + 1], lda,
                       &work[j + 1], &c__1, &c_b22, &a[j * a_dim1 + 1], &c__1, (ftnlen)12);
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
                    work[i__ + (jj - j) * ldwork] = a[i__ + jj * a_dim1];
                    a[i__ + jj * a_dim1] = 0.;
                }
            }
            if (j + jb <= *n) {
                i__2 = *n - j - jb + 1;
                dgemm_((char *)"No transpose", (char *)"No transpose", n, &jb, &i__2, &c_b20,
                       &a[(j + jb) * a_dim1 + 1], lda, &work[j + jb], &ldwork, &c_b22,
                       &a[j * a_dim1 + 1], lda, (ftnlen)12, (ftnlen)12);
            }
            dtrsm_((char *)"Right", (char *)"Lower", (char *)"No transpose", (char *)"Unit", n, &jb, &c_b22, &work[j], &ldwork,
                   &a[j * a_dim1 + 1], lda, (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
        }
    }
    for (j = *n - 1; j >= 1; --j) {
        jp = ipiv[j];
        if (jp != j) {
            dswap_(n, &a[j * a_dim1 + 1], &c__1, &a[jp * a_dim1 + 1], &c__1);
        }
    }
    work[1] = (doublereal)iws;
    return 0;
}
#ifdef __cplusplus
}
#endif
