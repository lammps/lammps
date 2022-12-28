#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b16 = 1.;
static doublereal c_b19 = -1.;
int dgetrf_(integer *m, integer *n, doublereal *a, integer *lda, integer *ipiv, integer *info)
{
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    integer i__, j, jb, nb;
    extern int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
                      integer *, doublereal *, integer *, doublereal *, doublereal *, integer *,
                      ftnlen, ftnlen);
    integer iinfo;
    extern int dtrsm_(char *, char *, char *, char *, integer *, integer *, doublereal *,
                      doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen,
                      ftnlen),
        xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    extern int dlaswp_(integer *, doublereal *, integer *, integer *, integer *, integer *,
                       integer *),
        dgetrf2_(integer *, integer *, doublereal *, integer *, integer *, integer *);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    *info = 0;
    if (*m < 0) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*lda < max(1, *m)) {
        *info = -4;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DGETRF", &i__1, (ftnlen)6);
        return 0;
    }
    if (*m == 0 || *n == 0) {
        return 0;
    }
    nb = ilaenv_(&c__1, (char *)"DGETRF", (char *)" ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
    if (nb <= 1 || nb >= min(*m, *n)) {
        dgetrf2_(m, n, &a[a_offset], lda, &ipiv[1], info);
    } else {
        i__1 = min(*m, *n);
        i__2 = nb;
        for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
            i__3 = min(*m, *n) - j + 1;
            jb = min(i__3, nb);
            i__3 = *m - j + 1;
            dgetrf2_(&i__3, &jb, &a[j + j * a_dim1], lda, &ipiv[j], &iinfo);
            if (*info == 0 && iinfo > 0) {
                *info = iinfo + j - 1;
            }
            i__4 = *m, i__5 = j + jb - 1;
            i__3 = min(i__4, i__5);
            for (i__ = j; i__ <= i__3; ++i__) {
                ipiv[i__] = j - 1 + ipiv[i__];
            }
            i__3 = j - 1;
            i__4 = j + jb - 1;
            dlaswp_(&i__3, &a[a_offset], lda, &j, &i__4, &ipiv[1], &c__1);
            if (j + jb <= *n) {
                i__3 = *n - j - jb + 1;
                i__4 = j + jb - 1;
                dlaswp_(&i__3, &a[(j + jb) * a_dim1 + 1], lda, &j, &i__4, &ipiv[1], &c__1);
                i__3 = *n - j - jb + 1;
                dtrsm_((char *)"Left", (char *)"Lower", (char *)"No transpose", (char *)"Unit", &jb, &i__3, &c_b16,
                       &a[j + j * a_dim1], lda, &a[j + (j + jb) * a_dim1], lda, (ftnlen)4,
                       (ftnlen)5, (ftnlen)12, (ftnlen)4);
                if (j + jb <= *m) {
                    i__3 = *m - j - jb + 1;
                    i__4 = *n - j - jb + 1;
                    dgemm_((char *)"No transpose", (char *)"No transpose", &i__3, &i__4, &jb, &c_b19,
                           &a[j + jb + j * a_dim1], lda, &a[j + (j + jb) * a_dim1], lda, &c_b16,
                           &a[j + jb + (j + jb) * a_dim1], lda, (ftnlen)12, (ftnlen)12);
                }
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
