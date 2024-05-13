#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static doublereal c_b18 = 1.;
static doublereal c_b22 = -1.;
int dtrtri_(char *uplo, char *diag, integer *n, doublereal *a, integer *lda, integer *info,
            ftnlen uplo_len, ftnlen diag_len)
{
    address a__1[2];
    integer a_dim1, a_offset, i__1, i__2[2], i__3, i__4, i__5;
    char ch__1[2];
    int s_lmp_cat(char *, char **, integer *, integer *, ftnlen);
    integer j, jb, nb, nn;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int dtrmm_(char *, char *, char *, char *, integer *, integer *, doublereal *,
                      doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen,
                      ftnlen),
        dtrsm_(char *, char *, char *, char *, integer *, integer *, doublereal *, doublereal *,
               integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    logical upper;
    extern int dtrti2_(char *, char *, integer *, doublereal *, integer *, integer *, ftnlen,
                       ftnlen),
        xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    logical nounit;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    *info = 0;
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    nounit = lsame_(diag, (char *)"N", (ftnlen)1, (ftnlen)1);
    if (!upper && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (!nounit && !lsame_(diag, (char *)"U", (ftnlen)1, (ftnlen)1)) {
        *info = -2;
    } else if (*n < 0) {
        *info = -3;
    } else if (*lda < max(1, *n)) {
        *info = -5;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DTRTRI", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    if (nounit) {
        i__1 = *n;
        for (*info = 1; *info <= i__1; ++(*info)) {
            if (a[*info + *info * a_dim1] == 0.) {
                return 0;
            }
        }
        *info = 0;
    }
    i__2[0] = 1, a__1[0] = uplo;
    i__2[1] = 1, a__1[1] = diag;
    s_lmp_cat(ch__1, a__1, i__2, &c__2, (ftnlen)2);
    nb = ilaenv_(&c__1, (char *)"DTRTRI", ch__1, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)2);
    if (nb <= 1 || nb >= *n) {
        dtrti2_(uplo, diag, n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)1);
    } else {
        if (upper) {
            i__1 = *n;
            i__3 = nb;
            for (j = 1; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {
                i__4 = nb, i__5 = *n - j + 1;
                jb = min(i__4, i__5);
                i__4 = j - 1;
                dtrmm_((char *)"Left", (char *)"Upper", (char *)"No transpose", diag, &i__4, &jb, &c_b18, &a[a_offset], lda,
                       &a[j * a_dim1 + 1], lda, (ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)1);
                i__4 = j - 1;
                dtrsm_((char *)"Right", (char *)"Upper", (char *)"No transpose", diag, &i__4, &jb, &c_b22,
                       &a[j + j * a_dim1], lda, &a[j * a_dim1 + 1], lda, (ftnlen)5, (ftnlen)5,
                       (ftnlen)12, (ftnlen)1);
                dtrti2_((char *)"Upper", diag, &jb, &a[j + j * a_dim1], lda, info, (ftnlen)5, (ftnlen)1);
            }
        } else {
            nn = (*n - 1) / nb * nb + 1;
            i__3 = -nb;
            for (j = nn; i__3 < 0 ? j >= 1 : j <= 1; j += i__3) {
                i__1 = nb, i__4 = *n - j + 1;
                jb = min(i__1, i__4);
                if (j + jb <= *n) {
                    i__1 = *n - j - jb + 1;
                    dtrmm_((char *)"Left", (char *)"Lower", (char *)"No transpose", diag, &i__1, &jb, &c_b18,
                           &a[j + jb + (j + jb) * a_dim1], lda, &a[j + jb + j * a_dim1], lda,
                           (ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)1);
                    i__1 = *n - j - jb + 1;
                    dtrsm_((char *)"Right", (char *)"Lower", (char *)"No transpose", diag, &i__1, &jb, &c_b22,
                           &a[j + j * a_dim1], lda, &a[j + jb + j * a_dim1], lda, (ftnlen)5,
                           (ftnlen)5, (ftnlen)12, (ftnlen)1);
                }
                dtrti2_((char *)"Lower", diag, &jb, &a[j + j * a_dim1], lda, info, (ftnlen)5, (ftnlen)1);
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
