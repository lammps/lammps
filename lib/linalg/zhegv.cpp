#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublecomplex c_b1 = {1., 0.};
static integer c__1 = 1;
static integer c_n1 = -1;
int zhegv_(integer *itype, char *jobz, char *uplo, integer *n, doublecomplex *a, integer *lda,
           doublecomplex *b, integer *ldb, doublereal *w, doublecomplex *work, integer *lwork,
           doublereal *rwork, integer *info, ftnlen jobz_len, ftnlen uplo_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    integer nb, neig;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int zheev_(char *, char *, integer *, doublecomplex *, integer *, doublereal *,
                      doublecomplex *, integer *, doublereal *, integer *, ftnlen, ftnlen);
    char trans[1];
    logical upper, wantz;
    extern int ztrmm_(char *, char *, char *, char *, integer *, integer *, doublecomplex *,
                      doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, ftnlen,
                      ftnlen, ftnlen),
        ztrsm_(char *, char *, char *, char *, integer *, integer *, doublecomplex *,
               doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, ftnlen, ftnlen,
               ftnlen),
        xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    extern int zhegst_(integer *, char *, integer *, doublecomplex *, integer *, doublecomplex *,
                       integer *, integer *, ftnlen);
    integer lwkopt;
    logical lquery;
    extern int zpotrf_(char *, integer *, doublecomplex *, integer *, integer *, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --w;
    --work;
    --rwork;
    wantz = lsame_(jobz, (char *)"V", (ftnlen)1, (ftnlen)1);
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    lquery = *lwork == -1;
    *info = 0;
    if (*itype < 1 || *itype > 3) {
        *info = -1;
    } else if (!(wantz || lsame_(jobz, (char *)"N", (ftnlen)1, (ftnlen)1))) {
        *info = -2;
    } else if (!(upper || lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1))) {
        *info = -3;
    } else if (*n < 0) {
        *info = -4;
    } else if (*lda < max(1, *n)) {
        *info = -6;
    } else if (*ldb < max(1, *n)) {
        *info = -8;
    }
    if (*info == 0) {
        nb = ilaenv_(&c__1, (char *)"ZHETRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
        i__1 = 1, i__2 = (nb + 1) * *n;
        lwkopt = max(i__1, i__2);
        work[1].r = (doublereal)lwkopt, work[1].i = 0.;
        i__1 = 1, i__2 = (*n << 1) - 1;
        if (*lwork < max(i__1, i__2) && !lquery) {
            *info = -11;
        }
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"ZHEGV ", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    zpotrf_(uplo, n, &b[b_offset], ldb, info, (ftnlen)1);
    if (*info != 0) {
        *info = *n + *info;
        return 0;
    }
    zhegst_(itype, uplo, n, &a[a_offset], lda, &b[b_offset], ldb, info, (ftnlen)1);
    zheev_(jobz, uplo, n, &a[a_offset], lda, &w[1], &work[1], lwork, &rwork[1], info, (ftnlen)1,
           (ftnlen)1);
    if (wantz) {
        neig = *n;
        if (*info > 0) {
            neig = *info - 1;
        }
        if (*itype == 1 || *itype == 2) {
            if (upper) {
                *(unsigned char *)trans = 'N';
            } else {
                *(unsigned char *)trans = 'C';
            }
            ztrsm_((char *)"L", uplo, trans, (char *)"N", n, &neig, &c_b1, &b[b_offset], ldb, &a[a_offset], lda,
                   (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
        } else if (*itype == 3) {
            if (upper) {
                *(unsigned char *)trans = 'C';
            } else {
                *(unsigned char *)trans = 'N';
            }
            ztrmm_((char *)"L", uplo, trans, (char *)"N", n, &neig, &c_b1, &b[b_offset], ldb, &a[a_offset], lda,
                   (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
        }
    }
    work[1].r = (doublereal)lwkopt, work[1].i = 0.;
    return 0;
}
#ifdef __cplusplus
}
#endif
