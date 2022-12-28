#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b16 = 1.;
int dsygv_(integer *itype, char *jobz, char *uplo, integer *n, doublereal *a, integer *lda,
           doublereal *b, integer *ldb, doublereal *w, doublereal *work, integer *lwork,
           integer *info, ftnlen jobz_len, ftnlen uplo_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    integer nb, neig;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int dtrmm_(char *, char *, char *, char *, integer *, integer *, doublereal *,
                      doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen,
                      ftnlen);
    char trans[1];
    extern int dtrsm_(char *, char *, char *, char *, integer *, integer *, doublereal *,
                      doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen,
                      ftnlen);
    logical upper;
    extern int dsyev_(char *, char *, integer *, doublereal *, integer *, doublereal *,
                      doublereal *, integer *, integer *, ftnlen, ftnlen);
    logical wantz;
    extern int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    extern int dpotrf_(char *, integer *, doublereal *, integer *, integer *, ftnlen);
    integer lwkmin;
    extern int dsygst_(integer *, char *, integer *, doublereal *, integer *, doublereal *,
                       integer *, integer *, ftnlen);
    integer lwkopt;
    logical lquery;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --w;
    --work;
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
        i__1 = 1, i__2 = *n * 3 - 1;
        lwkmin = max(i__1, i__2);
        nb = ilaenv_(&c__1, (char *)"DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
        i__1 = lwkmin, i__2 = (nb + 2) * *n;
        lwkopt = max(i__1, i__2);
        work[1] = (doublereal)lwkopt;
        if (*lwork < lwkmin && !lquery) {
            *info = -11;
        }
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DSYGV ", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    dpotrf_(uplo, n, &b[b_offset], ldb, info, (ftnlen)1);
    if (*info != 0) {
        *info = *n + *info;
        return 0;
    }
    dsygst_(itype, uplo, n, &a[a_offset], lda, &b[b_offset], ldb, info, (ftnlen)1);
    dsyev_(jobz, uplo, n, &a[a_offset], lda, &w[1], &work[1], lwork, info, (ftnlen)1, (ftnlen)1);
    if (wantz) {
        neig = *n;
        if (*info > 0) {
            neig = *info - 1;
        }
        if (*itype == 1 || *itype == 2) {
            if (upper) {
                *(unsigned char *)trans = 'N';
            } else {
                *(unsigned char *)trans = 'T';
            }
            dtrsm_((char *)"Left", uplo, trans, (char *)"Non-unit", n, &neig, &c_b16, &b[b_offset], ldb,
                   &a[a_offset], lda, (ftnlen)4, (ftnlen)1, (ftnlen)1, (ftnlen)8);
        } else if (*itype == 3) {
            if (upper) {
                *(unsigned char *)trans = 'T';
            } else {
                *(unsigned char *)trans = 'N';
            }
            dtrmm_((char *)"Left", uplo, trans, (char *)"Non-unit", n, &neig, &c_b16, &b[b_offset], ldb,
                   &a[a_offset], lda, (ftnlen)4, (ftnlen)1, (ftnlen)1, (ftnlen)8);
        }
    }
    work[1] = (doublereal)lwkopt;
    return 0;
}
#ifdef __cplusplus
}
#endif
