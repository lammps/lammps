#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b11 = 1.;
int dsygvd_(integer *itype, char *jobz, char *uplo, integer *n, doublereal *a, integer *lda,
            doublereal *b, integer *ldb, doublereal *w, doublereal *work, integer *lwork,
            integer *iwork, integer *liwork, integer *info, ftnlen jobz_len, ftnlen uplo_len)
{
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    doublereal d__1, d__2;
    integer lopt;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int dtrmm_(char *, char *, char *, char *, integer *, integer *, doublereal *,
                      doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen,
                      ftnlen);
    integer lwmin;
    char trans[1];
    integer liopt;
    extern int dtrsm_(char *, char *, char *, char *, integer *, integer *, doublereal *,
                      doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen,
                      ftnlen);
    logical upper, wantz;
    extern int xerbla_(char *, integer *, ftnlen),
        dpotrf_(char *, integer *, doublereal *, integer *, integer *, ftnlen);
    integer liwmin;
    extern int dsyevd_(char *, char *, integer *, doublereal *, integer *, doublereal *,
                       doublereal *, integer *, integer *, integer *, integer *, ftnlen, ftnlen),
        dsygst_(integer *, char *, integer *, doublereal *, integer *, doublereal *, integer *,
                integer *, ftnlen);
    logical lquery;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --w;
    --work;
    --iwork;
    wantz = lsame_(jobz, (char *)"V", (ftnlen)1, (ftnlen)1);
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    lquery = *lwork == -1 || *liwork == -1;
    *info = 0;
    if (*n <= 1) {
        liwmin = 1;
        lwmin = 1;
    } else if (wantz) {
        liwmin = *n * 5 + 3;
        i__1 = *n;
        lwmin = *n * 6 + 1 + (i__1 * i__1 << 1);
    } else {
        liwmin = 1;
        lwmin = (*n << 1) + 1;
    }
    lopt = lwmin;
    liopt = liwmin;
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
        work[1] = (doublereal)lopt;
        iwork[1] = liopt;
        if (*lwork < lwmin && !lquery) {
            *info = -11;
        } else if (*liwork < liwmin && !lquery) {
            *info = -13;
        }
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DSYGVD", &i__1, (ftnlen)6);
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
    dsyevd_(jobz, uplo, n, &a[a_offset], lda, &w[1], &work[1], lwork, &iwork[1], liwork, info,
            (ftnlen)1, (ftnlen)1);
    d__1 = (doublereal)lopt;
    lopt = (integer)max(d__1, work[1]);
    d__1 = (doublereal)liopt, d__2 = (doublereal)iwork[1];
    liopt = (integer)max(d__1, d__2);
    if (wantz && *info == 0) {
        if (*itype == 1 || *itype == 2) {
            if (upper) {
                *(unsigned char *)trans = 'N';
            } else {
                *(unsigned char *)trans = 'T';
            }
            dtrsm_((char *)"Left", uplo, trans, (char *)"Non-unit", n, n, &c_b11, &b[b_offset], ldb, &a[a_offset],
                   lda, (ftnlen)4, (ftnlen)1, (ftnlen)1, (ftnlen)8);
        } else if (*itype == 3) {
            if (upper) {
                *(unsigned char *)trans = 'T';
            } else {
                *(unsigned char *)trans = 'N';
            }
            dtrmm_((char *)"Left", uplo, trans, (char *)"Non-unit", n, n, &c_b11, &b[b_offset], ldb, &a[a_offset],
                   lda, (ftnlen)4, (ftnlen)1, (ftnlen)1, (ftnlen)8);
        }
    }
    work[1] = (doublereal)lopt;
    iwork[1] = liopt;
    return 0;
}
#ifdef __cplusplus
}
#endif
