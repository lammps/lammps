#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c_n1 = -1;
int dorgbr_(char *vect, integer *m, integer *n, integer *k, doublereal *a, integer *lda,
            doublereal *tau, doublereal *work, integer *lwork, integer *info, ftnlen vect_len)
{
    integer a_dim1, a_offset, i__1, i__2, i__3;
    integer i__, j, mn;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer iinfo;
    logical wantq;
    extern int xerbla_(char *, integer *, ftnlen),
        dorglq_(integer *, integer *, integer *, doublereal *, integer *, doublereal *,
                doublereal *, integer *, integer *),
        dorgqr_(integer *, integer *, integer *, doublereal *, integer *, doublereal *,
                doublereal *, integer *, integer *);
    integer lwkopt;
    logical lquery;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    *info = 0;
    wantq = lsame_(vect, (char *)"Q", (ftnlen)1, (ftnlen)1);
    mn = min(*m, *n);
    lquery = *lwork == -1;
    if (!wantq && !lsame_(vect, (char *)"P", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (*m < 0) {
        *info = -2;
    } else if (*n < 0 || wantq && (*n > *m || *n < min(*m, *k)) ||
               !wantq && (*m > *n || *m < min(*n, *k))) {
        *info = -3;
    } else if (*k < 0) {
        *info = -4;
    } else if (*lda < max(1, *m)) {
        *info = -6;
    } else if (*lwork < max(1, mn) && !lquery) {
        *info = -9;
    }
    if (*info == 0) {
        work[1] = 1.;
        if (wantq) {
            if (*m >= *k) {
                dorgqr_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, &iinfo);
            } else {
                if (*m > 1) {
                    i__1 = *m - 1;
                    i__2 = *m - 1;
                    i__3 = *m - 1;
                    dorgqr_(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], &work[1], &c_n1,
                            &iinfo);
                }
            }
        } else {
            if (*k < *n) {
                dorglq_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, &iinfo);
            } else {
                if (*n > 1) {
                    i__1 = *n - 1;
                    i__2 = *n - 1;
                    i__3 = *n - 1;
                    dorglq_(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], &work[1], &c_n1,
                            &iinfo);
                }
            }
        }
        lwkopt = (integer)work[1];
        lwkopt = max(lwkopt, mn);
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DORGBR", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        work[1] = (doublereal)lwkopt;
        return 0;
    }
    if (*m == 0 || *n == 0) {
        work[1] = 1.;
        return 0;
    }
    if (wantq) {
        if (*m >= *k) {
            dorgqr_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], lwork, &iinfo);
        } else {
            for (j = *m; j >= 2; --j) {
                a[j * a_dim1 + 1] = 0.;
                i__1 = *m;
                for (i__ = j + 1; i__ <= i__1; ++i__) {
                    a[i__ + j * a_dim1] = a[i__ + (j - 1) * a_dim1];
                }
            }
            a[a_dim1 + 1] = 1.;
            i__1 = *m;
            for (i__ = 2; i__ <= i__1; ++i__) {
                a[i__ + a_dim1] = 0.;
            }
            if (*m > 1) {
                i__1 = *m - 1;
                i__2 = *m - 1;
                i__3 = *m - 1;
                dorgqr_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[1], &work[1], lwork,
                        &iinfo);
            }
        }
    } else {
        if (*k < *n) {
            dorglq_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], lwork, &iinfo);
        } else {
            a[a_dim1 + 1] = 1.;
            i__1 = *n;
            for (i__ = 2; i__ <= i__1; ++i__) {
                a[i__ + a_dim1] = 0.;
            }
            i__1 = *n;
            for (j = 2; j <= i__1; ++j) {
                for (i__ = j - 1; i__ >= 2; --i__) {
                    a[i__ + j * a_dim1] = a[i__ - 1 + j * a_dim1];
                }
                a[j * a_dim1 + 1] = 0.;
            }
            if (*n > 1) {
                i__1 = *n - 1;
                i__2 = *n - 1;
                i__3 = *n - 1;
                dorglq_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[1], &work[1], lwork,
                        &iinfo);
            }
        }
    }
    work[1] = (doublereal)lwkopt;
    return 0;
}
#ifdef __cplusplus
}
#endif
