#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dsyconv_(char *uplo, char *way, integer *n, doublereal *a, integer *lda, integer *ipiv,
             doublereal *e, integer *info, ftnlen uplo_len, ftnlen way_len)
{
    integer a_dim1, a_offset, i__1;
    integer i__, j, ip;
    doublereal temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    logical upper;
    extern int xerbla_(char *, integer *, ftnlen);
    logical convert;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    --e;
    *info = 0;
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    convert = lsame_(way, (char *)"C", (ftnlen)1, (ftnlen)1);
    if (!upper && !lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (!convert && !lsame_(way, (char *)"R", (ftnlen)1, (ftnlen)1)) {
        *info = -2;
    } else if (*n < 0) {
        *info = -3;
    } else if (*lda < max(1, *n)) {
        *info = -5;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DSYCONV", &i__1, (ftnlen)7);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    if (upper) {
        if (convert) {
            i__ = *n;
            e[1] = 0.;
            while (i__ > 1) {
                if (ipiv[i__] < 0) {
                    e[i__] = a[i__ - 1 + i__ * a_dim1];
                    e[i__ - 1] = 0.;
                    a[i__ - 1 + i__ * a_dim1] = 0.;
                    --i__;
                } else {
                    e[i__] = 0.;
                }
                --i__;
            }
            i__ = *n;
            while (i__ >= 1) {
                if (ipiv[i__] > 0) {
                    ip = ipiv[i__];
                    if (i__ < *n) {
                        i__1 = *n;
                        for (j = i__ + 1; j <= i__1; ++j) {
                            temp = a[ip + j * a_dim1];
                            a[ip + j * a_dim1] = a[i__ + j * a_dim1];
                            a[i__ + j * a_dim1] = temp;
                        }
                    }
                } else {
                    ip = -ipiv[i__];
                    if (i__ < *n) {
                        i__1 = *n;
                        for (j = i__ + 1; j <= i__1; ++j) {
                            temp = a[ip + j * a_dim1];
                            a[ip + j * a_dim1] = a[i__ - 1 + j * a_dim1];
                            a[i__ - 1 + j * a_dim1] = temp;
                        }
                    }
                    --i__;
                }
                --i__;
            }
        } else {
            i__ = 1;
            while (i__ <= *n) {
                if (ipiv[i__] > 0) {
                    ip = ipiv[i__];
                    if (i__ < *n) {
                        i__1 = *n;
                        for (j = i__ + 1; j <= i__1; ++j) {
                            temp = a[ip + j * a_dim1];
                            a[ip + j * a_dim1] = a[i__ + j * a_dim1];
                            a[i__ + j * a_dim1] = temp;
                        }
                    }
                } else {
                    ip = -ipiv[i__];
                    ++i__;
                    if (i__ < *n) {
                        i__1 = *n;
                        for (j = i__ + 1; j <= i__1; ++j) {
                            temp = a[ip + j * a_dim1];
                            a[ip + j * a_dim1] = a[i__ - 1 + j * a_dim1];
                            a[i__ - 1 + j * a_dim1] = temp;
                        }
                    }
                }
                ++i__;
            }
            i__ = *n;
            while (i__ > 1) {
                if (ipiv[i__] < 0) {
                    a[i__ - 1 + i__ * a_dim1] = e[i__];
                    --i__;
                }
                --i__;
            }
        }
    } else {
        if (convert) {
            i__ = 1;
            e[*n] = 0.;
            while (i__ <= *n) {
                if (i__ < *n && ipiv[i__] < 0) {
                    e[i__] = a[i__ + 1 + i__ * a_dim1];
                    e[i__ + 1] = 0.;
                    a[i__ + 1 + i__ * a_dim1] = 0.;
                    ++i__;
                } else {
                    e[i__] = 0.;
                }
                ++i__;
            }
            i__ = 1;
            while (i__ <= *n) {
                if (ipiv[i__] > 0) {
                    ip = ipiv[i__];
                    if (i__ > 1) {
                        i__1 = i__ - 1;
                        for (j = 1; j <= i__1; ++j) {
                            temp = a[ip + j * a_dim1];
                            a[ip + j * a_dim1] = a[i__ + j * a_dim1];
                            a[i__ + j * a_dim1] = temp;
                        }
                    }
                } else {
                    ip = -ipiv[i__];
                    if (i__ > 1) {
                        i__1 = i__ - 1;
                        for (j = 1; j <= i__1; ++j) {
                            temp = a[ip + j * a_dim1];
                            a[ip + j * a_dim1] = a[i__ + 1 + j * a_dim1];
                            a[i__ + 1 + j * a_dim1] = temp;
                        }
                    }
                    ++i__;
                }
                ++i__;
            }
        } else {
            i__ = *n;
            while (i__ >= 1) {
                if (ipiv[i__] > 0) {
                    ip = ipiv[i__];
                    if (i__ > 1) {
                        i__1 = i__ - 1;
                        for (j = 1; j <= i__1; ++j) {
                            temp = a[i__ + j * a_dim1];
                            a[i__ + j * a_dim1] = a[ip + j * a_dim1];
                            a[ip + j * a_dim1] = temp;
                        }
                    }
                } else {
                    ip = -ipiv[i__];
                    --i__;
                    if (i__ > 1) {
                        i__1 = i__ - 1;
                        for (j = 1; j <= i__1; ++j) {
                            temp = a[i__ + 1 + j * a_dim1];
                            a[i__ + 1 + j * a_dim1] = a[ip + j * a_dim1];
                            a[ip + j * a_dim1] = temp;
                        }
                    }
                }
                --i__;
            }
            i__ = 1;
            while (i__ <= *n - 1) {
                if (ipiv[i__] < 0) {
                    a[i__ + 1 + i__ * a_dim1] = e[i__];
                    ++i__;
                }
                ++i__;
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
