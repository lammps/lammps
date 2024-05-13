#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlasr_(char *side, char *pivot, char *direct, integer *m, integer *n, doublereal *c__,
           doublereal *s, doublereal *a, integer *lda, ftnlen side_len, ftnlen pivot_len,
           ftnlen direct_len)
{
    integer a_dim1, a_offset, i__1, i__2;
    integer i__, j, info;
    doublereal temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    doublereal ctemp, stemp;
    extern int xerbla_(char *, integer *, ftnlen);
    --c__;
    --s;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    info = 0;
    if (!(lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1) || lsame_(side, (char *)"R", (ftnlen)1, (ftnlen)1))) {
        info = 1;
    } else if (!(lsame_(pivot, (char *)"V", (ftnlen)1, (ftnlen)1) ||
                 lsame_(pivot, (char *)"T", (ftnlen)1, (ftnlen)1) ||
                 lsame_(pivot, (char *)"B", (ftnlen)1, (ftnlen)1))) {
        info = 2;
    } else if (!(lsame_(direct, (char *)"F", (ftnlen)1, (ftnlen)1) ||
                 lsame_(direct, (char *)"B", (ftnlen)1, (ftnlen)1))) {
        info = 3;
    } else if (*m < 0) {
        info = 4;
    } else if (*n < 0) {
        info = 5;
    } else if (*lda < max(1, *m)) {
        info = 9;
    }
    if (info != 0) {
        xerbla_((char *)"DLASR ", &info, (ftnlen)6);
        return 0;
    }
    if (*m == 0 || *n == 0) {
        return 0;
    }
    if (lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        if (lsame_(pivot, (char *)"V", (ftnlen)1, (ftnlen)1)) {
            if (lsame_(direct, (char *)"F", (ftnlen)1, (ftnlen)1)) {
                i__1 = *m - 1;
                for (j = 1; j <= i__1; ++j) {
                    ctemp = c__[j];
                    stemp = s[j];
                    if (ctemp != 1. || stemp != 0.) {
                        i__2 = *n;
                        for (i__ = 1; i__ <= i__2; ++i__) {
                            temp = a[j + 1 + i__ * a_dim1];
                            a[j + 1 + i__ * a_dim1] = ctemp * temp - stemp * a[j + i__ * a_dim1];
                            a[j + i__ * a_dim1] = stemp * temp + ctemp * a[j + i__ * a_dim1];
                        }
                    }
                }
            } else if (lsame_(direct, (char *)"B", (ftnlen)1, (ftnlen)1)) {
                for (j = *m - 1; j >= 1; --j) {
                    ctemp = c__[j];
                    stemp = s[j];
                    if (ctemp != 1. || stemp != 0.) {
                        i__1 = *n;
                        for (i__ = 1; i__ <= i__1; ++i__) {
                            temp = a[j + 1 + i__ * a_dim1];
                            a[j + 1 + i__ * a_dim1] = ctemp * temp - stemp * a[j + i__ * a_dim1];
                            a[j + i__ * a_dim1] = stemp * temp + ctemp * a[j + i__ * a_dim1];
                        }
                    }
                }
            }
        } else if (lsame_(pivot, (char *)"T", (ftnlen)1, (ftnlen)1)) {
            if (lsame_(direct, (char *)"F", (ftnlen)1, (ftnlen)1)) {
                i__1 = *m;
                for (j = 2; j <= i__1; ++j) {
                    ctemp = c__[j - 1];
                    stemp = s[j - 1];
                    if (ctemp != 1. || stemp != 0.) {
                        i__2 = *n;
                        for (i__ = 1; i__ <= i__2; ++i__) {
                            temp = a[j + i__ * a_dim1];
                            a[j + i__ * a_dim1] = ctemp * temp - stemp * a[i__ * a_dim1 + 1];
                            a[i__ * a_dim1 + 1] = stemp * temp + ctemp * a[i__ * a_dim1 + 1];
                        }
                    }
                }
            } else if (lsame_(direct, (char *)"B", (ftnlen)1, (ftnlen)1)) {
                for (j = *m; j >= 2; --j) {
                    ctemp = c__[j - 1];
                    stemp = s[j - 1];
                    if (ctemp != 1. || stemp != 0.) {
                        i__1 = *n;
                        for (i__ = 1; i__ <= i__1; ++i__) {
                            temp = a[j + i__ * a_dim1];
                            a[j + i__ * a_dim1] = ctemp * temp - stemp * a[i__ * a_dim1 + 1];
                            a[i__ * a_dim1 + 1] = stemp * temp + ctemp * a[i__ * a_dim1 + 1];
                        }
                    }
                }
            }
        } else if (lsame_(pivot, (char *)"B", (ftnlen)1, (ftnlen)1)) {
            if (lsame_(direct, (char *)"F", (ftnlen)1, (ftnlen)1)) {
                i__1 = *m - 1;
                for (j = 1; j <= i__1; ++j) {
                    ctemp = c__[j];
                    stemp = s[j];
                    if (ctemp != 1. || stemp != 0.) {
                        i__2 = *n;
                        for (i__ = 1; i__ <= i__2; ++i__) {
                            temp = a[j + i__ * a_dim1];
                            a[j + i__ * a_dim1] = stemp * a[*m + i__ * a_dim1] + ctemp * temp;
                            a[*m + i__ * a_dim1] = ctemp * a[*m + i__ * a_dim1] - stemp * temp;
                        }
                    }
                }
            } else if (lsame_(direct, (char *)"B", (ftnlen)1, (ftnlen)1)) {
                for (j = *m - 1; j >= 1; --j) {
                    ctemp = c__[j];
                    stemp = s[j];
                    if (ctemp != 1. || stemp != 0.) {
                        i__1 = *n;
                        for (i__ = 1; i__ <= i__1; ++i__) {
                            temp = a[j + i__ * a_dim1];
                            a[j + i__ * a_dim1] = stemp * a[*m + i__ * a_dim1] + ctemp * temp;
                            a[*m + i__ * a_dim1] = ctemp * a[*m + i__ * a_dim1] - stemp * temp;
                        }
                    }
                }
            }
        }
    } else if (lsame_(side, (char *)"R", (ftnlen)1, (ftnlen)1)) {
        if (lsame_(pivot, (char *)"V", (ftnlen)1, (ftnlen)1)) {
            if (lsame_(direct, (char *)"F", (ftnlen)1, (ftnlen)1)) {
                i__1 = *n - 1;
                for (j = 1; j <= i__1; ++j) {
                    ctemp = c__[j];
                    stemp = s[j];
                    if (ctemp != 1. || stemp != 0.) {
                        i__2 = *m;
                        for (i__ = 1; i__ <= i__2; ++i__) {
                            temp = a[i__ + (j + 1) * a_dim1];
                            a[i__ + (j + 1) * a_dim1] = ctemp * temp - stemp * a[i__ + j * a_dim1];
                            a[i__ + j * a_dim1] = stemp * temp + ctemp * a[i__ + j * a_dim1];
                        }
                    }
                }
            } else if (lsame_(direct, (char *)"B", (ftnlen)1, (ftnlen)1)) {
                for (j = *n - 1; j >= 1; --j) {
                    ctemp = c__[j];
                    stemp = s[j];
                    if (ctemp != 1. || stemp != 0.) {
                        i__1 = *m;
                        for (i__ = 1; i__ <= i__1; ++i__) {
                            temp = a[i__ + (j + 1) * a_dim1];
                            a[i__ + (j + 1) * a_dim1] = ctemp * temp - stemp * a[i__ + j * a_dim1];
                            a[i__ + j * a_dim1] = stemp * temp + ctemp * a[i__ + j * a_dim1];
                        }
                    }
                }
            }
        } else if (lsame_(pivot, (char *)"T", (ftnlen)1, (ftnlen)1)) {
            if (lsame_(direct, (char *)"F", (ftnlen)1, (ftnlen)1)) {
                i__1 = *n;
                for (j = 2; j <= i__1; ++j) {
                    ctemp = c__[j - 1];
                    stemp = s[j - 1];
                    if (ctemp != 1. || stemp != 0.) {
                        i__2 = *m;
                        for (i__ = 1; i__ <= i__2; ++i__) {
                            temp = a[i__ + j * a_dim1];
                            a[i__ + j * a_dim1] = ctemp * temp - stemp * a[i__ + a_dim1];
                            a[i__ + a_dim1] = stemp * temp + ctemp * a[i__ + a_dim1];
                        }
                    }
                }
            } else if (lsame_(direct, (char *)"B", (ftnlen)1, (ftnlen)1)) {
                for (j = *n; j >= 2; --j) {
                    ctemp = c__[j - 1];
                    stemp = s[j - 1];
                    if (ctemp != 1. || stemp != 0.) {
                        i__1 = *m;
                        for (i__ = 1; i__ <= i__1; ++i__) {
                            temp = a[i__ + j * a_dim1];
                            a[i__ + j * a_dim1] = ctemp * temp - stemp * a[i__ + a_dim1];
                            a[i__ + a_dim1] = stemp * temp + ctemp * a[i__ + a_dim1];
                        }
                    }
                }
            }
        } else if (lsame_(pivot, (char *)"B", (ftnlen)1, (ftnlen)1)) {
            if (lsame_(direct, (char *)"F", (ftnlen)1, (ftnlen)1)) {
                i__1 = *n - 1;
                for (j = 1; j <= i__1; ++j) {
                    ctemp = c__[j];
                    stemp = s[j];
                    if (ctemp != 1. || stemp != 0.) {
                        i__2 = *m;
                        for (i__ = 1; i__ <= i__2; ++i__) {
                            temp = a[i__ + j * a_dim1];
                            a[i__ + j * a_dim1] = stemp * a[i__ + *n * a_dim1] + ctemp * temp;
                            a[i__ + *n * a_dim1] = ctemp * a[i__ + *n * a_dim1] - stemp * temp;
                        }
                    }
                }
            } else if (lsame_(direct, (char *)"B", (ftnlen)1, (ftnlen)1)) {
                for (j = *n - 1; j >= 1; --j) {
                    ctemp = c__[j];
                    stemp = s[j];
                    if (ctemp != 1. || stemp != 0.) {
                        i__1 = *m;
                        for (i__ = 1; i__ <= i__1; ++i__) {
                            temp = a[i__ + j * a_dim1];
                            a[i__ + j * a_dim1] = stemp * a[i__ + *n * a_dim1] + ctemp * temp;
                            a[i__ + *n * a_dim1] = ctemp * a[i__ + *n * a_dim1] - stemp * temp;
                        }
                    }
                }
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
