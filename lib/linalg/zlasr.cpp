#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int zlasr_(char *side, char *pivot, char *direct, integer *m, integer *n, doublereal *c__,
           doublereal *s, doublecomplex *a, integer *lda, ftnlen side_len, ftnlen pivot_len,
           ftnlen direct_len)
{
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2, z__3;
    integer i__, j, info;
    doublecomplex temp;
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
        xerbla_((char *)"ZLASR ", &info, (ftnlen)6);
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
                            i__3 = j + 1 + i__ * a_dim1;
                            temp.r = a[i__3].r, temp.i = a[i__3].i;
                            i__3 = j + 1 + i__ * a_dim1;
                            z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
                            i__4 = j + i__ * a_dim1;
                            z__3.r = stemp * a[i__4].r, z__3.i = stemp * a[i__4].i;
                            z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
                            a[i__3].r = z__1.r, a[i__3].i = z__1.i;
                            i__3 = j + i__ * a_dim1;
                            z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
                            i__4 = j + i__ * a_dim1;
                            z__3.r = ctemp * a[i__4].r, z__3.i = ctemp * a[i__4].i;
                            z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                            a[i__3].r = z__1.r, a[i__3].i = z__1.i;
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
                            i__2 = j + 1 + i__ * a_dim1;
                            temp.r = a[i__2].r, temp.i = a[i__2].i;
                            i__2 = j + 1 + i__ * a_dim1;
                            z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
                            i__3 = j + i__ * a_dim1;
                            z__3.r = stemp * a[i__3].r, z__3.i = stemp * a[i__3].i;
                            z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
                            a[i__2].r = z__1.r, a[i__2].i = z__1.i;
                            i__2 = j + i__ * a_dim1;
                            z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
                            i__3 = j + i__ * a_dim1;
                            z__3.r = ctemp * a[i__3].r, z__3.i = ctemp * a[i__3].i;
                            z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                            a[i__2].r = z__1.r, a[i__2].i = z__1.i;
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
                            i__3 = j + i__ * a_dim1;
                            temp.r = a[i__3].r, temp.i = a[i__3].i;
                            i__3 = j + i__ * a_dim1;
                            z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
                            i__4 = i__ * a_dim1 + 1;
                            z__3.r = stemp * a[i__4].r, z__3.i = stemp * a[i__4].i;
                            z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
                            a[i__3].r = z__1.r, a[i__3].i = z__1.i;
                            i__3 = i__ * a_dim1 + 1;
                            z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
                            i__4 = i__ * a_dim1 + 1;
                            z__3.r = ctemp * a[i__4].r, z__3.i = ctemp * a[i__4].i;
                            z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                            a[i__3].r = z__1.r, a[i__3].i = z__1.i;
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
                            i__2 = j + i__ * a_dim1;
                            temp.r = a[i__2].r, temp.i = a[i__2].i;
                            i__2 = j + i__ * a_dim1;
                            z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
                            i__3 = i__ * a_dim1 + 1;
                            z__3.r = stemp * a[i__3].r, z__3.i = stemp * a[i__3].i;
                            z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
                            a[i__2].r = z__1.r, a[i__2].i = z__1.i;
                            i__2 = i__ * a_dim1 + 1;
                            z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
                            i__3 = i__ * a_dim1 + 1;
                            z__3.r = ctemp * a[i__3].r, z__3.i = ctemp * a[i__3].i;
                            z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                            a[i__2].r = z__1.r, a[i__2].i = z__1.i;
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
                            i__3 = j + i__ * a_dim1;
                            temp.r = a[i__3].r, temp.i = a[i__3].i;
                            i__3 = j + i__ * a_dim1;
                            i__4 = *m + i__ * a_dim1;
                            z__2.r = stemp * a[i__4].r, z__2.i = stemp * a[i__4].i;
                            z__3.r = ctemp * temp.r, z__3.i = ctemp * temp.i;
                            z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                            a[i__3].r = z__1.r, a[i__3].i = z__1.i;
                            i__3 = *m + i__ * a_dim1;
                            i__4 = *m + i__ * a_dim1;
                            z__2.r = ctemp * a[i__4].r, z__2.i = ctemp * a[i__4].i;
                            z__3.r = stemp * temp.r, z__3.i = stemp * temp.i;
                            z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
                            a[i__3].r = z__1.r, a[i__3].i = z__1.i;
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
                            i__2 = j + i__ * a_dim1;
                            temp.r = a[i__2].r, temp.i = a[i__2].i;
                            i__2 = j + i__ * a_dim1;
                            i__3 = *m + i__ * a_dim1;
                            z__2.r = stemp * a[i__3].r, z__2.i = stemp * a[i__3].i;
                            z__3.r = ctemp * temp.r, z__3.i = ctemp * temp.i;
                            z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                            a[i__2].r = z__1.r, a[i__2].i = z__1.i;
                            i__2 = *m + i__ * a_dim1;
                            i__3 = *m + i__ * a_dim1;
                            z__2.r = ctemp * a[i__3].r, z__2.i = ctemp * a[i__3].i;
                            z__3.r = stemp * temp.r, z__3.i = stemp * temp.i;
                            z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
                            a[i__2].r = z__1.r, a[i__2].i = z__1.i;
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
                            i__3 = i__ + (j + 1) * a_dim1;
                            temp.r = a[i__3].r, temp.i = a[i__3].i;
                            i__3 = i__ + (j + 1) * a_dim1;
                            z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
                            i__4 = i__ + j * a_dim1;
                            z__3.r = stemp * a[i__4].r, z__3.i = stemp * a[i__4].i;
                            z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
                            a[i__3].r = z__1.r, a[i__3].i = z__1.i;
                            i__3 = i__ + j * a_dim1;
                            z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
                            i__4 = i__ + j * a_dim1;
                            z__3.r = ctemp * a[i__4].r, z__3.i = ctemp * a[i__4].i;
                            z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                            a[i__3].r = z__1.r, a[i__3].i = z__1.i;
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
                            i__2 = i__ + (j + 1) * a_dim1;
                            temp.r = a[i__2].r, temp.i = a[i__2].i;
                            i__2 = i__ + (j + 1) * a_dim1;
                            z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
                            i__3 = i__ + j * a_dim1;
                            z__3.r = stemp * a[i__3].r, z__3.i = stemp * a[i__3].i;
                            z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
                            a[i__2].r = z__1.r, a[i__2].i = z__1.i;
                            i__2 = i__ + j * a_dim1;
                            z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
                            i__3 = i__ + j * a_dim1;
                            z__3.r = ctemp * a[i__3].r, z__3.i = ctemp * a[i__3].i;
                            z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                            a[i__2].r = z__1.r, a[i__2].i = z__1.i;
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
                            i__3 = i__ + j * a_dim1;
                            temp.r = a[i__3].r, temp.i = a[i__3].i;
                            i__3 = i__ + j * a_dim1;
                            z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
                            i__4 = i__ + a_dim1;
                            z__3.r = stemp * a[i__4].r, z__3.i = stemp * a[i__4].i;
                            z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
                            a[i__3].r = z__1.r, a[i__3].i = z__1.i;
                            i__3 = i__ + a_dim1;
                            z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
                            i__4 = i__ + a_dim1;
                            z__3.r = ctemp * a[i__4].r, z__3.i = ctemp * a[i__4].i;
                            z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                            a[i__3].r = z__1.r, a[i__3].i = z__1.i;
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
                            i__2 = i__ + j * a_dim1;
                            temp.r = a[i__2].r, temp.i = a[i__2].i;
                            i__2 = i__ + j * a_dim1;
                            z__2.r = ctemp * temp.r, z__2.i = ctemp * temp.i;
                            i__3 = i__ + a_dim1;
                            z__3.r = stemp * a[i__3].r, z__3.i = stemp * a[i__3].i;
                            z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
                            a[i__2].r = z__1.r, a[i__2].i = z__1.i;
                            i__2 = i__ + a_dim1;
                            z__2.r = stemp * temp.r, z__2.i = stemp * temp.i;
                            i__3 = i__ + a_dim1;
                            z__3.r = ctemp * a[i__3].r, z__3.i = ctemp * a[i__3].i;
                            z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                            a[i__2].r = z__1.r, a[i__2].i = z__1.i;
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
                            i__3 = i__ + j * a_dim1;
                            temp.r = a[i__3].r, temp.i = a[i__3].i;
                            i__3 = i__ + j * a_dim1;
                            i__4 = i__ + *n * a_dim1;
                            z__2.r = stemp * a[i__4].r, z__2.i = stemp * a[i__4].i;
                            z__3.r = ctemp * temp.r, z__3.i = ctemp * temp.i;
                            z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                            a[i__3].r = z__1.r, a[i__3].i = z__1.i;
                            i__3 = i__ + *n * a_dim1;
                            i__4 = i__ + *n * a_dim1;
                            z__2.r = ctemp * a[i__4].r, z__2.i = ctemp * a[i__4].i;
                            z__3.r = stemp * temp.r, z__3.i = stemp * temp.i;
                            z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
                            a[i__3].r = z__1.r, a[i__3].i = z__1.i;
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
                            i__2 = i__ + j * a_dim1;
                            temp.r = a[i__2].r, temp.i = a[i__2].i;
                            i__2 = i__ + j * a_dim1;
                            i__3 = i__ + *n * a_dim1;
                            z__2.r = stemp * a[i__3].r, z__2.i = stemp * a[i__3].i;
                            z__3.r = ctemp * temp.r, z__3.i = ctemp * temp.i;
                            z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
                            a[i__2].r = z__1.r, a[i__2].i = z__1.i;
                            i__2 = i__ + *n * a_dim1;
                            i__3 = i__ + *n * a_dim1;
                            z__2.r = ctemp * a[i__3].r, z__2.i = ctemp * a[i__3].i;
                            z__3.r = stemp * temp.r, z__3.i = stemp * temp.i;
                            z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
                            a[i__2].r = z__1.r, a[i__2].i = z__1.i;
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
