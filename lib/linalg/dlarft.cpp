#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static doublereal c_b7 = 1.;
int dlarft_(char *direct, char *storev, integer *n, integer *k, doublereal *v, integer *ldv,
            doublereal *tau, doublereal *t, integer *ldt, ftnlen direct_len, ftnlen storev_len)
{
    integer t_dim1, t_offset, v_dim1, v_offset, i__1, i__2, i__3;
    doublereal d__1;
    integer i__, j, prevlastv;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int dgemv_(char *, integer *, integer *, doublereal *, doublereal *, integer *,
                      doublereal *, integer *, doublereal *, doublereal *, integer *, ftnlen);
    integer lastv;
    extern int dtrmv_(char *, char *, char *, integer *, doublereal *, integer *, doublereal *,
                      integer *, ftnlen, ftnlen, ftnlen);
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --tau;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    if (*n == 0) {
        return 0;
    }
    if (lsame_(direct, (char *)"F", (ftnlen)1, (ftnlen)1)) {
        prevlastv = *n;
        i__1 = *k;
        for (i__ = 1; i__ <= i__1; ++i__) {
            prevlastv = max(i__, prevlastv);
            if (tau[i__] == 0.) {
                i__2 = i__;
                for (j = 1; j <= i__2; ++j) {
                    t[j + i__ * t_dim1] = 0.;
                }
            } else {
                if (lsame_(storev, (char *)"C", (ftnlen)1, (ftnlen)1)) {
                    i__2 = i__ + 1;
                    for (lastv = *n; lastv >= i__2; --lastv) {
                        if (v[lastv + i__ * v_dim1] != 0.) {
                            goto L219;
                        }
                    }
                L219:
                    i__2 = i__ - 1;
                    for (j = 1; j <= i__2; ++j) {
                        t[j + i__ * t_dim1] = -tau[i__] * v[i__ + j * v_dim1];
                    }
                    j = min(lastv, prevlastv);
                    i__2 = j - i__;
                    i__3 = i__ - 1;
                    d__1 = -tau[i__];
                    dgemv_((char *)"Transpose", &i__2, &i__3, &d__1, &v[i__ + 1 + v_dim1], ldv,
                           &v[i__ + 1 + i__ * v_dim1], &c__1, &c_b7, &t[i__ * t_dim1 + 1], &c__1,
                           (ftnlen)9);
                } else {
                    i__2 = i__ + 1;
                    for (lastv = *n; lastv >= i__2; --lastv) {
                        if (v[i__ + lastv * v_dim1] != 0.) {
                            goto L235;
                        }
                    }
                L235:
                    i__2 = i__ - 1;
                    for (j = 1; j <= i__2; ++j) {
                        t[j + i__ * t_dim1] = -tau[i__] * v[j + i__ * v_dim1];
                    }
                    j = min(lastv, prevlastv);
                    i__2 = i__ - 1;
                    i__3 = j - i__;
                    d__1 = -tau[i__];
                    dgemv_((char *)"No transpose", &i__2, &i__3, &d__1, &v[(i__ + 1) * v_dim1 + 1], ldv,
                           &v[i__ + (i__ + 1) * v_dim1], ldv, &c_b7, &t[i__ * t_dim1 + 1], &c__1,
                           (ftnlen)12);
                }
                i__2 = i__ - 1;
                dtrmv_((char *)"Upper", (char *)"No transpose", (char *)"Non-unit", &i__2, &t[t_offset], ldt,
                       &t[i__ * t_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
                t[i__ + i__ * t_dim1] = tau[i__];
                if (i__ > 1) {
                    prevlastv = max(prevlastv, lastv);
                } else {
                    prevlastv = lastv;
                }
            }
        }
    } else {
        prevlastv = 1;
        for (i__ = *k; i__ >= 1; --i__) {
            if (tau[i__] == 0.) {
                i__1 = *k;
                for (j = i__; j <= i__1; ++j) {
                    t[j + i__ * t_dim1] = 0.;
                }
            } else {
                if (i__ < *k) {
                    if (lsame_(storev, (char *)"C", (ftnlen)1, (ftnlen)1)) {
                        i__1 = i__ - 1;
                        for (lastv = 1; lastv <= i__1; ++lastv) {
                            if (v[lastv + i__ * v_dim1] != 0.) {
                                goto L280;
                            }
                        }
                    L280:
                        i__1 = *k;
                        for (j = i__ + 1; j <= i__1; ++j) {
                            t[j + i__ * t_dim1] = -tau[i__] * v[*n - *k + i__ + j * v_dim1];
                        }
                        j = max(lastv, prevlastv);
                        i__1 = *n - *k + i__ - j;
                        i__2 = *k - i__;
                        d__1 = -tau[i__];
                        dgemv_((char *)"Transpose", &i__1, &i__2, &d__1, &v[j + (i__ + 1) * v_dim1], ldv,
                               &v[j + i__ * v_dim1], &c__1, &c_b7, &t[i__ + 1 + i__ * t_dim1],
                               &c__1, (ftnlen)9);
                    } else {
                        i__1 = i__ - 1;
                        for (lastv = 1; lastv <= i__1; ++lastv) {
                            if (v[i__ + lastv * v_dim1] != 0.) {
                                goto L296;
                            }
                        }
                    L296:
                        i__1 = *k;
                        for (j = i__ + 1; j <= i__1; ++j) {
                            t[j + i__ * t_dim1] = -tau[i__] * v[j + (*n - *k + i__) * v_dim1];
                        }
                        j = max(lastv, prevlastv);
                        i__1 = *k - i__;
                        i__2 = *n - *k + i__ - j;
                        d__1 = -tau[i__];
                        dgemv_((char *)"No transpose", &i__1, &i__2, &d__1, &v[i__ + 1 + j * v_dim1], ldv,
                               &v[i__ + j * v_dim1], ldv, &c_b7, &t[i__ + 1 + i__ * t_dim1], &c__1,
                               (ftnlen)12);
                    }
                    i__1 = *k - i__;
                    dtrmv_((char *)"Lower", (char *)"No transpose", (char *)"Non-unit", &i__1,
                           &t[i__ + 1 + (i__ + 1) * t_dim1], ldt, &t[i__ + 1 + i__ * t_dim1], &c__1,
                           (ftnlen)5, (ftnlen)12, (ftnlen)8);
                    if (i__ > 1) {
                        prevlastv = min(prevlastv, lastv);
                    } else {
                        prevlastv = lastv;
                    }
                }
                t[i__ + i__ * t_dim1] = tau[i__];
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
