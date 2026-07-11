#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublecomplex c_b1 = {1., 0.};
static integer c__1 = 1;
int zlarft_(char *direct, char *storev, integer *n, integer *k, doublecomplex *v, integer *ldv,
            doublecomplex *tau, doublecomplex *t, integer *ldt, ftnlen direct_len,
            ftnlen storev_len)
{
    integer t_dim1, t_offset, v_dim1, v_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3;
    void d_lmp_cnjg(doublecomplex *, doublecomplex *);
    integer i__, j, prevlastv;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int zgemm_(char *, char *, integer *, integer *, integer *, doublecomplex *,
                      doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *,
                      doublecomplex *, integer *, ftnlen, ftnlen),
        zgemv_(char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *,
               doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    integer lastv;
    extern int ztrmv_(char *, char *, char *, integer *, doublecomplex *, integer *,
                      doublecomplex *, integer *, ftnlen, ftnlen, ftnlen);
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
            prevlastv = max(prevlastv, i__);
            i__2 = i__;
            if (tau[i__2].r == 0. && tau[i__2].i == 0.) {
                i__2 = i__;
                for (j = 1; j <= i__2; ++j) {
                    i__3 = j + i__ * t_dim1;
                    t[i__3].r = 0., t[i__3].i = 0.;
                }
            } else {
                if (lsame_(storev, (char *)"C", (ftnlen)1, (ftnlen)1)) {
                    i__2 = i__ + 1;
                    for (lastv = *n; lastv >= i__2; --lastv) {
                        i__3 = lastv + i__ * v_dim1;
                        if (v[i__3].r != 0. || v[i__3].i != 0.) {
                            goto L220;
                        }
                    }
                L220:
                    i__2 = i__ - 1;
                    for (j = 1; j <= i__2; ++j) {
                        i__3 = j + i__ * t_dim1;
                        i__4 = i__;
                        z__2.r = -tau[i__4].r, z__2.i = -tau[i__4].i;
                        d_lmp_cnjg(&z__3, &v[i__ + j * v_dim1]);
                        z__1.r = z__2.r * z__3.r - z__2.i * z__3.i,
                        z__1.i = z__2.r * z__3.i + z__2.i * z__3.r;
                        t[i__3].r = z__1.r, t[i__3].i = z__1.i;
                    }
                    j = min(lastv, prevlastv);
                    i__2 = j - i__;
                    i__3 = i__ - 1;
                    i__4 = i__;
                    z__1.r = -tau[i__4].r, z__1.i = -tau[i__4].i;
                    zgemv_((char *)"Conjugate transpose", &i__2, &i__3, &z__1, &v[i__ + 1 + v_dim1], ldv,
                           &v[i__ + 1 + i__ * v_dim1], &c__1, &c_b1, &t[i__ * t_dim1 + 1], &c__1,
                           (ftnlen)19);
                } else {
                    i__2 = i__ + 1;
                    for (lastv = *n; lastv >= i__2; --lastv) {
                        i__3 = i__ + lastv * v_dim1;
                        if (v[i__3].r != 0. || v[i__3].i != 0.) {
                            goto L236;
                        }
                    }
                L236:
                    i__2 = i__ - 1;
                    for (j = 1; j <= i__2; ++j) {
                        i__3 = j + i__ * t_dim1;
                        i__4 = i__;
                        z__2.r = -tau[i__4].r, z__2.i = -tau[i__4].i;
                        i__5 = j + i__ * v_dim1;
                        z__1.r = z__2.r * v[i__5].r - z__2.i * v[i__5].i,
                        z__1.i = z__2.r * v[i__5].i + z__2.i * v[i__5].r;
                        t[i__3].r = z__1.r, t[i__3].i = z__1.i;
                    }
                    j = min(lastv, prevlastv);
                    i__2 = i__ - 1;
                    i__3 = j - i__;
                    i__4 = i__;
                    z__1.r = -tau[i__4].r, z__1.i = -tau[i__4].i;
                    zgemm_((char *)"N", (char *)"C", &i__2, &c__1, &i__3, &z__1, &v[(i__ + 1) * v_dim1 + 1], ldv,
                           &v[i__ + (i__ + 1) * v_dim1], ldv, &c_b1, &t[i__ * t_dim1 + 1], ldt,
                           (ftnlen)1, (ftnlen)1);
                }
                i__2 = i__ - 1;
                ztrmv_((char *)"Upper", (char *)"No transpose", (char *)"Non-unit", &i__2, &t[t_offset], ldt,
                       &t[i__ * t_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
                i__2 = i__ + i__ * t_dim1;
                i__3 = i__;
                t[i__2].r = tau[i__3].r, t[i__2].i = tau[i__3].i;
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
            i__1 = i__;
            if (tau[i__1].r == 0. && tau[i__1].i == 0.) {
                i__1 = *k;
                for (j = i__; j <= i__1; ++j) {
                    i__2 = j + i__ * t_dim1;
                    t[i__2].r = 0., t[i__2].i = 0.;
                }
            } else {
                if (i__ < *k) {
                    if (lsame_(storev, (char *)"C", (ftnlen)1, (ftnlen)1)) {
                        i__1 = i__ - 1;
                        for (lastv = 1; lastv <= i__1; ++lastv) {
                            i__2 = lastv + i__ * v_dim1;
                            if (v[i__2].r != 0. || v[i__2].i != 0.) {
                                goto L281;
                            }
                        }
                    L281:
                        i__1 = *k;
                        for (j = i__ + 1; j <= i__1; ++j) {
                            i__2 = j + i__ * t_dim1;
                            i__3 = i__;
                            z__2.r = -tau[i__3].r, z__2.i = -tau[i__3].i;
                            d_lmp_cnjg(&z__3, &v[*n - *k + i__ + j * v_dim1]);
                            z__1.r = z__2.r * z__3.r - z__2.i * z__3.i,
                            z__1.i = z__2.r * z__3.i + z__2.i * z__3.r;
                            t[i__2].r = z__1.r, t[i__2].i = z__1.i;
                        }
                        j = max(lastv, prevlastv);
                        i__1 = *n - *k + i__ - j;
                        i__2 = *k - i__;
                        i__3 = i__;
                        z__1.r = -tau[i__3].r, z__1.i = -tau[i__3].i;
                        zgemv_((char *)"Conjugate transpose", &i__1, &i__2, &z__1,
                               &v[j + (i__ + 1) * v_dim1], ldv, &v[j + i__ * v_dim1], &c__1, &c_b1,
                               &t[i__ + 1 + i__ * t_dim1], &c__1, (ftnlen)19);
                    } else {
                        i__1 = i__ - 1;
                        for (lastv = 1; lastv <= i__1; ++lastv) {
                            i__2 = i__ + lastv * v_dim1;
                            if (v[i__2].r != 0. || v[i__2].i != 0.) {
                                goto L297;
                            }
                        }
                    L297:
                        i__1 = *k;
                        for (j = i__ + 1; j <= i__1; ++j) {
                            i__2 = j + i__ * t_dim1;
                            i__3 = i__;
                            z__2.r = -tau[i__3].r, z__2.i = -tau[i__3].i;
                            i__4 = j + (*n - *k + i__) * v_dim1;
                            z__1.r = z__2.r * v[i__4].r - z__2.i * v[i__4].i,
                            z__1.i = z__2.r * v[i__4].i + z__2.i * v[i__4].r;
                            t[i__2].r = z__1.r, t[i__2].i = z__1.i;
                        }
                        j = max(lastv, prevlastv);
                        i__1 = *k - i__;
                        i__2 = *n - *k + i__ - j;
                        i__3 = i__;
                        z__1.r = -tau[i__3].r, z__1.i = -tau[i__3].i;
                        zgemm_((char *)"N", (char *)"C", &i__1, &c__1, &i__2, &z__1, &v[i__ + 1 + j * v_dim1], ldv,
                               &v[i__ + j * v_dim1], ldv, &c_b1, &t[i__ + 1 + i__ * t_dim1], ldt,
                               (ftnlen)1, (ftnlen)1);
                    }
                    i__1 = *k - i__;
                    ztrmv_((char *)"Lower", (char *)"No transpose", (char *)"Non-unit", &i__1,
                           &t[i__ + 1 + (i__ + 1) * t_dim1], ldt, &t[i__ + 1 + i__ * t_dim1], &c__1,
                           (ftnlen)5, (ftnlen)12, (ftnlen)8);
                    if (i__ > 1) {
                        prevlastv = min(prevlastv, lastv);
                    } else {
                        prevlastv = lastv;
                    }
                }
                i__1 = i__ + i__ * t_dim1;
                i__2 = i__;
                t[i__1].r = tau[i__2].r, t[i__1].i = tau[i__2].i;
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
