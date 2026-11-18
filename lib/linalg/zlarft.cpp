#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublecomplex c_b8 = {1., 0.};
static doublecomplex c_b17 = {-1., 0.};
int zlarft_(char *direct, char *storev, integer *n, integer *k, doublecomplex *v, integer *ldv,
            doublecomplex *tau, doublecomplex *t, integer *ldt, ftnlen direct_len,
            ftnlen storev_len)
{
    integer t_dim1, t_offset, v_dim1, v_offset, i__1, i__2, i__3;
    doublecomplex z__1;
    void d_lmp_cnjg(doublecomplex *, doublecomplex *);
    integer i__, j, l;
    logical lq, ql, qr, dirf, colv;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int zgemm_(char *, char *, integer *, integer *, integer *, doublecomplex *,
                      doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *,
                      doublecomplex *, integer *, ftnlen, ftnlen),
        ztrmm_(char *, char *, char *, char *, integer *, integer *, doublecomplex *,
               doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, ftnlen, ftnlen,
               ftnlen),
        zlacpy_(char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *,
                integer *, ftnlen);
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --tau;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    if (*n == 0 || *k == 0) {
        return 0;
    }
    if (*n == 1 || *k == 1) {
        i__1 = t_dim1 + 1;
        t[i__1].r = tau[1].r, t[i__1].i = tau[1].i;
        return 0;
    }
    l = *k / 2;
    dirf = lsame_(direct, (char *)"F", (ftnlen)1, (ftnlen)1);
    colv = lsame_(storev, (char *)"C", (ftnlen)1, (ftnlen)1);
    qr = dirf && colv;
    lq = dirf && !colv;
    ql = !dirf && colv;
    if (qr) {
        zlarft_(direct, storev, n, &l, &v[v_offset], ldv, &tau[1], &t[t_offset], ldt, (ftnlen)1,
                (ftnlen)1);
        i__1 = *n - l;
        i__2 = *k - l;
        zlarft_(direct, storev, &i__1, &i__2, &v[l + 1 + (l + 1) * v_dim1], ldv, &tau[l + 1],
                &t[l + 1 + (l + 1) * t_dim1], ldt, (ftnlen)1, (ftnlen)1);
        i__1 = l;
        for (j = 1; j <= i__1; ++j) {
            i__2 = *k - l;
            for (i__ = 1; i__ <= i__2; ++i__) {
                i__3 = j + (l + i__) * t_dim1;
                d_lmp_cnjg(&z__1, &v[l + i__ + j * v_dim1]);
                t[i__3].r = z__1.r, t[i__3].i = z__1.i;
            }
        }
        i__1 = *k - l;
        ztrmm_((char *)"R", (char *)"L", (char *)"N", (char *)"U", &l, &i__1, &c_b8, &v[l + 1 + (l + 1) * v_dim1], ldv,
               &t[(l + 1) * t_dim1 + 1], ldt, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
        i__1 = *k - l;
        i__2 = *n - *k;
        zgemm_((char *)"C", (char *)"N", &l, &i__1, &i__2, &c_b8, &v[*k + 1 + v_dim1], ldv,
               &v[*k + 1 + (l + 1) * v_dim1], ldv, &c_b8, &t[(l + 1) * t_dim1 + 1], ldt, (ftnlen)1,
               (ftnlen)1);
        i__1 = *k - l;
        ztrmm_((char *)"L", (char *)"U", (char *)"N", (char *)"N", &l, &i__1, &c_b17, &t[t_offset], ldt, &t[(l + 1) * t_dim1 + 1],
               ldt, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
        i__1 = *k - l;
        ztrmm_((char *)"R", (char *)"U", (char *)"N", (char *)"N", &l, &i__1, &c_b8, &t[l + 1 + (l + 1) * t_dim1], ldt,
               &t[(l + 1) * t_dim1 + 1], ldt, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
    } else if (lq) {
        zlarft_(direct, storev, n, &l, &v[v_offset], ldv, &tau[1], &t[t_offset], ldt, (ftnlen)1,
                (ftnlen)1);
        i__1 = *n - l;
        i__2 = *k - l;
        zlarft_(direct, storev, &i__1, &i__2, &v[l + 1 + (l + 1) * v_dim1], ldv, &tau[l + 1],
                &t[l + 1 + (l + 1) * t_dim1], ldt, (ftnlen)1, (ftnlen)1);
        i__1 = *k - l;
        zlacpy_((char *)"A", &l, &i__1, &v[(l + 1) * v_dim1 + 1], ldv, &t[(l + 1) * t_dim1 + 1], ldt,
                (ftnlen)1);
        i__1 = *k - l;
        ztrmm_((char *)"R", (char *)"U", (char *)"C", (char *)"U", &l, &i__1, &c_b8, &v[l + 1 + (l + 1) * v_dim1], ldv,
               &t[(l + 1) * t_dim1 + 1], ldt, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
        i__1 = *k - l;
        i__2 = *n - *k;
        zgemm_((char *)"N", (char *)"C", &l, &i__1, &i__2, &c_b8, &v[(*k + 1) * v_dim1 + 1], ldv,
               &v[l + 1 + (*k + 1) * v_dim1], ldv, &c_b8, &t[(l + 1) * t_dim1 + 1], ldt, (ftnlen)1,
               (ftnlen)1);
        i__1 = *k - l;
        ztrmm_((char *)"L", (char *)"U", (char *)"N", (char *)"N", &l, &i__1, &c_b17, &t[t_offset], ldt, &t[(l + 1) * t_dim1 + 1],
               ldt, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
        i__1 = *k - l;
        ztrmm_((char *)"R", (char *)"U", (char *)"N", (char *)"N", &l, &i__1, &c_b8, &t[l + 1 + (l + 1) * t_dim1], ldt,
               &t[(l + 1) * t_dim1 + 1], ldt, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
    } else if (ql) {
        i__1 = *n - l;
        i__2 = *k - l;
        zlarft_(direct, storev, &i__1, &i__2, &v[v_offset], ldv, &tau[1], &t[t_offset], ldt,
                (ftnlen)1, (ftnlen)1);
        zlarft_(direct, storev, n, &l, &v[(*k - l + 1) * v_dim1 + 1], ldv, &tau[*k - l + 1],
                &t[*k - l + 1 + (*k - l + 1) * t_dim1], ldt, (ftnlen)1, (ftnlen)1);
        i__1 = *k - l;
        for (j = 1; j <= i__1; ++j) {
            i__2 = l;
            for (i__ = 1; i__ <= i__2; ++i__) {
                i__3 = *k - l + i__ + j * t_dim1;
                d_lmp_cnjg(&z__1, &v[*n - *k + j + (*k - l + i__) * v_dim1]);
                t[i__3].r = z__1.r, t[i__3].i = z__1.i;
            }
        }
        i__1 = *k - l;
        ztrmm_((char *)"R", (char *)"U", (char *)"N", (char *)"U", &l, &i__1, &c_b8, &v[*n - *k + 1 + v_dim1], ldv,
               &t[*k - l + 1 + t_dim1], ldt, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
        i__1 = *k - l;
        i__2 = *n - *k;
        zgemm_((char *)"C", (char *)"N", &l, &i__1, &i__2, &c_b8, &v[(*k - l + 1) * v_dim1 + 1], ldv, &v[v_offset],
               ldv, &c_b8, &t[*k - l + 1 + t_dim1], ldt, (ftnlen)1, (ftnlen)1);
        i__1 = *k - l;
        ztrmm_((char *)"L", (char *)"L", (char *)"N", (char *)"N", &l, &i__1, &c_b17, &t[*k - l + 1 + (*k - l + 1) * t_dim1], ldt,
               &t[*k - l + 1 + t_dim1], ldt, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
        i__1 = *k - l;
        ztrmm_((char *)"R", (char *)"L", (char *)"N", (char *)"N", &l, &i__1, &c_b8, &t[t_offset], ldt, &t[*k - l + 1 + t_dim1],
               ldt, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
    } else {
        i__1 = *n - l;
        i__2 = *k - l;
        zlarft_(direct, storev, &i__1, &i__2, &v[v_offset], ldv, &tau[1], &t[t_offset], ldt,
                (ftnlen)1, (ftnlen)1);
        zlarft_(direct, storev, n, &l, &v[*k - l + 1 + v_dim1], ldv, &tau[*k - l + 1],
                &t[*k - l + 1 + (*k - l + 1) * t_dim1], ldt, (ftnlen)1, (ftnlen)1);
        i__1 = *k - l;
        zlacpy_((char *)"A", &l, &i__1, &v[*k - l + 1 + (*n - *k + 1) * v_dim1], ldv,
                &t[*k - l + 1 + t_dim1], ldt, (ftnlen)1);
        i__1 = *k - l;
        ztrmm_((char *)"R", (char *)"L", (char *)"C", (char *)"U", &l, &i__1, &c_b8, &v[(*n - *k + 1) * v_dim1 + 1], ldv,
               &t[*k - l + 1 + t_dim1], ldt, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
        i__1 = *k - l;
        i__2 = *n - *k;
        zgemm_((char *)"N", (char *)"C", &l, &i__1, &i__2, &c_b8, &v[*k - l + 1 + v_dim1], ldv, &v[v_offset], ldv,
               &c_b8, &t[*k - l + 1 + t_dim1], ldt, (ftnlen)1, (ftnlen)1);
        i__1 = *k - l;
        ztrmm_((char *)"L", (char *)"L", (char *)"No tranpose", (char *)"N", &l, &i__1, &c_b17,
               &t[*k - l + 1 + (*k - l + 1) * t_dim1], ldt, &t[*k - l + 1 + t_dim1], ldt, (ftnlen)1,
               (ftnlen)1, (ftnlen)11, (ftnlen)1);
        i__1 = *k - l;
        ztrmm_((char *)"R", (char *)"L", (char *)"No tranpose", (char *)"N", &l, &i__1, &c_b8, &t[t_offset], ldt,
               &t[*k - l + 1 + t_dim1], ldt, (ftnlen)1, (ftnlen)1, (ftnlen)11, (ftnlen)1);
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
