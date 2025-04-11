#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublecomplex c_b1 = {1., 0.};
static integer c__1 = 1;
int zlarfb_(char *side, char *trans, char *direct, char *storev, integer *m, integer *n, integer *k,
            doublecomplex *v, integer *ldv, doublecomplex *t, integer *ldt, doublecomplex *c__,
            integer *ldc, doublecomplex *work, integer *ldwork, ftnlen side_len, ftnlen trans_len,
            ftnlen direct_len, ftnlen storev_len)
{
    integer c_dim1, c_offset, t_dim1, t_offset, v_dim1, v_offset, work_dim1, work_offset, i__1,
        i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2;
    void d_lmp_cnjg(doublecomplex *, doublecomplex *);
    integer i__, j;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int zgemm_(char *, char *, integer *, integer *, integer *, doublecomplex *,
                      doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *,
                      doublecomplex *, integer *, ftnlen, ftnlen),
        zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, integer *),
        ztrmm_(char *, char *, char *, char *, integer *, integer *, doublecomplex *,
               doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, ftnlen, ftnlen,
               ftnlen),
        zlacgv_(integer *, doublecomplex *, integer *);
    char transt[1];
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    work_dim1 = *ldwork;
    work_offset = 1 + work_dim1;
    work -= work_offset;
    if (*m <= 0 || *n <= 0) {
        return 0;
    }
    if (lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1)) {
        *(unsigned char *)transt = 'C';
    } else {
        *(unsigned char *)transt = 'N';
    }
    if (lsame_(storev, (char *)"C", (ftnlen)1, (ftnlen)1)) {
        if (lsame_(direct, (char *)"F", (ftnlen)1, (ftnlen)1)) {
            if (lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1)) {
                i__1 = *k;
                for (j = 1; j <= i__1; ++j) {
                    zcopy_(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1], &c__1);
                    zlacgv_(n, &work[j * work_dim1 + 1], &c__1);
                }
                ztrmm_((char *)"R", (char *)"L", (char *)"N", (char *)"U", n, k, &c_b1, &v[v_offset], ldv, &work[work_offset],
                       ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                if (*m > *k) {
                    i__1 = *m - *k;
                    zgemm_((char *)"C", (char *)"N", n, k, &i__1, &c_b1, &c__[*k + 1 + c_dim1], ldc,
                           &v[*k + 1 + v_dim1], ldv, &c_b1, &work[work_offset], ldwork, (ftnlen)1,
                           (ftnlen)1);
                }
                ztrmm_((char *)"R", (char *)"U", transt, (char *)"N", n, k, &c_b1, &t[t_offset], ldt, &work[work_offset],
                       ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                if (*m > *k) {
                    i__1 = *m - *k;
                    z__1.r = -1., z__1.i = -0.;
                    zgemm_((char *)"N", (char *)"C", &i__1, n, k, &z__1, &v[*k + 1 + v_dim1], ldv,
                           &work[work_offset], ldwork, &c_b1, &c__[*k + 1 + c_dim1], ldc, (ftnlen)1,
                           (ftnlen)1);
                }
                ztrmm_((char *)"R", (char *)"L", (char *)"C", (char *)"U", n, k, &c_b1, &v[v_offset], ldv, &work[work_offset],
                       ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                i__1 = *k;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *n;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = j + i__ * c_dim1;
                        i__4 = j + i__ * c_dim1;
                        d_lmp_cnjg(&z__2, &work[i__ + j * work_dim1]);
                        z__1.r = c__[i__4].r - z__2.r, z__1.i = c__[i__4].i - z__2.i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                }
            } else if (lsame_(side, (char *)"R", (ftnlen)1, (ftnlen)1)) {
                i__1 = *k;
                for (j = 1; j <= i__1; ++j) {
                    zcopy_(m, &c__[j * c_dim1 + 1], &c__1, &work[j * work_dim1 + 1], &c__1);
                }
                ztrmm_((char *)"R", (char *)"L", (char *)"N", (char *)"U", m, k, &c_b1, &v[v_offset], ldv, &work[work_offset],
                       ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                if (*n > *k) {
                    i__1 = *n - *k;
                    zgemm_((char *)"N", (char *)"N", m, k, &i__1, &c_b1, &c__[(*k + 1) * c_dim1 + 1], ldc,
                           &v[*k + 1 + v_dim1], ldv, &c_b1, &work[work_offset], ldwork, (ftnlen)1,
                           (ftnlen)1);
                }
                ztrmm_((char *)"R", (char *)"U", trans, (char *)"N", m, k, &c_b1, &t[t_offset], ldt, &work[work_offset],
                       ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                if (*n > *k) {
                    i__1 = *n - *k;
                    z__1.r = -1., z__1.i = -0.;
                    zgemm_((char *)"N", (char *)"C", m, &i__1, k, &z__1, &work[work_offset], ldwork,
                           &v[*k + 1 + v_dim1], ldv, &c_b1, &c__[(*k + 1) * c_dim1 + 1], ldc,
                           (ftnlen)1, (ftnlen)1);
                }
                ztrmm_((char *)"R", (char *)"L", (char *)"C", (char *)"U", m, k, &c_b1, &v[v_offset], ldv, &work[work_offset],
                       ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                i__1 = *k;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + j * c_dim1;
                        i__5 = i__ + j * work_dim1;
                        z__1.r = c__[i__4].r - work[i__5].r, z__1.i = c__[i__4].i - work[i__5].i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                }
            }
        } else {
            if (lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1)) {
                i__1 = *k;
                for (j = 1; j <= i__1; ++j) {
                    zcopy_(n, &c__[*m - *k + j + c_dim1], ldc, &work[j * work_dim1 + 1], &c__1);
                    zlacgv_(n, &work[j * work_dim1 + 1], &c__1);
                }
                ztrmm_((char *)"R", (char *)"U", (char *)"N", (char *)"U", n, k, &c_b1, &v[*m - *k + 1 + v_dim1], ldv,
                       &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                if (*m > *k) {
                    i__1 = *m - *k;
                    zgemm_((char *)"C", (char *)"N", n, k, &i__1, &c_b1, &c__[c_offset], ldc, &v[v_offset], ldv,
                           &c_b1, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1);
                }
                ztrmm_((char *)"R", (char *)"L", transt, (char *)"N", n, k, &c_b1, &t[t_offset], ldt, &work[work_offset],
                       ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                if (*m > *k) {
                    i__1 = *m - *k;
                    z__1.r = -1., z__1.i = -0.;
                    zgemm_((char *)"N", (char *)"C", &i__1, n, k, &z__1, &v[v_offset], ldv, &work[work_offset],
                           ldwork, &c_b1, &c__[c_offset], ldc, (ftnlen)1, (ftnlen)1);
                }
                ztrmm_((char *)"R", (char *)"U", (char *)"C", (char *)"U", n, k, &c_b1, &v[*m - *k + 1 + v_dim1], ldv,
                       &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                i__1 = *k;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *n;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = *m - *k + j + i__ * c_dim1;
                        i__4 = *m - *k + j + i__ * c_dim1;
                        d_lmp_cnjg(&z__2, &work[i__ + j * work_dim1]);
                        z__1.r = c__[i__4].r - z__2.r, z__1.i = c__[i__4].i - z__2.i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                }
            } else if (lsame_(side, (char *)"R", (ftnlen)1, (ftnlen)1)) {
                i__1 = *k;
                for (j = 1; j <= i__1; ++j) {
                    zcopy_(m, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &work[j * work_dim1 + 1],
                           &c__1);
                }
                ztrmm_((char *)"R", (char *)"U", (char *)"N", (char *)"U", m, k, &c_b1, &v[*n - *k + 1 + v_dim1], ldv,
                       &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                if (*n > *k) {
                    i__1 = *n - *k;
                    zgemm_((char *)"N", (char *)"N", m, k, &i__1, &c_b1, &c__[c_offset], ldc, &v[v_offset], ldv,
                           &c_b1, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1);
                }
                ztrmm_((char *)"R", (char *)"L", trans, (char *)"N", m, k, &c_b1, &t[t_offset], ldt, &work[work_offset],
                       ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                if (*n > *k) {
                    i__1 = *n - *k;
                    z__1.r = -1., z__1.i = -0.;
                    zgemm_((char *)"N", (char *)"C", m, &i__1, k, &z__1, &work[work_offset], ldwork, &v[v_offset],
                           ldv, &c_b1, &c__[c_offset], ldc, (ftnlen)1, (ftnlen)1);
                }
                ztrmm_((char *)"R", (char *)"U", (char *)"C", (char *)"U", m, k, &c_b1, &v[*n - *k + 1 + v_dim1], ldv,
                       &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                i__1 = *k;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + (*n - *k + j) * c_dim1;
                        i__4 = i__ + (*n - *k + j) * c_dim1;
                        i__5 = i__ + j * work_dim1;
                        z__1.r = c__[i__4].r - work[i__5].r, z__1.i = c__[i__4].i - work[i__5].i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                }
            }
        }
    } else if (lsame_(storev, (char *)"R", (ftnlen)1, (ftnlen)1)) {
        if (lsame_(direct, (char *)"F", (ftnlen)1, (ftnlen)1)) {
            if (lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1)) {
                i__1 = *k;
                for (j = 1; j <= i__1; ++j) {
                    zcopy_(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1], &c__1);
                    zlacgv_(n, &work[j * work_dim1 + 1], &c__1);
                }
                ztrmm_((char *)"R", (char *)"U", (char *)"C", (char *)"U", n, k, &c_b1, &v[v_offset], ldv, &work[work_offset],
                       ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                if (*m > *k) {
                    i__1 = *m - *k;
                    zgemm_((char *)"C", (char *)"C", n, k, &i__1, &c_b1, &c__[*k + 1 + c_dim1], ldc,
                           &v[(*k + 1) * v_dim1 + 1], ldv, &c_b1, &work[work_offset], ldwork,
                           (ftnlen)1, (ftnlen)1);
                }
                ztrmm_((char *)"R", (char *)"U", transt, (char *)"N", n, k, &c_b1, &t[t_offset], ldt, &work[work_offset],
                       ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                if (*m > *k) {
                    i__1 = *m - *k;
                    z__1.r = -1., z__1.i = -0.;
                    zgemm_((char *)"C", (char *)"C", &i__1, n, k, &z__1, &v[(*k + 1) * v_dim1 + 1], ldv,
                           &work[work_offset], ldwork, &c_b1, &c__[*k + 1 + c_dim1], ldc, (ftnlen)1,
                           (ftnlen)1);
                }
                ztrmm_((char *)"R", (char *)"U", (char *)"N", (char *)"U", n, k, &c_b1, &v[v_offset], ldv, &work[work_offset],
                       ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                i__1 = *k;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *n;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = j + i__ * c_dim1;
                        i__4 = j + i__ * c_dim1;
                        d_lmp_cnjg(&z__2, &work[i__ + j * work_dim1]);
                        z__1.r = c__[i__4].r - z__2.r, z__1.i = c__[i__4].i - z__2.i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                }
            } else if (lsame_(side, (char *)"R", (ftnlen)1, (ftnlen)1)) {
                i__1 = *k;
                for (j = 1; j <= i__1; ++j) {
                    zcopy_(m, &c__[j * c_dim1 + 1], &c__1, &work[j * work_dim1 + 1], &c__1);
                }
                ztrmm_((char *)"R", (char *)"U", (char *)"C", (char *)"U", m, k, &c_b1, &v[v_offset], ldv, &work[work_offset],
                       ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                if (*n > *k) {
                    i__1 = *n - *k;
                    zgemm_((char *)"N", (char *)"C", m, k, &i__1, &c_b1, &c__[(*k + 1) * c_dim1 + 1], ldc,
                           &v[(*k + 1) * v_dim1 + 1], ldv, &c_b1, &work[work_offset], ldwork,
                           (ftnlen)1, (ftnlen)1);
                }
                ztrmm_((char *)"R", (char *)"U", trans, (char *)"N", m, k, &c_b1, &t[t_offset], ldt, &work[work_offset],
                       ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                if (*n > *k) {
                    i__1 = *n - *k;
                    z__1.r = -1., z__1.i = -0.;
                    zgemm_((char *)"N", (char *)"N", m, &i__1, k, &z__1, &work[work_offset], ldwork,
                           &v[(*k + 1) * v_dim1 + 1], ldv, &c_b1, &c__[(*k + 1) * c_dim1 + 1], ldc,
                           (ftnlen)1, (ftnlen)1);
                }
                ztrmm_((char *)"R", (char *)"U", (char *)"N", (char *)"U", m, k, &c_b1, &v[v_offset], ldv, &work[work_offset],
                       ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                i__1 = *k;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + j * c_dim1;
                        i__5 = i__ + j * work_dim1;
                        z__1.r = c__[i__4].r - work[i__5].r, z__1.i = c__[i__4].i - work[i__5].i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                }
            }
        } else {
            if (lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1)) {
                i__1 = *k;
                for (j = 1; j <= i__1; ++j) {
                    zcopy_(n, &c__[*m - *k + j + c_dim1], ldc, &work[j * work_dim1 + 1], &c__1);
                    zlacgv_(n, &work[j * work_dim1 + 1], &c__1);
                }
                ztrmm_((char *)"R", (char *)"L", (char *)"C", (char *)"U", n, k, &c_b1, &v[(*m - *k + 1) * v_dim1 + 1], ldv,
                       &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                if (*m > *k) {
                    i__1 = *m - *k;
                    zgemm_((char *)"C", (char *)"C", n, k, &i__1, &c_b1, &c__[c_offset], ldc, &v[v_offset], ldv,
                           &c_b1, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1);
                }
                ztrmm_((char *)"R", (char *)"L", transt, (char *)"N", n, k, &c_b1, &t[t_offset], ldt, &work[work_offset],
                       ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                if (*m > *k) {
                    i__1 = *m - *k;
                    z__1.r = -1., z__1.i = -0.;
                    zgemm_((char *)"C", (char *)"C", &i__1, n, k, &z__1, &v[v_offset], ldv, &work[work_offset],
                           ldwork, &c_b1, &c__[c_offset], ldc, (ftnlen)1, (ftnlen)1);
                }
                ztrmm_((char *)"R", (char *)"L", (char *)"N", (char *)"U", n, k, &c_b1, &v[(*m - *k + 1) * v_dim1 + 1], ldv,
                       &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                i__1 = *k;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *n;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = *m - *k + j + i__ * c_dim1;
                        i__4 = *m - *k + j + i__ * c_dim1;
                        d_lmp_cnjg(&z__2, &work[i__ + j * work_dim1]);
                        z__1.r = c__[i__4].r - z__2.r, z__1.i = c__[i__4].i - z__2.i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
                    }
                }
            } else if (lsame_(side, (char *)"R", (ftnlen)1, (ftnlen)1)) {
                i__1 = *k;
                for (j = 1; j <= i__1; ++j) {
                    zcopy_(m, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &work[j * work_dim1 + 1],
                           &c__1);
                }
                ztrmm_((char *)"R", (char *)"L", (char *)"C", (char *)"U", m, k, &c_b1, &v[(*n - *k + 1) * v_dim1 + 1], ldv,
                       &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                if (*n > *k) {
                    i__1 = *n - *k;
                    zgemm_((char *)"N", (char *)"C", m, k, &i__1, &c_b1, &c__[c_offset], ldc, &v[v_offset], ldv,
                           &c_b1, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1);
                }
                ztrmm_((char *)"R", (char *)"L", trans, (char *)"N", m, k, &c_b1, &t[t_offset], ldt, &work[work_offset],
                       ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                if (*n > *k) {
                    i__1 = *n - *k;
                    z__1.r = -1., z__1.i = -0.;
                    zgemm_((char *)"N", (char *)"N", m, &i__1, k, &z__1, &work[work_offset], ldwork, &v[v_offset],
                           ldv, &c_b1, &c__[c_offset], ldc, (ftnlen)1, (ftnlen)1);
                }
                ztrmm_((char *)"R", (char *)"L", (char *)"N", (char *)"U", m, k, &c_b1, &v[(*n - *k + 1) * v_dim1 + 1], ldv,
                       &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
                i__1 = *k;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + (*n - *k + j) * c_dim1;
                        i__4 = i__ + (*n - *k + j) * c_dim1;
                        i__5 = i__ + j * work_dim1;
                        z__1.r = c__[i__4].r - work[i__5].r, z__1.i = c__[i__4].i - work[i__5].i;
                        c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
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
