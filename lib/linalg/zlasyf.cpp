#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublecomplex c_b1 = {1., 0.};
static integer c__1 = 1;
int zlasyf_(char *uplo, integer *n, integer *nb, integer *kb, doublecomplex *a, integer *lda,
            integer *ipiv, doublecomplex *w, integer *ldw, integer *info, ftnlen uplo_len)
{
    integer a_dim1, a_offset, w_dim1, w_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3;
    double sqrt(doublereal), d_lmp_imag(doublecomplex *);
    void z_lmp_div(doublecomplex *, doublecomplex *, doublecomplex *);
    integer j, k;
    doublecomplex t, r1, d11, d21, d22;
    integer jj, kk, jp, kp, kw, kkw, imax, jmax;
    doublereal alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int zscal_(integer *, doublecomplex *, doublecomplex *, integer *);
    integer kstep;
    extern int zgemv_(char *, integer *, integer *, doublecomplex *, doublecomplex *, integer *,
                      doublecomplex *, integer *, doublecomplex *, doublecomplex *, integer *,
                      ftnlen),
        zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, integer *),
        zswap_(integer *, doublecomplex *, integer *, doublecomplex *, integer *);
    doublereal absakk, colmax;
    extern integer izamax_(integer *, doublecomplex *, integer *);
    doublereal rowmax;
    extern int zgemmtr_(char *, char *, char *, integer *, integer *, doublecomplex *,
                        doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *,
                        doublecomplex *, integer *, ftnlen, ftnlen, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    w_dim1 = *ldw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    *info = 0;
    alpha = (sqrt(17.) + 1.) / 8.;
    if (lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1)) {
        k = *n;
    L10:
        kw = *nb + k - *n;
        if (k <= *n - *nb + 1 && *nb < *n || k < 1) {
            goto L30;
        }
        zcopy_(&k, &a[k * a_dim1 + 1], &c__1, &w[kw * w_dim1 + 1], &c__1);
        if (k < *n) {
            i__1 = *n - k;
            z__1.r = -1., z__1.i = -0.;
            zgemv_((char *)"N", &k, &i__1, &z__1, &a[(k + 1) * a_dim1 + 1], lda, &w[k + (kw + 1) * w_dim1],
                   ldw, &c_b1, &w[kw * w_dim1 + 1], &c__1, (ftnlen)1);
        }
        kstep = 1;
        i__1 = k + kw * w_dim1;
        absakk = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_lmp_imag(&w[k + kw * w_dim1]), abs(d__2));
        if (k > 1) {
            i__1 = k - 1;
            imax = izamax_(&i__1, &w[kw * w_dim1 + 1], &c__1);
            i__1 = imax + kw * w_dim1;
            colmax =
                (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_lmp_imag(&w[imax + kw * w_dim1]), abs(d__2));
        } else {
            colmax = 0.;
        }
        if (max(absakk, colmax) == 0.) {
            if (*info == 0) {
                *info = k;
            }
            kp = k;
        } else {
            if (absakk >= alpha * colmax) {
                kp = k;
            } else {
                zcopy_(&imax, &a[imax * a_dim1 + 1], &c__1, &w[(kw - 1) * w_dim1 + 1], &c__1);
                i__1 = k - imax;
                zcopy_(&i__1, &a[imax + (imax + 1) * a_dim1], lda, &w[imax + 1 + (kw - 1) * w_dim1],
                       &c__1);
                if (k < *n) {
                    i__1 = *n - k;
                    z__1.r = -1., z__1.i = -0.;
                    zgemv_((char *)"N", &k, &i__1, &z__1, &a[(k + 1) * a_dim1 + 1], lda,
                           &w[imax + (kw + 1) * w_dim1], ldw, &c_b1, &w[(kw - 1) * w_dim1 + 1],
                           &c__1, (ftnlen)1);
                }
                i__1 = k - imax;
                jmax = imax + izamax_(&i__1, &w[imax + 1 + (kw - 1) * w_dim1], &c__1);
                i__1 = jmax + (kw - 1) * w_dim1;
                rowmax = (d__1 = w[i__1].r, abs(d__1)) +
                         (d__2 = d_lmp_imag(&w[jmax + (kw - 1) * w_dim1]), abs(d__2));
                if (imax > 1) {
                    i__1 = imax - 1;
                    jmax = izamax_(&i__1, &w[(kw - 1) * w_dim1 + 1], &c__1);
                    i__1 = jmax + (kw - 1) * w_dim1;
                    d__3 = rowmax, d__4 = (d__1 = w[i__1].r, abs(d__1)) +
                                          (d__2 = d_lmp_imag(&w[jmax + (kw - 1) * w_dim1]), abs(d__2));
                    rowmax = max(d__3, d__4);
                }
                if (absakk >= alpha * colmax * (colmax / rowmax)) {
                    kp = k;
                } else {
                    i__1 = imax + (kw - 1) * w_dim1;
                    if ((d__1 = w[i__1].r, abs(d__1)) +
                            (d__2 = d_lmp_imag(&w[imax + (kw - 1) * w_dim1]), abs(d__2)) >=
                        alpha * rowmax) {
                        kp = imax;
                        zcopy_(&k, &w[(kw - 1) * w_dim1 + 1], &c__1, &w[kw * w_dim1 + 1], &c__1);
                    } else {
                        kp = imax;
                        kstep = 2;
                    }
                }
            }
            kk = k - kstep + 1;
            kkw = *nb + kk - *n;
            if (kp != kk) {
                i__1 = kp + kp * a_dim1;
                i__2 = kk + kk * a_dim1;
                a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
                i__1 = kk - 1 - kp;
                zcopy_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (kp + 1) * a_dim1], lda);
                if (kp > 1) {
                    i__1 = kp - 1;
                    zcopy_(&i__1, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &c__1);
                }
                if (k < *n) {
                    i__1 = *n - k;
                    zswap_(&i__1, &a[kk + (k + 1) * a_dim1], lda, &a[kp + (k + 1) * a_dim1], lda);
                }
                i__1 = *n - kk + 1;
                zswap_(&i__1, &w[kk + kkw * w_dim1], ldw, &w[kp + kkw * w_dim1], ldw);
            }
            if (kstep == 1) {
                zcopy_(&k, &w[kw * w_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
                z_lmp_div(&z__1, &c_b1, &a[k + k * a_dim1]);
                r1.r = z__1.r, r1.i = z__1.i;
                i__1 = k - 1;
                zscal_(&i__1, &r1, &a[k * a_dim1 + 1], &c__1);
            } else {
                if (k > 2) {
                    i__1 = k - 1 + kw * w_dim1;
                    d21.r = w[i__1].r, d21.i = w[i__1].i;
                    z_lmp_div(&z__1, &w[k + kw * w_dim1], &d21);
                    d11.r = z__1.r, d11.i = z__1.i;
                    z_lmp_div(&z__1, &w[k - 1 + (kw - 1) * w_dim1], &d21);
                    d22.r = z__1.r, d22.i = z__1.i;
                    z__3.r = d11.r * d22.r - d11.i * d22.i, z__3.i = d11.r * d22.i + d11.i * d22.r;
                    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
                    z_lmp_div(&z__1, &c_b1, &z__2);
                    t.r = z__1.r, t.i = z__1.i;
                    z_lmp_div(&z__1, &t, &d21);
                    d21.r = z__1.r, d21.i = z__1.i;
                    i__1 = k - 2;
                    for (j = 1; j <= i__1; ++j) {
                        i__2 = j + (k - 1) * a_dim1;
                        i__3 = j + (kw - 1) * w_dim1;
                        z__3.r = d11.r * w[i__3].r - d11.i * w[i__3].i,
                        z__3.i = d11.r * w[i__3].i + d11.i * w[i__3].r;
                        i__4 = j + kw * w_dim1;
                        z__2.r = z__3.r - w[i__4].r, z__2.i = z__3.i - w[i__4].i;
                        z__1.r = d21.r * z__2.r - d21.i * z__2.i,
                        z__1.i = d21.r * z__2.i + d21.i * z__2.r;
                        a[i__2].r = z__1.r, a[i__2].i = z__1.i;
                        i__2 = j + k * a_dim1;
                        i__3 = j + kw * w_dim1;
                        z__3.r = d22.r * w[i__3].r - d22.i * w[i__3].i,
                        z__3.i = d22.r * w[i__3].i + d22.i * w[i__3].r;
                        i__4 = j + (kw - 1) * w_dim1;
                        z__2.r = z__3.r - w[i__4].r, z__2.i = z__3.i - w[i__4].i;
                        z__1.r = d21.r * z__2.r - d21.i * z__2.i,
                        z__1.i = d21.r * z__2.i + d21.i * z__2.r;
                        a[i__2].r = z__1.r, a[i__2].i = z__1.i;
                    }
                }
                i__1 = k - 1 + (k - 1) * a_dim1;
                i__2 = k - 1 + (kw - 1) * w_dim1;
                a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
                i__1 = k - 1 + k * a_dim1;
                i__2 = k - 1 + kw * w_dim1;
                a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
                i__1 = k + k * a_dim1;
                i__2 = k + kw * w_dim1;
                a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
            }
        }
        if (kstep == 1) {
            ipiv[k] = kp;
        } else {
            ipiv[k] = -kp;
            ipiv[k - 1] = -kp;
        }
        k -= kstep;
        goto L10;
    L30:
        i__1 = *n - k;
        z__1.r = -1., z__1.i = -0.;
        zgemmtr_((char *)"U", (char *)"N", (char *)"T", &k, &i__1, &z__1, &a[(k + 1) * a_dim1 + 1], lda,
                 &w[(kw + 1) * w_dim1 + 1], ldw, &c_b1, &a[a_dim1 + 1], lda, (ftnlen)1, (ftnlen)1,
                 (ftnlen)1);
        j = k + 1;
    L60:
        jj = j;
        jp = ipiv[j];
        if (jp < 0) {
            jp = -jp;
            ++j;
        }
        ++j;
        if (jp != jj && j <= *n) {
            i__1 = *n - j + 1;
            zswap_(&i__1, &a[jp + j * a_dim1], lda, &a[jj + j * a_dim1], lda);
        }
        if (j < *n) {
            goto L60;
        }
        *kb = *n - k;
    } else {
        k = 1;
    L70:
        if (k >= *nb && *nb < *n || k > *n) {
            goto L90;
        }
        i__1 = *n - k + 1;
        zcopy_(&i__1, &a[k + k * a_dim1], &c__1, &w[k + k * w_dim1], &c__1);
        i__1 = *n - k + 1;
        i__2 = k - 1;
        z__1.r = -1., z__1.i = -0.;
        zgemv_((char *)"N", &i__1, &i__2, &z__1, &a[k + a_dim1], lda, &w[k + w_dim1], ldw, &c_b1,
               &w[k + k * w_dim1], &c__1, (ftnlen)1);
        kstep = 1;
        i__1 = k + k * w_dim1;
        absakk = (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_lmp_imag(&w[k + k * w_dim1]), abs(d__2));
        if (k < *n) {
            i__1 = *n - k;
            imax = k + izamax_(&i__1, &w[k + 1 + k * w_dim1], &c__1);
            i__1 = imax + k * w_dim1;
            colmax =
                (d__1 = w[i__1].r, abs(d__1)) + (d__2 = d_lmp_imag(&w[imax + k * w_dim1]), abs(d__2));
        } else {
            colmax = 0.;
        }
        if (max(absakk, colmax) == 0.) {
            if (*info == 0) {
                *info = k;
            }
            kp = k;
        } else {
            if (absakk >= alpha * colmax) {
                kp = k;
            } else {
                i__1 = imax - k;
                zcopy_(&i__1, &a[imax + k * a_dim1], lda, &w[k + (k + 1) * w_dim1], &c__1);
                i__1 = *n - imax + 1;
                zcopy_(&i__1, &a[imax + imax * a_dim1], &c__1, &w[imax + (k + 1) * w_dim1], &c__1);
                i__1 = *n - k + 1;
                i__2 = k - 1;
                z__1.r = -1., z__1.i = -0.;
                zgemv_((char *)"N", &i__1, &i__2, &z__1, &a[k + a_dim1], lda, &w[imax + w_dim1], ldw, &c_b1,
                       &w[k + (k + 1) * w_dim1], &c__1, (ftnlen)1);
                i__1 = imax - k;
                jmax = k - 1 + izamax_(&i__1, &w[k + (k + 1) * w_dim1], &c__1);
                i__1 = jmax + (k + 1) * w_dim1;
                rowmax = (d__1 = w[i__1].r, abs(d__1)) +
                         (d__2 = d_lmp_imag(&w[jmax + (k + 1) * w_dim1]), abs(d__2));
                if (imax < *n) {
                    i__1 = *n - imax;
                    jmax = imax + izamax_(&i__1, &w[imax + 1 + (k + 1) * w_dim1], &c__1);
                    i__1 = jmax + (k + 1) * w_dim1;
                    d__3 = rowmax, d__4 = (d__1 = w[i__1].r, abs(d__1)) +
                                          (d__2 = d_lmp_imag(&w[jmax + (k + 1) * w_dim1]), abs(d__2));
                    rowmax = max(d__3, d__4);
                }
                if (absakk >= alpha * colmax * (colmax / rowmax)) {
                    kp = k;
                } else {
                    i__1 = imax + (k + 1) * w_dim1;
                    if ((d__1 = w[i__1].r, abs(d__1)) +
                            (d__2 = d_lmp_imag(&w[imax + (k + 1) * w_dim1]), abs(d__2)) >=
                        alpha * rowmax) {
                        kp = imax;
                        i__1 = *n - k + 1;
                        zcopy_(&i__1, &w[k + (k + 1) * w_dim1], &c__1, &w[k + k * w_dim1], &c__1);
                    } else {
                        kp = imax;
                        kstep = 2;
                    }
                }
            }
            kk = k + kstep - 1;
            if (kp != kk) {
                i__1 = kp + kp * a_dim1;
                i__2 = kk + kk * a_dim1;
                a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
                i__1 = kp - kk - 1;
                zcopy_(&i__1, &a[kk + 1 + kk * a_dim1], &c__1, &a[kp + (kk + 1) * a_dim1], lda);
                if (kp < *n) {
                    i__1 = *n - kp;
                    zcopy_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + 1 + kp * a_dim1], &c__1);
                }
                if (k > 1) {
                    i__1 = k - 1;
                    zswap_(&i__1, &a[kk + a_dim1], lda, &a[kp + a_dim1], lda);
                }
                zswap_(&kk, &w[kk + w_dim1], ldw, &w[kp + w_dim1], ldw);
            }
            if (kstep == 1) {
                i__1 = *n - k + 1;
                zcopy_(&i__1, &w[k + k * w_dim1], &c__1, &a[k + k * a_dim1], &c__1);
                if (k < *n) {
                    z_lmp_div(&z__1, &c_b1, &a[k + k * a_dim1]);
                    r1.r = z__1.r, r1.i = z__1.i;
                    i__1 = *n - k;
                    zscal_(&i__1, &r1, &a[k + 1 + k * a_dim1], &c__1);
                }
            } else {
                if (k < *n - 1) {
                    i__1 = k + 1 + k * w_dim1;
                    d21.r = w[i__1].r, d21.i = w[i__1].i;
                    z_lmp_div(&z__1, &w[k + 1 + (k + 1) * w_dim1], &d21);
                    d11.r = z__1.r, d11.i = z__1.i;
                    z_lmp_div(&z__1, &w[k + k * w_dim1], &d21);
                    d22.r = z__1.r, d22.i = z__1.i;
                    z__3.r = d11.r * d22.r - d11.i * d22.i, z__3.i = d11.r * d22.i + d11.i * d22.r;
                    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
                    z_lmp_div(&z__1, &c_b1, &z__2);
                    t.r = z__1.r, t.i = z__1.i;
                    z_lmp_div(&z__1, &t, &d21);
                    d21.r = z__1.r, d21.i = z__1.i;
                    i__1 = *n;
                    for (j = k + 2; j <= i__1; ++j) {
                        i__2 = j + k * a_dim1;
                        i__3 = j + k * w_dim1;
                        z__3.r = d11.r * w[i__3].r - d11.i * w[i__3].i,
                        z__3.i = d11.r * w[i__3].i + d11.i * w[i__3].r;
                        i__4 = j + (k + 1) * w_dim1;
                        z__2.r = z__3.r - w[i__4].r, z__2.i = z__3.i - w[i__4].i;
                        z__1.r = d21.r * z__2.r - d21.i * z__2.i,
                        z__1.i = d21.r * z__2.i + d21.i * z__2.r;
                        a[i__2].r = z__1.r, a[i__2].i = z__1.i;
                        i__2 = j + (k + 1) * a_dim1;
                        i__3 = j + (k + 1) * w_dim1;
                        z__3.r = d22.r * w[i__3].r - d22.i * w[i__3].i,
                        z__3.i = d22.r * w[i__3].i + d22.i * w[i__3].r;
                        i__4 = j + k * w_dim1;
                        z__2.r = z__3.r - w[i__4].r, z__2.i = z__3.i - w[i__4].i;
                        z__1.r = d21.r * z__2.r - d21.i * z__2.i,
                        z__1.i = d21.r * z__2.i + d21.i * z__2.r;
                        a[i__2].r = z__1.r, a[i__2].i = z__1.i;
                    }
                }
                i__1 = k + k * a_dim1;
                i__2 = k + k * w_dim1;
                a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
                i__1 = k + 1 + k * a_dim1;
                i__2 = k + 1 + k * w_dim1;
                a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
                i__1 = k + 1 + (k + 1) * a_dim1;
                i__2 = k + 1 + (k + 1) * w_dim1;
                a[i__1].r = w[i__2].r, a[i__1].i = w[i__2].i;
            }
        }
        if (kstep == 1) {
            ipiv[k] = kp;
        } else {
            ipiv[k] = -kp;
            ipiv[k + 1] = -kp;
        }
        k += kstep;
        goto L70;
    L90:
        i__1 = *n - k + 1;
        i__2 = k - 1;
        z__1.r = -1., z__1.i = -0.;
        zgemmtr_((char *)"L", (char *)"N", (char *)"T", &i__1, &i__2, &z__1, &a[k + a_dim1], lda, &w[k + w_dim1], ldw,
                 &c_b1, &a[k + k * a_dim1], lda, (ftnlen)1, (ftnlen)1, (ftnlen)1);
        j = k - 1;
    L120:
        jj = j;
        jp = ipiv[j];
        if (jp < 0) {
            jp = -jp;
            --j;
        }
        --j;
        if (jp != jj && j >= 1) {
            zswap_(&j, &a[jp + a_dim1], lda, &a[jj + a_dim1], lda);
        }
        if (j > 1) {
            goto L120;
        }
        *kb = k - 1;
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
