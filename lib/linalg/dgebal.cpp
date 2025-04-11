#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
int dgebal_(char *job, integer *n, doublereal *a, integer *lda, integer *ilo, integer *ihi,
            doublereal *scale, integer *info, ftnlen job_len)
{
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;
    doublereal c__, f, g;
    integer i__, j, k, l;
    doublereal r__, s, ca, ra;
    integer ica, ira;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern int dscal_(integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
    doublereal sfmin1, sfmin2, sfmax1, sfmax2;
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern logical disnan_(doublereal *);
    extern int xerbla_(char *, integer *, ftnlen);
    logical noconv, canswap;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --scale;
    *info = 0;
    if (!lsame_(job, (char *)"N", (ftnlen)1, (ftnlen)1) && !lsame_(job, (char *)"P", (ftnlen)1, (ftnlen)1) &&
        !lsame_(job, (char *)"S", (ftnlen)1, (ftnlen)1) && !lsame_(job, (char *)"B", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*lda < max(1, *n)) {
        *info = -4;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DGEBAL", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        *ilo = 1;
        *ihi = 0;
        return 0;
    }
    if (lsame_(job, (char *)"N", (ftnlen)1, (ftnlen)1)) {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            scale[i__] = 1.;
        }
        *ilo = 1;
        *ihi = *n;
        return 0;
    }
    k = 1;
    l = *n;
    if (!lsame_(job, (char *)"S", (ftnlen)1, (ftnlen)1)) {
        noconv = TRUE_;
        while (noconv) {
            noconv = FALSE_;
            for (i__ = l; i__ >= 1; --i__) {
                canswap = TRUE_;
                i__1 = l;
                for (j = 1; j <= i__1; ++j) {
                    if (i__ != j && a[i__ + j * a_dim1] != 0.) {
                        canswap = FALSE_;
                        goto L100;
                    }
                }
            L100:
                if (canswap) {
                    scale[l] = (doublereal)i__;
                    if (i__ != l) {
                        dswap_(&l, &a[i__ * a_dim1 + 1], &c__1, &a[l * a_dim1 + 1], &c__1);
                        i__1 = *n - k + 1;
                        dswap_(&i__1, &a[i__ + k * a_dim1], lda, &a[l + k * a_dim1], lda);
                    }
                    noconv = TRUE_;
                    if (l == 1) {
                        *ilo = 1;
                        *ihi = 1;
                        return 0;
                    }
                    --l;
                }
            }
        }
        noconv = TRUE_;
        while (noconv) {
            noconv = FALSE_;
            i__1 = l;
            for (j = k; j <= i__1; ++j) {
                canswap = TRUE_;
                i__2 = l;
                for (i__ = k; i__ <= i__2; ++i__) {
                    if (i__ != j && a[i__ + j * a_dim1] != 0.) {
                        canswap = FALSE_;
                        goto L200;
                    }
                }
            L200:
                if (canswap) {
                    scale[k] = (doublereal)j;
                    if (j != k) {
                        dswap_(&l, &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
                        i__2 = *n - k + 1;
                        dswap_(&i__2, &a[j + k * a_dim1], lda, &a[k + k * a_dim1], lda);
                    }
                    noconv = TRUE_;
                    ++k;
                }
            }
        }
    }
    i__1 = l;
    for (i__ = k; i__ <= i__1; ++i__) {
        scale[i__] = 1.;
    }
    if (lsame_(job, (char *)"P", (ftnlen)1, (ftnlen)1)) {
        *ilo = k;
        *ihi = l;
        return 0;
    }
    sfmin1 = dlamch_((char *)"S", (ftnlen)1) / dlamch_((char *)"P", (ftnlen)1);
    sfmax1 = 1. / sfmin1;
    sfmin2 = sfmin1 * 2.;
    sfmax2 = 1. / sfmin2;
    noconv = TRUE_;
    while (noconv) {
        noconv = FALSE_;
        i__1 = l;
        for (i__ = k; i__ <= i__1; ++i__) {
            i__2 = l - k + 1;
            c__ = dnrm2_(&i__2, &a[k + i__ * a_dim1], &c__1);
            i__2 = l - k + 1;
            r__ = dnrm2_(&i__2, &a[i__ + k * a_dim1], lda);
            ica = idamax_(&l, &a[i__ * a_dim1 + 1], &c__1);
            ca = (d__1 = a[ica + i__ * a_dim1], abs(d__1));
            i__2 = *n - k + 1;
            ira = idamax_(&i__2, &a[i__ + k * a_dim1], lda);
            ra = (d__1 = a[i__ + (ira + k - 1) * a_dim1], abs(d__1));
            if (c__ == 0. || r__ == 0.) {
                goto L300;
            }
            d__1 = c__ + ca + r__ + ra;
            if (disnan_(&d__1)) {
                *info = -3;
                i__2 = -(*info);
                xerbla_((char *)"DGEBAL", &i__2, (ftnlen)6);
                return 0;
            }
            g = r__ / 2.;
            f = 1.;
            s = c__ + r__;
            for (;;) {
                d__1 = max(f, c__);
                d__2 = min(r__, g);
                if (!(c__ < g && max(d__1, ca) < sfmax2 && min(d__2, ra) > sfmin2)) break;
                f *= 2.;
                c__ *= 2.;
                ca *= 2.;
                r__ /= 2.;
                g /= 2.;
                ra /= 2.;
            }
            g = c__ / 2.;
            for (;;) {
                d__1 = min(f, c__), d__1 = min(d__1, g);
                if (!(g >= r__ && max(r__, ra) < sfmax2 && min(d__1, ca) > sfmin2)) break;
                f /= 2.;
                c__ /= 2.;
                g /= 2.;
                ca /= 2.;
                r__ *= 2.;
                ra *= 2.;
            }
            if (c__ + r__ >= s * .95) {
                goto L300;
            }
            if (f < 1. && scale[i__] < 1.) {
                if (f * scale[i__] <= sfmin1) {
                    goto L300;
                }
            }
            if (f > 1. && scale[i__] > 1.) {
                if (scale[i__] >= sfmax1 / f) {
                    goto L300;
                }
            }
            g = 1. / f;
            scale[i__] *= f;
            noconv = TRUE_;
            i__2 = *n - k + 1;
            dscal_(&i__2, &g, &a[i__ + k * a_dim1], lda);
            dscal_(&l, &f, &a[i__ * a_dim1 + 1], &c__1);
        L300:;
        }
    }
    *ilo = k;
    *ihi = l;
    return 0;
}
#ifdef __cplusplus
}
#endif
