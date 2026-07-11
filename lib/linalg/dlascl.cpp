#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlascl_(char *type__, integer *kl, integer *ku, doublereal *cfrom, doublereal *cto, integer *m,
            integer *n, doublereal *a, integer *lda, integer *info, ftnlen type_len)
{
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    integer i__, j, k1, k2, k3, k4;
    doublereal mul, cto1;
    logical done;
    doublereal ctoc;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer itype;
    doublereal cfrom1;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal cfromc;
    extern logical disnan_(doublereal *);
    extern int xerbla_(char *, integer *, ftnlen);
    doublereal bignum, smlnum;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    *info = 0;
    if (lsame_(type__, (char *)"G", (ftnlen)1, (ftnlen)1)) {
        itype = 0;
    } else if (lsame_(type__, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        itype = 1;
    } else if (lsame_(type__, (char *)"U", (ftnlen)1, (ftnlen)1)) {
        itype = 2;
    } else if (lsame_(type__, (char *)"H", (ftnlen)1, (ftnlen)1)) {
        itype = 3;
    } else if (lsame_(type__, (char *)"B", (ftnlen)1, (ftnlen)1)) {
        itype = 4;
    } else if (lsame_(type__, (char *)"Q", (ftnlen)1, (ftnlen)1)) {
        itype = 5;
    } else if (lsame_(type__, (char *)"Z", (ftnlen)1, (ftnlen)1)) {
        itype = 6;
    } else {
        itype = -1;
    }
    if (itype == -1) {
        *info = -1;
    } else if (*cfrom == 0. || disnan_(cfrom)) {
        *info = -4;
    } else if (disnan_(cto)) {
        *info = -5;
    } else if (*m < 0) {
        *info = -6;
    } else if (*n < 0 || itype == 4 && *n != *m || itype == 5 && *n != *m) {
        *info = -7;
    } else if (itype <= 3 && *lda < max(1, *m)) {
        *info = -9;
    } else if (itype >= 4) {
        i__1 = *m - 1;
        if (*kl < 0 || *kl > max(i__1, 0)) {
            *info = -2;
        } else {
            i__1 = *n - 1;
            if (*ku < 0 || *ku > max(i__1, 0) || (itype == 4 || itype == 5) && *kl != *ku) {
                *info = -3;
            } else if (itype == 4 && *lda < *kl + 1 || itype == 5 && *lda < *ku + 1 ||
                       itype == 6 && *lda < (*kl << 1) + *ku + 1) {
                *info = -9;
            }
        }
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLASCL", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0 || *m == 0) {
        return 0;
    }
    smlnum = dlamch_((char *)"S", (ftnlen)1);
    bignum = 1. / smlnum;
    cfromc = *cfrom;
    ctoc = *cto;
L10:
    cfrom1 = cfromc * smlnum;
    if (cfrom1 == cfromc) {
        mul = ctoc / cfromc;
        done = TRUE_;
        cto1 = ctoc;
    } else {
        cto1 = ctoc / bignum;
        if (cto1 == ctoc) {
            mul = ctoc;
            done = TRUE_;
            cfromc = 1.;
        } else if (abs(cfrom1) > abs(ctoc) && ctoc != 0.) {
            mul = smlnum;
            done = FALSE_;
            cfromc = cfrom1;
        } else if (abs(cto1) > abs(cfromc)) {
            mul = bignum;
            done = FALSE_;
            ctoc = cto1;
        } else {
            mul = ctoc / cfromc;
            done = TRUE_;
            if (mul == 1.) {
                return 0;
            }
        }
    }
    if (itype == 0) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__2 = *m;
            for (i__ = 1; i__ <= i__2; ++i__) {
                a[i__ + j * a_dim1] *= mul;
            }
        }
    } else if (itype == 1) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__2 = *m;
            for (i__ = j; i__ <= i__2; ++i__) {
                a[i__ + j * a_dim1] *= mul;
            }
        }
    } else if (itype == 2) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__2 = min(j, *m);
            for (i__ = 1; i__ <= i__2; ++i__) {
                a[i__ + j * a_dim1] *= mul;
            }
        }
    } else if (itype == 3) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__3 = j + 1;
            i__2 = min(i__3, *m);
            for (i__ = 1; i__ <= i__2; ++i__) {
                a[i__ + j * a_dim1] *= mul;
            }
        }
    } else if (itype == 4) {
        k3 = *kl + 1;
        k4 = *n + 1;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__3 = k3, i__4 = k4 - j;
            i__2 = min(i__3, i__4);
            for (i__ = 1; i__ <= i__2; ++i__) {
                a[i__ + j * a_dim1] *= mul;
            }
        }
    } else if (itype == 5) {
        k1 = *ku + 2;
        k3 = *ku + 1;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__2 = k1 - j;
            i__3 = k3;
            for (i__ = max(i__2, 1); i__ <= i__3; ++i__) {
                a[i__ + j * a_dim1] *= mul;
            }
        }
    } else if (itype == 6) {
        k1 = *kl + *ku + 2;
        k2 = *kl + 1;
        k3 = (*kl << 1) + *ku + 1;
        k4 = *kl + *ku + 1 + *m;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__3 = k1 - j;
            i__4 = k3, i__5 = k4 - j;
            i__2 = min(i__4, i__5);
            for (i__ = max(i__3, k2); i__ <= i__2; ++i__) {
                a[i__ + j * a_dim1] *= mul;
            }
        }
    }
    if (!done) {
        goto L10;
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
