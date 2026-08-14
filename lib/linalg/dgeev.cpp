#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;
int dgeev_(char *jobvl, char *jobvr, integer *n, doublereal *a, integer *lda, doublereal *wr,
           doublereal *wi, doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr,
           doublereal *work, integer *lwork, integer *info, ftnlen jobvl_len, ftnlen jobvr_len)
{
    integer a_dim1, a_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;
    double sqrt(doublereal);
    integer i__, k;
    doublereal r__, cs, sn;
    integer ihi;
    doublereal scl;
    integer ilo;
    doublereal dum[1], eps;
    integer lwork_trevc__, ibal;
    char side[1];
    doublereal anrm;
    integer ierr, itau;
    extern int drot_(integer *, doublereal *, integer *, doublereal *, integer *, doublereal *,
                     doublereal *);
    integer iwrk, nout;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern int dscal_(integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern doublereal dlapy2_(doublereal *, doublereal *);
    extern int dgebak_(char *, char *, integer *, integer *, integer *, doublereal *, integer *,
                       doublereal *, integer *, integer *, ftnlen, ftnlen),
        dgebal_(char *, integer *, doublereal *, integer *, integer *, integer *, doublereal *,
                integer *, ftnlen);
    logical scalea;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal cscale;
    extern doublereal dlange_(char *, integer *, integer *, doublereal *, integer *, doublereal *,
                              ftnlen);
    extern int dgehrd_(integer *, integer *, integer *, doublereal *, integer *, doublereal *,
                       doublereal *, integer *, integer *),
        dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *,
                doublereal *, integer *, integer *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern int dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *,
                       integer *, ftnlen),
        dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *),
        xerbla_(char *, integer *, ftnlen);
    logical select[1];
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    doublereal bignum;
    extern int dorghr_(integer *, integer *, integer *, doublereal *, integer *, doublereal *,
                       doublereal *, integer *, integer *),
        dhseqr_(char *, char *, integer *, integer *, integer *, doublereal *, integer *,
                doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *,
                integer *, ftnlen, ftnlen);
    integer minwrk, maxwrk;
    logical wantvl;
    doublereal smlnum;
    integer hswork;
    logical lquery, wantvr;
    extern int dtrevc3_(char *, char *, logical *, integer *, doublereal *, integer *, doublereal *,
                        integer *, doublereal *, integer *, integer *, integer *, doublereal *,
                        integer *, integer *, ftnlen, ftnlen);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --wr;
    --wi;
    vl_dim1 = *ldvl;
    vl_offset = 1 + vl_dim1;
    vl -= vl_offset;
    vr_dim1 = *ldvr;
    vr_offset = 1 + vr_dim1;
    vr -= vr_offset;
    --work;
    *info = 0;
    lquery = *lwork == -1;
    wantvl = lsame_(jobvl, (char *)"V", (ftnlen)1, (ftnlen)1);
    wantvr = lsame_(jobvr, (char *)"V", (ftnlen)1, (ftnlen)1);
    if (!wantvl && !lsame_(jobvl, (char *)"N", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (!wantvr && !lsame_(jobvr, (char *)"N", (ftnlen)1, (ftnlen)1)) {
        *info = -2;
    } else if (*n < 0) {
        *info = -3;
    } else if (*lda < max(1, *n)) {
        *info = -5;
    } else if (*ldvl < 1 || wantvl && *ldvl < *n) {
        *info = -9;
    } else if (*ldvr < 1 || wantvr && *ldvr < *n) {
        *info = -11;
    }
    if (*info == 0) {
        if (*n == 0) {
            minwrk = 1;
            maxwrk = 1;
        } else {
            maxwrk = (*n << 1) +
                     *n * ilaenv_(&c__1, (char *)"DGEHRD", (char *)" ", n, &c__1, n, &c__0, (ftnlen)6, (ftnlen)1);
            if (wantvl) {
                minwrk = *n << 2;
                i__1 = maxwrk,
                i__2 = (*n << 1) + (*n - 1) * ilaenv_(&c__1, (char *)"DORGHR", (char *)" ", n, &c__1, n, &c_n1,
                                                      (ftnlen)6, (ftnlen)1);
                maxwrk = max(i__1, i__2);
                dhseqr_((char *)"S", (char *)"V", n, &c__1, n, &a[a_offset], lda, &wr[1], &wi[1], &vl[vl_offset],
                        ldvl, &work[1], &c_n1, info, (ftnlen)1, (ftnlen)1);
                hswork = (integer)work[1];
                i__1 = maxwrk, i__2 = *n + 1, i__1 = max(i__1, i__2), i__2 = *n + hswork;
                maxwrk = max(i__1, i__2);
                dtrevc3_((char *)"L", (char *)"B", select, n, &a[a_offset], lda, &vl[vl_offset], ldvl,
                         &vr[vr_offset], ldvr, n, &nout, &work[1], &c_n1, &ierr, (ftnlen)1,
                         (ftnlen)1);
                lwork_trevc__ = (integer)work[1];
                i__1 = maxwrk, i__2 = *n + lwork_trevc__;
                maxwrk = max(i__1, i__2);
                i__1 = maxwrk, i__2 = *n << 2;
                maxwrk = max(i__1, i__2);
            } else if (wantvr) {
                minwrk = *n << 2;
                i__1 = maxwrk,
                i__2 = (*n << 1) + (*n - 1) * ilaenv_(&c__1, (char *)"DORGHR", (char *)" ", n, &c__1, n, &c_n1,
                                                      (ftnlen)6, (ftnlen)1);
                maxwrk = max(i__1, i__2);
                dhseqr_((char *)"S", (char *)"V", n, &c__1, n, &a[a_offset], lda, &wr[1], &wi[1], &vr[vr_offset],
                        ldvr, &work[1], &c_n1, info, (ftnlen)1, (ftnlen)1);
                hswork = (integer)work[1];
                i__1 = maxwrk, i__2 = *n + 1, i__1 = max(i__1, i__2), i__2 = *n + hswork;
                maxwrk = max(i__1, i__2);
                dtrevc3_((char *)"R", (char *)"B", select, n, &a[a_offset], lda, &vl[vl_offset], ldvl,
                         &vr[vr_offset], ldvr, n, &nout, &work[1], &c_n1, &ierr, (ftnlen)1,
                         (ftnlen)1);
                lwork_trevc__ = (integer)work[1];
                i__1 = maxwrk, i__2 = *n + lwork_trevc__;
                maxwrk = max(i__1, i__2);
                i__1 = maxwrk, i__2 = *n << 2;
                maxwrk = max(i__1, i__2);
            } else {
                minwrk = *n * 3;
                dhseqr_((char *)"E", (char *)"N", n, &c__1, n, &a[a_offset], lda, &wr[1], &wi[1], &vr[vr_offset],
                        ldvr, &work[1], &c_n1, info, (ftnlen)1, (ftnlen)1);
                hswork = (integer)work[1];
                i__1 = maxwrk, i__2 = *n + 1, i__1 = max(i__1, i__2), i__2 = *n + hswork;
                maxwrk = max(i__1, i__2);
            }
            maxwrk = max(maxwrk, minwrk);
        }
        work[1] = (doublereal)maxwrk;
        if (*lwork < minwrk && !lquery) {
            *info = -13;
        }
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DGEEV ", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    eps = dlamch_((char *)"P", (ftnlen)1);
    smlnum = dlamch_((char *)"S", (ftnlen)1);
    bignum = 1. / smlnum;
    smlnum = sqrt(smlnum) / eps;
    bignum = 1. / smlnum;
    anrm = dlange_((char *)"M", n, n, &a[a_offset], lda, dum, (ftnlen)1);
    scalea = FALSE_;
    if (anrm > 0. && anrm < smlnum) {
        scalea = TRUE_;
        cscale = smlnum;
    } else if (anrm > bignum) {
        scalea = TRUE_;
        cscale = bignum;
    }
    if (scalea) {
        dlascl_((char *)"G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, &ierr, (ftnlen)1);
    }
    ibal = 1;
    dgebal_((char *)"B", n, &a[a_offset], lda, &ilo, &ihi, &work[ibal], &ierr, (ftnlen)1);
    itau = ibal + *n;
    iwrk = itau + *n;
    i__1 = *lwork - iwrk + 1;
    dgehrd_(n, &ilo, &ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1, &ierr);
    if (wantvl) {
        *(unsigned char *)side = 'L';
        dlacpy_((char *)"L", n, n, &a[a_offset], lda, &vl[vl_offset], ldvl, (ftnlen)1);
        i__1 = *lwork - iwrk + 1;
        dorghr_(n, &ilo, &ihi, &vl[vl_offset], ldvl, &work[itau], &work[iwrk], &i__1, &ierr);
        iwrk = itau;
        i__1 = *lwork - iwrk + 1;
        dhseqr_((char *)"S", (char *)"V", n, &ilo, &ihi, &a[a_offset], lda, &wr[1], &wi[1], &vl[vl_offset], ldvl,
                &work[iwrk], &i__1, info, (ftnlen)1, (ftnlen)1);
        if (wantvr) {
            *(unsigned char *)side = 'B';
            dlacpy_((char *)"F", n, n, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, (ftnlen)1);
        }
    } else if (wantvr) {
        *(unsigned char *)side = 'R';
        dlacpy_((char *)"L", n, n, &a[a_offset], lda, &vr[vr_offset], ldvr, (ftnlen)1);
        i__1 = *lwork - iwrk + 1;
        dorghr_(n, &ilo, &ihi, &vr[vr_offset], ldvr, &work[itau], &work[iwrk], &i__1, &ierr);
        iwrk = itau;
        i__1 = *lwork - iwrk + 1;
        dhseqr_((char *)"S", (char *)"V", n, &ilo, &ihi, &a[a_offset], lda, &wr[1], &wi[1], &vr[vr_offset], ldvr,
                &work[iwrk], &i__1, info, (ftnlen)1, (ftnlen)1);
    } else {
        iwrk = itau;
        i__1 = *lwork - iwrk + 1;
        dhseqr_((char *)"E", (char *)"N", n, &ilo, &ihi, &a[a_offset], lda, &wr[1], &wi[1], &vr[vr_offset], ldvr,
                &work[iwrk], &i__1, info, (ftnlen)1, (ftnlen)1);
    }
    if (*info != 0) {
        goto L50;
    }
    if (wantvl || wantvr) {
        i__1 = *lwork - iwrk + 1;
        dtrevc3_(side, (char *)"B", select, n, &a[a_offset], lda, &vl[vl_offset], ldvl, &vr[vr_offset],
                 ldvr, n, &nout, &work[iwrk], &i__1, &ierr, (ftnlen)1, (ftnlen)1);
    }
    if (wantvl) {
        dgebak_((char *)"B", (char *)"L", n, &ilo, &ihi, &work[ibal], n, &vl[vl_offset], ldvl, &ierr, (ftnlen)1,
                (ftnlen)1);
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            if (wi[i__] == 0.) {
                scl = 1. / dnrm2_(n, &vl[i__ * vl_dim1 + 1], &c__1);
                dscal_(n, &scl, &vl[i__ * vl_dim1 + 1], &c__1);
            } else if (wi[i__] > 0.) {
                d__1 = dnrm2_(n, &vl[i__ * vl_dim1 + 1], &c__1);
                d__2 = dnrm2_(n, &vl[(i__ + 1) * vl_dim1 + 1], &c__1);
                scl = 1. / dlapy2_(&d__1, &d__2);
                dscal_(n, &scl, &vl[i__ * vl_dim1 + 1], &c__1);
                dscal_(n, &scl, &vl[(i__ + 1) * vl_dim1 + 1], &c__1);
                i__2 = *n;
                for (k = 1; k <= i__2; ++k) {
                    d__1 = vl[k + i__ * vl_dim1];
                    d__2 = vl[k + (i__ + 1) * vl_dim1];
                    work[iwrk + k - 1] = d__1 * d__1 + d__2 * d__2;
                }
                k = idamax_(n, &work[iwrk], &c__1);
                dlartg_(&vl[k + i__ * vl_dim1], &vl[k + (i__ + 1) * vl_dim1], &cs, &sn, &r__);
                drot_(n, &vl[i__ * vl_dim1 + 1], &c__1, &vl[(i__ + 1) * vl_dim1 + 1], &c__1, &cs,
                      &sn);
                vl[k + (i__ + 1) * vl_dim1] = 0.;
            }
        }
    }
    if (wantvr) {
        dgebak_((char *)"B", (char *)"R", n, &ilo, &ihi, &work[ibal], n, &vr[vr_offset], ldvr, &ierr, (ftnlen)1,
                (ftnlen)1);
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            if (wi[i__] == 0.) {
                scl = 1. / dnrm2_(n, &vr[i__ * vr_dim1 + 1], &c__1);
                dscal_(n, &scl, &vr[i__ * vr_dim1 + 1], &c__1);
            } else if (wi[i__] > 0.) {
                d__1 = dnrm2_(n, &vr[i__ * vr_dim1 + 1], &c__1);
                d__2 = dnrm2_(n, &vr[(i__ + 1) * vr_dim1 + 1], &c__1);
                scl = 1. / dlapy2_(&d__1, &d__2);
                dscal_(n, &scl, &vr[i__ * vr_dim1 + 1], &c__1);
                dscal_(n, &scl, &vr[(i__ + 1) * vr_dim1 + 1], &c__1);
                i__2 = *n;
                for (k = 1; k <= i__2; ++k) {
                    d__1 = vr[k + i__ * vr_dim1];
                    d__2 = vr[k + (i__ + 1) * vr_dim1];
                    work[iwrk + k - 1] = d__1 * d__1 + d__2 * d__2;
                }
                k = idamax_(n, &work[iwrk], &c__1);
                dlartg_(&vr[k + i__ * vr_dim1], &vr[k + (i__ + 1) * vr_dim1], &cs, &sn, &r__);
                drot_(n, &vr[i__ * vr_dim1 + 1], &c__1, &vr[(i__ + 1) * vr_dim1 + 1], &c__1, &cs,
                      &sn);
                vr[k + (i__ + 1) * vr_dim1] = 0.;
            }
        }
    }
L50:
    if (scalea) {
        i__1 = *n - *info;
        i__3 = *n - *info;
        i__2 = max(i__3, 1);
        dlascl_((char *)"G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wr[*info + 1], &i__2, &ierr,
                (ftnlen)1);
        i__1 = *n - *info;
        i__3 = *n - *info;
        i__2 = max(i__3, 1);
        dlascl_((char *)"G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[*info + 1], &i__2, &ierr,
                (ftnlen)1);
        if (*info > 0) {
            i__1 = ilo - 1;
            dlascl_((char *)"G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wr[1], n, &ierr, (ftnlen)1);
            i__1 = ilo - 1;
            dlascl_((char *)"G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[1], n, &ierr, (ftnlen)1);
        }
    }
    work[1] = (doublereal)maxwrk;
    return 0;
}
#ifdef __cplusplus
}
#endif
