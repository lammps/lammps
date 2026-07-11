#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b11 = 0.;
static doublereal c_b12 = 1.;
static integer c__12 = 12;
static integer c__2 = 2;
static integer c__49 = 49;
int dhseqr_(char *job, char *compz, integer *n, integer *ilo, integer *ihi, doublereal *h__,
            integer *ldh, doublereal *wr, doublereal *wi, doublereal *z__, integer *ldz,
            doublereal *work, integer *lwork, integer *info, ftnlen job_len, ftnlen compz_len)
{
    address a__1[2];
    integer h_dim1, h_offset, z_dim1, z_offset, i__1, i__2[2], i__3;
    doublereal d__1;
    char ch__1[2];
    int s_lmp_cat(char *, char **, integer *, integer *, ftnlen);
    integer i__;
    doublereal hl[2401];
    integer kbot, nmin;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    logical initz;
    doublereal workl[49];
    logical wantt, wantz;
    extern int dlaqr0_(logical *, logical *, integer *, integer *, integer *, doublereal *,
                       integer *, doublereal *, doublereal *, integer *, integer *, doublereal *,
                       integer *, doublereal *, integer *, integer *),
        dlahqr_(logical *, logical *, integer *, integer *, integer *, doublereal *, integer *,
                doublereal *, doublereal *, integer *, integer *, doublereal *, integer *,
                integer *),
        dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *,
                ftnlen),
        dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    extern int xerbla_(char *, integer *, ftnlen);
    logical lquery;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --wr;
    --wi;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    wantt = lsame_(job, (char *)"S", (ftnlen)1, (ftnlen)1);
    initz = lsame_(compz, (char *)"I", (ftnlen)1, (ftnlen)1);
    wantz = initz || lsame_(compz, (char *)"V", (ftnlen)1, (ftnlen)1);
    work[1] = (doublereal)max(1, *n);
    lquery = *lwork == -1;
    *info = 0;
    if (!lsame_(job, (char *)"E", (ftnlen)1, (ftnlen)1) && !wantt) {
        *info = -1;
    } else if (!lsame_(compz, (char *)"N", (ftnlen)1, (ftnlen)1) && !wantz) {
        *info = -2;
    } else if (*n < 0) {
        *info = -3;
    } else if (*ilo < 1 || *ilo > max(1, *n)) {
        *info = -4;
    } else if (*ihi < min(*ilo, *n) || *ihi > *n) {
        *info = -5;
    } else if (*ldh < max(1, *n)) {
        *info = -7;
    } else if (*ldz < 1 || wantz && *ldz < max(1, *n)) {
        *info = -11;
    } else if (*lwork < max(1, *n) && !lquery) {
        *info = -13;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DHSEQR", &i__1, (ftnlen)6);
        return 0;
    } else if (*n == 0) {
        return 0;
    } else if (lquery) {
        dlaqr0_(&wantt, &wantz, n, ilo, ihi, &h__[h_offset], ldh, &wr[1], &wi[1], ilo, ihi,
                &z__[z_offset], ldz, &work[1], lwork, info);
        d__1 = (doublereal)max(1, *n);
        work[1] = max(d__1, work[1]);
        return 0;
    } else {
        i__1 = *ilo - 1;
        for (i__ = 1; i__ <= i__1; ++i__) {
            wr[i__] = h__[i__ + i__ * h_dim1];
            wi[i__] = 0.;
        }
        i__1 = *n;
        for (i__ = *ihi + 1; i__ <= i__1; ++i__) {
            wr[i__] = h__[i__ + i__ * h_dim1];
            wi[i__] = 0.;
        }
        if (initz) {
            dlaset_((char *)"A", n, n, &c_b11, &c_b12, &z__[z_offset], ldz, (ftnlen)1);
        }
        if (*ilo == *ihi) {
            wr[*ilo] = h__[*ilo + *ilo * h_dim1];
            wi[*ilo] = 0.;
            return 0;
        }
        i__2[0] = 1, a__1[0] = job;
        i__2[1] = 1, a__1[1] = compz;
        s_lmp_cat(ch__1, a__1, i__2, &c__2, (ftnlen)2);
        nmin = ilaenv_(&c__12, (char *)"DHSEQR", ch__1, n, ilo, ihi, lwork, (ftnlen)6, (ftnlen)2);
        nmin = max(15, nmin);
        if (*n > nmin) {
            dlaqr0_(&wantt, &wantz, n, ilo, ihi, &h__[h_offset], ldh, &wr[1], &wi[1], ilo, ihi,
                    &z__[z_offset], ldz, &work[1], lwork, info);
        } else {
            dlahqr_(&wantt, &wantz, n, ilo, ihi, &h__[h_offset], ldh, &wr[1], &wi[1], ilo, ihi,
                    &z__[z_offset], ldz, info);
            if (*info > 0) {
                kbot = *info;
                if (*n >= 49) {
                    dlaqr0_(&wantt, &wantz, n, ilo, &kbot, &h__[h_offset], ldh, &wr[1], &wi[1], ilo,
                            ihi, &z__[z_offset], ldz, &work[1], lwork, info);
                } else {
                    dlacpy_((char *)"A", n, n, &h__[h_offset], ldh, hl, &c__49, (ftnlen)1);
                    hl[*n + 1 + *n * 49 - 50] = 0.;
                    i__1 = 49 - *n;
                    dlaset_((char *)"A", &c__49, &i__1, &c_b11, &c_b11, &hl[(*n + 1) * 49 - 49], &c__49,
                            (ftnlen)1);
                    dlaqr0_(&wantt, &wantz, &c__49, ilo, &kbot, hl, &c__49, &wr[1], &wi[1], ilo,
                            ihi, &z__[z_offset], ldz, workl, &c__49, info);
                    if (wantt || *info != 0) {
                        dlacpy_((char *)"A", n, n, hl, &c__49, &h__[h_offset], ldh, (ftnlen)1);
                    }
                }
            }
        }
        if ((wantt || *info != 0) && *n > 2) {
            i__1 = *n - 2;
            i__3 = *n - 2;
            dlaset_((char *)"L", &i__1, &i__3, &c_b11, &c_b11, &h__[h_dim1 + 3], ldh, (ftnlen)1);
        }
        d__1 = (doublereal)max(1, *n);
        work[1] = max(d__1, work[1]);
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
