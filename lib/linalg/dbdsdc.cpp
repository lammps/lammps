#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__9 = 9;
static integer c__0 = 0;
static doublereal c_b15 = 1.;
static integer c__1 = 1;
static doublereal c_b29 = 0.;
int dbdsdc_(char *uplo, char *compq, integer *n, doublereal *d__, doublereal *e, doublereal *u,
            integer *ldu, doublereal *vt, integer *ldvt, doublereal *q, integer *iq,
            doublereal *work, integer *iwork, integer *info, ftnlen uplo_len, ftnlen compq_len)
{
    integer u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2;
    doublereal d__1;
    double d_lmp_sign(doublereal *, doublereal *), log(doublereal);
    integer i__, j, k;
    doublereal p, r__;
    integer z__, ic, ii, kk;
    doublereal cs;
    integer is, iu;
    doublereal sn;
    integer nm1;
    doublereal eps;
    integer ivt, difl, difr, ierr, perm, mlvl, sqre;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern int dlasr_(char *, char *, char *, integer *, integer *, doublereal *, doublereal *,
                      doublereal *, integer *, ftnlen, ftnlen, ftnlen),
        dcopy_(integer *, doublereal *, integer *, doublereal *, integer *),
        dswap_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer poles, iuplo, nsize, start;
    extern int dlasd0_(integer *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                       doublereal *, integer *, integer *, integer *, doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern int dlasda_(integer *, integer *, integer *, integer *, doublereal *, doublereal *,
                       doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                       doublereal *, doublereal *, integer *, integer *, integer *, integer *,
                       doublereal *, doublereal *, doublereal *, doublereal *, integer *,
                       integer *),
        dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *,
                doublereal *, integer *, integer *, ftnlen),
        dlasdq_(char *, integer *, integer *, integer *, integer *, integer *, doublereal *,
                doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *,
                integer *, doublereal *, integer *, ftnlen),
        dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                ftnlen),
        dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *,
                           ftnlen, ftnlen);
    extern int xerbla_(char *, integer *, ftnlen);
    integer givcol;
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *, ftnlen);
    integer icompq;
    doublereal orgnrm;
    integer givnum, givptr, qstart, smlsiz, wstart, smlszp;
    --d__;
    --e;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --q;
    --iq;
    --work;
    --iwork;
    *info = 0;
    iuplo = 0;
    if (lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1)) {
        iuplo = 1;
    }
    if (lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        iuplo = 2;
    }
    if (lsame_(compq, (char *)"N", (ftnlen)1, (ftnlen)1)) {
        icompq = 0;
    } else if (lsame_(compq, (char *)"P", (ftnlen)1, (ftnlen)1)) {
        icompq = 1;
    } else if (lsame_(compq, (char *)"I", (ftnlen)1, (ftnlen)1)) {
        icompq = 2;
    } else {
        icompq = -1;
    }
    if (iuplo == 0) {
        *info = -1;
    } else if (icompq < 0) {
        *info = -2;
    } else if (*n < 0) {
        *info = -3;
    } else if (*ldu < 1 || icompq == 2 && *ldu < *n) {
        *info = -7;
    } else if (*ldvt < 1 || icompq == 2 && *ldvt < *n) {
        *info = -9;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DBDSDC", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n == 0) {
        return 0;
    }
    smlsiz = ilaenv_(&c__9, (char *)"DBDSDC", (char *)" ", &c__0, &c__0, &c__0, &c__0, (ftnlen)6, (ftnlen)1);
    if (*n == 1) {
        if (icompq == 1) {
            q[1] = d_lmp_sign(&c_b15, &d__[1]);
            q[smlsiz * *n + 1] = 1.;
        } else if (icompq == 2) {
            u[u_dim1 + 1] = d_lmp_sign(&c_b15, &d__[1]);
            vt[vt_dim1 + 1] = 1.;
        }
        d__[1] = abs(d__[1]);
        return 0;
    }
    nm1 = *n - 1;
    wstart = 1;
    qstart = 3;
    if (icompq == 1) {
        dcopy_(n, &d__[1], &c__1, &q[1], &c__1);
        i__1 = *n - 1;
        dcopy_(&i__1, &e[1], &c__1, &q[*n + 1], &c__1);
    }
    if (iuplo == 2) {
        qstart = 5;
        if (icompq == 2) {
            wstart = (*n << 1) - 1;
        }
        i__1 = *n - 1;
        for (i__ = 1; i__ <= i__1; ++i__) {
            dlartg_(&d__[i__], &e[i__], &cs, &sn, &r__);
            d__[i__] = r__;
            e[i__] = sn * d__[i__ + 1];
            d__[i__ + 1] = cs * d__[i__ + 1];
            if (icompq == 1) {
                q[i__ + (*n << 1)] = cs;
                q[i__ + *n * 3] = sn;
            } else if (icompq == 2) {
                work[i__] = cs;
                work[nm1 + i__] = -sn;
            }
        }
    }
    if (icompq == 0) {
        dlasdq_((char *)"U", &c__0, n, &c__0, &c__0, &c__0, &d__[1], &e[1], &vt[vt_offset], ldvt,
                &u[u_offset], ldu, &u[u_offset], ldu, &work[1], info, (ftnlen)1);
        goto L40;
    }
    if (*n <= smlsiz) {
        if (icompq == 2) {
            dlaset_((char *)"A", n, n, &c_b29, &c_b15, &u[u_offset], ldu, (ftnlen)1);
            dlaset_((char *)"A", n, n, &c_b29, &c_b15, &vt[vt_offset], ldvt, (ftnlen)1);
            dlasdq_((char *)"U", &c__0, n, n, n, &c__0, &d__[1], &e[1], &vt[vt_offset], ldvt, &u[u_offset],
                    ldu, &u[u_offset], ldu, &work[wstart], info, (ftnlen)1);
        } else if (icompq == 1) {
            iu = 1;
            ivt = iu + *n;
            dlaset_((char *)"A", n, n, &c_b29, &c_b15, &q[iu + (qstart - 1) * *n], n, (ftnlen)1);
            dlaset_((char *)"A", n, n, &c_b29, &c_b15, &q[ivt + (qstart - 1) * *n], n, (ftnlen)1);
            dlasdq_((char *)"U", &c__0, n, n, n, &c__0, &d__[1], &e[1], &q[ivt + (qstart - 1) * *n], n,
                    &q[iu + (qstart - 1) * *n], n, &q[iu + (qstart - 1) * *n], n, &work[wstart],
                    info, (ftnlen)1);
        }
        goto L40;
    }
    if (icompq == 2) {
        dlaset_((char *)"A", n, n, &c_b29, &c_b15, &u[u_offset], ldu, (ftnlen)1);
        dlaset_((char *)"A", n, n, &c_b29, &c_b15, &vt[vt_offset], ldvt, (ftnlen)1);
    }
    orgnrm = dlanst_((char *)"M", n, &d__[1], &e[1], (ftnlen)1);
    if (orgnrm == 0.) {
        return 0;
    }
    dlascl_((char *)"G", &c__0, &c__0, &orgnrm, &c_b15, n, &c__1, &d__[1], n, &ierr, (ftnlen)1);
    dlascl_((char *)"G", &c__0, &c__0, &orgnrm, &c_b15, &nm1, &c__1, &e[1], &nm1, &ierr, (ftnlen)1);
    eps = dlamch_((char *)"Epsilon", (ftnlen)7) * .9;
    mlvl = (integer)(log((doublereal)(*n) / (doublereal)(smlsiz + 1)) / log(2.)) + 1;
    smlszp = smlsiz + 1;
    if (icompq == 1) {
        iu = 1;
        ivt = smlsiz + 1;
        difl = ivt + smlszp;
        difr = difl + mlvl;
        z__ = difr + (mlvl << 1);
        ic = z__ + mlvl;
        is = ic + 1;
        poles = is + 1;
        givnum = poles + (mlvl << 1);
        k = 1;
        givptr = 2;
        perm = 3;
        givcol = perm + mlvl;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        if ((d__1 = d__[i__], abs(d__1)) < eps) {
            d__[i__] = d_lmp_sign(&eps, &d__[i__]);
        }
    }
    start = 1;
    sqre = 0;
    i__1 = nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
        if ((d__1 = e[i__], abs(d__1)) < eps || i__ == nm1) {
            if (i__ < nm1) {
                nsize = i__ - start + 1;
            } else if ((d__1 = e[i__], abs(d__1)) >= eps) {
                nsize = *n - start + 1;
            } else {
                nsize = i__ - start + 1;
                if (icompq == 2) {
                    u[*n + *n * u_dim1] = d_lmp_sign(&c_b15, &d__[*n]);
                    vt[*n + *n * vt_dim1] = 1.;
                } else if (icompq == 1) {
                    q[*n + (qstart - 1) * *n] = d_lmp_sign(&c_b15, &d__[*n]);
                    q[*n + (smlsiz + qstart - 1) * *n] = 1.;
                }
                d__[*n] = (d__1 = d__[*n], abs(d__1));
            }
            if (icompq == 2) {
                dlasd0_(&nsize, &sqre, &d__[start], &e[start], &u[start + start * u_dim1], ldu,
                        &vt[start + start * vt_dim1], ldvt, &smlsiz, &iwork[1], &work[wstart],
                        info);
            } else {
                dlasda_(&icompq, &smlsiz, &nsize, &sqre, &d__[start], &e[start],
                        &q[start + (iu + qstart - 2) * *n], n, &q[start + (ivt + qstart - 2) * *n],
                        &iq[start + k * *n], &q[start + (difl + qstart - 2) * *n],
                        &q[start + (difr + qstart - 2) * *n], &q[start + (z__ + qstart - 2) * *n],
                        &q[start + (poles + qstart - 2) * *n], &iq[start + givptr * *n],
                        &iq[start + givcol * *n], n, &iq[start + perm * *n],
                        &q[start + (givnum + qstart - 2) * *n], &q[start + (ic + qstart - 2) * *n],
                        &q[start + (is + qstart - 2) * *n], &work[wstart], &iwork[1], info);
            }
            if (*info != 0) {
                return 0;
            }
            start = i__ + 1;
        }
    }
    dlascl_((char *)"G", &c__0, &c__0, &c_b15, &orgnrm, n, &c__1, &d__[1], n, &ierr, (ftnlen)1);
L40:
    i__1 = *n;
    for (ii = 2; ii <= i__1; ++ii) {
        i__ = ii - 1;
        kk = i__;
        p = d__[i__];
        i__2 = *n;
        for (j = ii; j <= i__2; ++j) {
            if (d__[j] > p) {
                kk = j;
                p = d__[j];
            }
        }
        if (kk != i__) {
            d__[kk] = d__[i__];
            d__[i__] = p;
            if (icompq == 1) {
                iq[i__] = kk;
            } else if (icompq == 2) {
                dswap_(n, &u[i__ * u_dim1 + 1], &c__1, &u[kk * u_dim1 + 1], &c__1);
                dswap_(n, &vt[i__ + vt_dim1], ldvt, &vt[kk + vt_dim1], ldvt);
            }
        } else if (icompq == 1) {
            iq[i__] = i__;
        }
    }
    if (icompq == 1) {
        if (iuplo == 1) {
            iq[*n] = 1;
        } else {
            iq[*n] = 0;
        }
    }
    if (iuplo == 2 && icompq == 2) {
        dlasr_((char *)"L", (char *)"V", (char *)"B", n, n, &work[1], &work[*n], &u[u_offset], ldu, (ftnlen)1, (ftnlen)1,
               (ftnlen)1);
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
