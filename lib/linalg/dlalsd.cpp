#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static doublereal c_b6 = 0.;
static integer c__0 = 0;
static doublereal c_b11 = 1.;
int dlalsd_(char *uplo, integer *smlsiz, integer *n, integer *nrhs, doublereal *d__, doublereal *e,
            doublereal *b, integer *ldb, doublereal *rcond, integer *rank, doublereal *work,
            integer *iwork, integer *info, ftnlen uplo_len)
{
    integer b_dim1, b_offset, i__1, i__2;
    doublereal d__1;
    double log(doublereal), d_lmp_sign(doublereal *, doublereal *);
    integer c__, i__, j, k;
    doublereal r__;
    integer s, u, z__;
    doublereal cs;
    integer bx;
    doublereal sn;
    integer st, vt, nm1, st1;
    doublereal eps;
    integer iwk;
    doublereal tol;
    integer difl, difr;
    doublereal rcnd;
    integer perm, nsub;
    extern int drot_(integer *, doublereal *, integer *, doublereal *, integer *, doublereal *,
                     doublereal *);
    integer nlvl, sqre, bxst;
    extern int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
                      integer *, doublereal *, integer *, doublereal *, doublereal *, integer *,
                      ftnlen, ftnlen),
        dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer poles, sizei, nsize, nwork, icmpq1, icmpq2;
    extern doublereal dlamch_(char *, ftnlen);
    extern int dlasda_(integer *, integer *, integer *, integer *, doublereal *, doublereal *,
                       doublereal *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                       doublereal *, doublereal *, integer *, integer *, integer *, integer *,
                       doublereal *, doublereal *, doublereal *, doublereal *, integer *,
                       integer *),
        dlalsa_(integer *, integer *, integer *, integer *, doublereal *, integer *, doublereal *,
                integer *, doublereal *, integer *, doublereal *, integer *, doublereal *,
                doublereal *, doublereal *, doublereal *, integer *, integer *, integer *,
                integer *, doublereal *, doublereal *, doublereal *, doublereal *, integer *,
                integer *),
        dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *,
                doublereal *, integer *, integer *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern int dlasdq_(char *, integer *, integer *, integer *, integer *, integer *, doublereal *,
                       doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *,
                       integer *, doublereal *, integer *, ftnlen),
        dlacpy_(char *, integer *, integer *, doublereal *, integer *, doublereal *, integer *,
                ftnlen),
        dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *),
        dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                ftnlen),
        xerbla_(char *, integer *, ftnlen);
    integer givcol;
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *, ftnlen);
    extern int dlasrt_(char *, integer *, doublereal *, integer *, ftnlen);
    doublereal orgnrm;
    integer givnum, givptr, smlszp;
    --d__;
    --e;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --work;
    --iwork;
    *info = 0;
    if (*n < 0) {
        *info = -3;
    } else if (*nrhs < 1) {
        *info = -4;
    } else if (*ldb < 1 || *ldb < *n) {
        *info = -8;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLALSD", &i__1, (ftnlen)6);
        return 0;
    }
    eps = dlamch_((char *)"Epsilon", (ftnlen)7);
    if (*rcond <= 0. || *rcond >= 1.) {
        rcnd = eps;
    } else {
        rcnd = *rcond;
    }
    *rank = 0;
    if (*n == 0) {
        return 0;
    } else if (*n == 1) {
        if (d__[1] == 0.) {
            dlaset_((char *)"A", &c__1, nrhs, &c_b6, &c_b6, &b[b_offset], ldb, (ftnlen)1);
        } else {
            *rank = 1;
            dlascl_((char *)"G", &c__0, &c__0, &d__[1], &c_b11, &c__1, nrhs, &b[b_offset], ldb, info,
                    (ftnlen)1);
            d__[1] = abs(d__[1]);
        }
        return 0;
    }
    if (*(unsigned char *)uplo == 'L') {
        i__1 = *n - 1;
        for (i__ = 1; i__ <= i__1; ++i__) {
            dlartg_(&d__[i__], &e[i__], &cs, &sn, &r__);
            d__[i__] = r__;
            e[i__] = sn * d__[i__ + 1];
            d__[i__ + 1] = cs * d__[i__ + 1];
            if (*nrhs == 1) {
                drot_(&c__1, &b[i__ + b_dim1], &c__1, &b[i__ + 1 + b_dim1], &c__1, &cs, &sn);
            } else {
                work[(i__ << 1) - 1] = cs;
                work[i__ * 2] = sn;
            }
        }
        if (*nrhs > 1) {
            i__1 = *nrhs;
            for (i__ = 1; i__ <= i__1; ++i__) {
                i__2 = *n - 1;
                for (j = 1; j <= i__2; ++j) {
                    cs = work[(j << 1) - 1];
                    sn = work[j * 2];
                    drot_(&c__1, &b[j + i__ * b_dim1], &c__1, &b[j + 1 + i__ * b_dim1], &c__1, &cs,
                          &sn);
                }
            }
        }
    }
    nm1 = *n - 1;
    orgnrm = dlanst_((char *)"M", n, &d__[1], &e[1], (ftnlen)1);
    if (orgnrm == 0.) {
        dlaset_((char *)"A", n, nrhs, &c_b6, &c_b6, &b[b_offset], ldb, (ftnlen)1);
        return 0;
    }
    dlascl_((char *)"G", &c__0, &c__0, &orgnrm, &c_b11, n, &c__1, &d__[1], n, info, (ftnlen)1);
    dlascl_((char *)"G", &c__0, &c__0, &orgnrm, &c_b11, &nm1, &c__1, &e[1], &nm1, info, (ftnlen)1);
    if (*n <= *smlsiz) {
        nwork = *n * *n + 1;
        dlaset_((char *)"A", n, n, &c_b6, &c_b11, &work[1], n, (ftnlen)1);
        dlasdq_((char *)"U", &c__0, n, n, &c__0, nrhs, &d__[1], &e[1], &work[1], n, &work[1], n,
                &b[b_offset], ldb, &work[nwork], info, (ftnlen)1);
        if (*info != 0) {
            return 0;
        }
        tol = rcnd * (d__1 = d__[idamax_(n, &d__[1], &c__1)], abs(d__1));
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            if (d__[i__] <= tol) {
                dlaset_((char *)"A", &c__1, nrhs, &c_b6, &c_b6, &b[i__ + b_dim1], ldb, (ftnlen)1);
            } else {
                dlascl_((char *)"G", &c__0, &c__0, &d__[i__], &c_b11, &c__1, nrhs, &b[i__ + b_dim1], ldb,
                        info, (ftnlen)1);
                ++(*rank);
            }
        }
        dgemm_((char *)"T", (char *)"N", n, nrhs, n, &c_b11, &work[1], n, &b[b_offset], ldb, &c_b6, &work[nwork], n,
               (ftnlen)1, (ftnlen)1);
        dlacpy_((char *)"A", n, nrhs, &work[nwork], n, &b[b_offset], ldb, (ftnlen)1);
        dlascl_((char *)"G", &c__0, &c__0, &c_b11, &orgnrm, n, &c__1, &d__[1], n, info, (ftnlen)1);
        dlasrt_((char *)"D", n, &d__[1], info, (ftnlen)1);
        dlascl_((char *)"G", &c__0, &c__0, &orgnrm, &c_b11, n, nrhs, &b[b_offset], ldb, info, (ftnlen)1);
        return 0;
    }
    nlvl = (integer)(log((doublereal)(*n) / (doublereal)(*smlsiz + 1)) / log(2.)) + 1;
    smlszp = *smlsiz + 1;
    u = 1;
    vt = *smlsiz * *n + 1;
    difl = vt + smlszp * *n;
    difr = difl + nlvl * *n;
    z__ = difr + (nlvl * *n << 1);
    c__ = z__ + nlvl * *n;
    s = c__ + *n;
    poles = s + *n;
    givnum = poles + (nlvl << 1) * *n;
    bx = givnum + (nlvl << 1) * *n;
    nwork = bx + *n * *nrhs;
    sizei = *n + 1;
    k = sizei + *n;
    givptr = k + *n;
    perm = givptr + *n;
    givcol = perm + nlvl * *n;
    iwk = givcol + (nlvl * *n << 1);
    st = 1;
    sqre = 0;
    icmpq1 = 1;
    icmpq2 = 0;
    nsub = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        if ((d__1 = d__[i__], abs(d__1)) < eps) {
            d__[i__] = d_lmp_sign(&eps, &d__[i__]);
        }
    }
    i__1 = nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
        if ((d__1 = e[i__], abs(d__1)) < eps || i__ == nm1) {
            ++nsub;
            iwork[nsub] = st;
            if (i__ < nm1) {
                nsize = i__ - st + 1;
                iwork[sizei + nsub - 1] = nsize;
            } else if ((d__1 = e[i__], abs(d__1)) >= eps) {
                nsize = *n - st + 1;
                iwork[sizei + nsub - 1] = nsize;
            } else {
                nsize = i__ - st + 1;
                iwork[sizei + nsub - 1] = nsize;
                ++nsub;
                iwork[nsub] = *n;
                iwork[sizei + nsub - 1] = 1;
                dcopy_(nrhs, &b[*n + b_dim1], ldb, &work[bx + nm1], n);
            }
            st1 = st - 1;
            if (nsize == 1) {
                dcopy_(nrhs, &b[st + b_dim1], ldb, &work[bx + st1], n);
            } else if (nsize <= *smlsiz) {
                dlaset_((char *)"A", &nsize, &nsize, &c_b6, &c_b11, &work[vt + st1], n, (ftnlen)1);
                dlasdq_((char *)"U", &c__0, &nsize, &nsize, &c__0, nrhs, &d__[st], &e[st], &work[vt + st1],
                        n, &work[nwork], n, &b[st + b_dim1], ldb, &work[nwork], info, (ftnlen)1);
                if (*info != 0) {
                    return 0;
                }
                dlacpy_((char *)"A", &nsize, nrhs, &b[st + b_dim1], ldb, &work[bx + st1], n, (ftnlen)1);
            } else {
                dlasda_(&icmpq1, smlsiz, &nsize, &sqre, &d__[st], &e[st], &work[u + st1], n,
                        &work[vt + st1], &iwork[k + st1], &work[difl + st1], &work[difr + st1],
                        &work[z__ + st1], &work[poles + st1], &iwork[givptr + st1],
                        &iwork[givcol + st1], n, &iwork[perm + st1], &work[givnum + st1],
                        &work[c__ + st1], &work[s + st1], &work[nwork], &iwork[iwk], info);
                if (*info != 0) {
                    return 0;
                }
                bxst = bx + st1;
                dlalsa_(&icmpq2, smlsiz, &nsize, nrhs, &b[st + b_dim1], ldb, &work[bxst], n,
                        &work[u + st1], n, &work[vt + st1], &iwork[k + st1], &work[difl + st1],
                        &work[difr + st1], &work[z__ + st1], &work[poles + st1],
                        &iwork[givptr + st1], &iwork[givcol + st1], n, &iwork[perm + st1],
                        &work[givnum + st1], &work[c__ + st1], &work[s + st1], &work[nwork],
                        &iwork[iwk], info);
                if (*info != 0) {
                    return 0;
                }
            }
            st = i__ + 1;
        }
    }
    tol = rcnd * (d__1 = d__[idamax_(n, &d__[1], &c__1)], abs(d__1));
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        if ((d__1 = d__[i__], abs(d__1)) <= tol) {
            dlaset_((char *)"A", &c__1, nrhs, &c_b6, &c_b6, &work[bx + i__ - 1], n, (ftnlen)1);
        } else {
            ++(*rank);
            dlascl_((char *)"G", &c__0, &c__0, &d__[i__], &c_b11, &c__1, nrhs, &work[bx + i__ - 1], n, info,
                    (ftnlen)1);
        }
        d__[i__] = (d__1 = d__[i__], abs(d__1));
    }
    icmpq2 = 1;
    i__1 = nsub;
    for (i__ = 1; i__ <= i__1; ++i__) {
        st = iwork[i__];
        st1 = st - 1;
        nsize = iwork[sizei + i__ - 1];
        bxst = bx + st1;
        if (nsize == 1) {
            dcopy_(nrhs, &work[bxst], n, &b[st + b_dim1], ldb);
        } else if (nsize <= *smlsiz) {
            dgemm_((char *)"T", (char *)"N", &nsize, nrhs, &nsize, &c_b11, &work[vt + st1], n, &work[bxst], n,
                   &c_b6, &b[st + b_dim1], ldb, (ftnlen)1, (ftnlen)1);
        } else {
            dlalsa_(&icmpq2, smlsiz, &nsize, nrhs, &work[bxst], n, &b[st + b_dim1], ldb,
                    &work[u + st1], n, &work[vt + st1], &iwork[k + st1], &work[difl + st1],
                    &work[difr + st1], &work[z__ + st1], &work[poles + st1], &iwork[givptr + st1],
                    &iwork[givcol + st1], n, &iwork[perm + st1], &work[givnum + st1],
                    &work[c__ + st1], &work[s + st1], &work[nwork], &iwork[iwk], info);
            if (*info != 0) {
                return 0;
            }
        }
    }
    dlascl_((char *)"G", &c__0, &c__0, &c_b11, &orgnrm, n, &c__1, &d__[1], n, info, (ftnlen)1);
    dlasrt_((char *)"D", n, &d__[1], info, (ftnlen)1);
    dlascl_((char *)"G", &c__0, &c__0, &orgnrm, &c_b11, n, nrhs, &b[b_offset], ldb, info, (ftnlen)1);
    return 0;
}
#ifdef __cplusplus
}
#endif
