#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__0 = 0;
static doublereal c_b11 = 0.;
static doublereal c_b12 = 1.;
static integer c__1 = 1;
static integer c__2 = 2;
int dlasda_(integer *icompq, integer *smlsiz, integer *n, integer *sqre, doublereal *d__,
            doublereal *e, doublereal *u, integer *ldu, doublereal *vt, integer *k,
            doublereal *difl, doublereal *difr, doublereal *z__, doublereal *poles, integer *givptr,
            integer *givcol, integer *ldgcol, integer *perm, doublereal *givnum, doublereal *c__,
            doublereal *s, doublereal *work, integer *iwork, integer *info)
{
    integer givcol_dim1, givcol_offset, perm_dim1, perm_offset, difl_dim1, difl_offset, difr_dim1,
        difr_offset, givnum_dim1, givnum_offset, poles_dim1, poles_offset, u_dim1, u_offset,
        vt_dim1, vt_offset, z_dim1, z_offset, i__1, i__2;
    integer pow_lmp_ii(integer *, integer *);
    integer i__, j, m, i1, ic, lf, nd, ll, nl, vf, nr, vl, im1, ncc, nlf, nrf, vfi, iwk, vli, lvl,
        nru, ndb1, nlp1, lvl2, nrp1;
    doublereal beta;
    integer idxq, nlvl;
    doublereal alpha;
    integer inode, ndiml, ndimr, idxqi, itemp;
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    integer sqrei;
    extern int dlasd6_(integer *, integer *, integer *, integer *, doublereal *, doublereal *,
                       doublereal *, doublereal *, doublereal *, integer *, integer *, integer *,
                       integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                       doublereal *, doublereal *, integer *, doublereal *, doublereal *,
                       doublereal *, integer *, integer *);
    integer nwork1, nwork2;
    extern int dlasdq_(char *, integer *, integer *, integer *, integer *, integer *, doublereal *,
                       doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *,
                       integer *, doublereal *, integer *, ftnlen),
        dlasdt_(integer *, integer *, integer *, integer *, integer *, integer *, integer *),
        dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *,
                ftnlen),
        xerbla_(char *, integer *, ftnlen);
    integer smlszp;
    --d__;
    --e;
    givnum_dim1 = *ldu;
    givnum_offset = 1 + givnum_dim1;
    givnum -= givnum_offset;
    poles_dim1 = *ldu;
    poles_offset = 1 + poles_dim1;
    poles -= poles_offset;
    z_dim1 = *ldu;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    difr_dim1 = *ldu;
    difr_offset = 1 + difr_dim1;
    difr -= difr_offset;
    difl_dim1 = *ldu;
    difl_offset = 1 + difl_dim1;
    difl -= difl_offset;
    vt_dim1 = *ldu;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --k;
    --givptr;
    perm_dim1 = *ldgcol;
    perm_offset = 1 + perm_dim1;
    perm -= perm_offset;
    givcol_dim1 = *ldgcol;
    givcol_offset = 1 + givcol_dim1;
    givcol -= givcol_offset;
    --c__;
    --s;
    --work;
    --iwork;
    *info = 0;
    if (*icompq < 0 || *icompq > 1) {
        *info = -1;
    } else if (*smlsiz < 3) {
        *info = -2;
    } else if (*n < 0) {
        *info = -3;
    } else if (*sqre < 0 || *sqre > 1) {
        *info = -4;
    } else if (*ldu < *n + *sqre) {
        *info = -8;
    } else if (*ldgcol < *n) {
        *info = -17;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLASDA", &i__1, (ftnlen)6);
        return 0;
    }
    m = *n + *sqre;
    if (*n <= *smlsiz) {
        if (*icompq == 0) {
            dlasdq_((char *)"U", sqre, n, &c__0, &c__0, &c__0, &d__[1], &e[1], &vt[vt_offset], ldu,
                    &u[u_offset], ldu, &u[u_offset], ldu, &work[1], info, (ftnlen)1);
        } else {
            dlasdq_((char *)"U", sqre, n, &m, n, &c__0, &d__[1], &e[1], &vt[vt_offset], ldu, &u[u_offset],
                    ldu, &u[u_offset], ldu, &work[1], info, (ftnlen)1);
        }
        return 0;
    }
    inode = 1;
    ndiml = inode + *n;
    ndimr = ndiml + *n;
    idxq = ndimr + *n;
    iwk = idxq + *n;
    ncc = 0;
    nru = 0;
    smlszp = *smlsiz + 1;
    vf = 1;
    vl = vf + m;
    nwork1 = vl + m;
    nwork2 = nwork1 + smlszp * smlszp;
    dlasdt_(n, &nlvl, &nd, &iwork[inode], &iwork[ndiml], &iwork[ndimr], smlsiz);
    ndb1 = (nd + 1) / 2;
    i__1 = nd;
    for (i__ = ndb1; i__ <= i__1; ++i__) {
        i1 = i__ - 1;
        ic = iwork[inode + i1];
        nl = iwork[ndiml + i1];
        nlp1 = nl + 1;
        nr = iwork[ndimr + i1];
        nlf = ic - nl;
        nrf = ic + 1;
        idxqi = idxq + nlf - 2;
        vfi = vf + nlf - 1;
        vli = vl + nlf - 1;
        sqrei = 1;
        if (*icompq == 0) {
            dlaset_((char *)"A", &nlp1, &nlp1, &c_b11, &c_b12, &work[nwork1], &smlszp, (ftnlen)1);
            dlasdq_((char *)"U", &sqrei, &nl, &nlp1, &nru, &ncc, &d__[nlf], &e[nlf], &work[nwork1], &smlszp,
                    &work[nwork2], &nl, &work[nwork2], &nl, &work[nwork2], info, (ftnlen)1);
            itemp = nwork1 + nl * smlszp;
            dcopy_(&nlp1, &work[nwork1], &c__1, &work[vfi], &c__1);
            dcopy_(&nlp1, &work[itemp], &c__1, &work[vli], &c__1);
        } else {
            dlaset_((char *)"A", &nl, &nl, &c_b11, &c_b12, &u[nlf + u_dim1], ldu, (ftnlen)1);
            dlaset_((char *)"A", &nlp1, &nlp1, &c_b11, &c_b12, &vt[nlf + vt_dim1], ldu, (ftnlen)1);
            dlasdq_((char *)"U", &sqrei, &nl, &nlp1, &nl, &ncc, &d__[nlf], &e[nlf], &vt[nlf + vt_dim1], ldu,
                    &u[nlf + u_dim1], ldu, &u[nlf + u_dim1], ldu, &work[nwork1], info, (ftnlen)1);
            dcopy_(&nlp1, &vt[nlf + vt_dim1], &c__1, &work[vfi], &c__1);
            dcopy_(&nlp1, &vt[nlf + nlp1 * vt_dim1], &c__1, &work[vli], &c__1);
        }
        if (*info != 0) {
            return 0;
        }
        i__2 = nl;
        for (j = 1; j <= i__2; ++j) {
            iwork[idxqi + j] = j;
        }
        if (i__ == nd && *sqre == 0) {
            sqrei = 0;
        } else {
            sqrei = 1;
        }
        idxqi += nlp1;
        vfi += nlp1;
        vli += nlp1;
        nrp1 = nr + sqrei;
        if (*icompq == 0) {
            dlaset_((char *)"A", &nrp1, &nrp1, &c_b11, &c_b12, &work[nwork1], &smlszp, (ftnlen)1);
            dlasdq_((char *)"U", &sqrei, &nr, &nrp1, &nru, &ncc, &d__[nrf], &e[nrf], &work[nwork1], &smlszp,
                    &work[nwork2], &nr, &work[nwork2], &nr, &work[nwork2], info, (ftnlen)1);
            itemp = nwork1 + (nrp1 - 1) * smlszp;
            dcopy_(&nrp1, &work[nwork1], &c__1, &work[vfi], &c__1);
            dcopy_(&nrp1, &work[itemp], &c__1, &work[vli], &c__1);
        } else {
            dlaset_((char *)"A", &nr, &nr, &c_b11, &c_b12, &u[nrf + u_dim1], ldu, (ftnlen)1);
            dlaset_((char *)"A", &nrp1, &nrp1, &c_b11, &c_b12, &vt[nrf + vt_dim1], ldu, (ftnlen)1);
            dlasdq_((char *)"U", &sqrei, &nr, &nrp1, &nr, &ncc, &d__[nrf], &e[nrf], &vt[nrf + vt_dim1], ldu,
                    &u[nrf + u_dim1], ldu, &u[nrf + u_dim1], ldu, &work[nwork1], info, (ftnlen)1);
            dcopy_(&nrp1, &vt[nrf + vt_dim1], &c__1, &work[vfi], &c__1);
            dcopy_(&nrp1, &vt[nrf + nrp1 * vt_dim1], &c__1, &work[vli], &c__1);
        }
        if (*info != 0) {
            return 0;
        }
        i__2 = nr;
        for (j = 1; j <= i__2; ++j) {
            iwork[idxqi + j] = j;
        }
    }
    j = pow_lmp_ii(&c__2, &nlvl);
    for (lvl = nlvl; lvl >= 1; --lvl) {
        lvl2 = (lvl << 1) - 1;
        if (lvl == 1) {
            lf = 1;
            ll = 1;
        } else {
            i__1 = lvl - 1;
            lf = pow_lmp_ii(&c__2, &i__1);
            ll = (lf << 1) - 1;
        }
        i__1 = ll;
        for (i__ = lf; i__ <= i__1; ++i__) {
            im1 = i__ - 1;
            ic = iwork[inode + im1];
            nl = iwork[ndiml + im1];
            nr = iwork[ndimr + im1];
            nlf = ic - nl;
            nrf = ic + 1;
            if (i__ == ll) {
                sqrei = *sqre;
            } else {
                sqrei = 1;
            }
            vfi = vf + nlf - 1;
            vli = vl + nlf - 1;
            idxqi = idxq + nlf - 1;
            alpha = d__[ic];
            beta = e[ic];
            if (*icompq == 0) {
                dlasd6_(icompq, &nl, &nr, &sqrei, &d__[nlf], &work[vfi], &work[vli], &alpha, &beta,
                        &iwork[idxqi], &perm[perm_offset], &givptr[1], &givcol[givcol_offset],
                        ldgcol, &givnum[givnum_offset], ldu, &poles[poles_offset],
                        &difl[difl_offset], &difr[difr_offset], &z__[z_offset], &k[1], &c__[1],
                        &s[1], &work[nwork1], &iwork[iwk], info);
            } else {
                --j;
                dlasd6_(icompq, &nl, &nr, &sqrei, &d__[nlf], &work[vfi], &work[vli], &alpha, &beta,
                        &iwork[idxqi], &perm[nlf + lvl * perm_dim1], &givptr[j],
                        &givcol[nlf + lvl2 * givcol_dim1], ldgcol,
                        &givnum[nlf + lvl2 * givnum_dim1], ldu, &poles[nlf + lvl2 * poles_dim1],
                        &difl[nlf + lvl * difl_dim1], &difr[nlf + lvl2 * difr_dim1],
                        &z__[nlf + lvl * z_dim1], &k[j], &c__[j], &s[j], &work[nwork1], &iwork[iwk],
                        info);
            }
            if (*info != 0) {
                return 0;
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
