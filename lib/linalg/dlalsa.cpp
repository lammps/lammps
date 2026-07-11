#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static doublereal c_b7 = 1.;
static doublereal c_b8 = 0.;
static integer c__2 = 2;
int dlalsa_(integer *icompq, integer *smlsiz, integer *n, integer *nrhs, doublereal *b,
            integer *ldb, doublereal *bx, integer *ldbx, doublereal *u, integer *ldu,
            doublereal *vt, integer *k, doublereal *difl, doublereal *difr, doublereal *z__,
            doublereal *poles, integer *givptr, integer *givcol, integer *ldgcol, integer *perm,
            doublereal *givnum, doublereal *c__, doublereal *s, doublereal *work, integer *iwork,
            integer *info)
{
    integer givcol_dim1, givcol_offset, perm_dim1, perm_offset, b_dim1, b_offset, bx_dim1,
        bx_offset, difl_dim1, difl_offset, difr_dim1, difr_offset, givnum_dim1, givnum_offset,
        poles_dim1, poles_offset, u_dim1, u_offset, vt_dim1, vt_offset, z_dim1, z_offset, i__1,
        i__2;
    integer pow_lmp_ii(integer *, integer *);
    integer i__, j, i1, ic, lf, nd, ll, nl, nr, im1, nlf, nrf, lvl, ndb1, nlp1, lvl2, nrp1, nlvl,
        sqre;
    extern int dgemm_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
                      integer *, doublereal *, integer *, doublereal *, doublereal *, integer *,
                      ftnlen, ftnlen);
    integer inode, ndiml, ndimr;
    extern int dcopy_(integer *, doublereal *, integer *, doublereal *, integer *),
        dlals0_(integer *, integer *, integer *, integer *, integer *, doublereal *, integer *,
                doublereal *, integer *, integer *, integer *, integer *, integer *, doublereal *,
                integer *, doublereal *, doublereal *, doublereal *, doublereal *, integer *,
                doublereal *, doublereal *, doublereal *, integer *),
        dlasdt_(integer *, integer *, integer *, integer *, integer *, integer *, integer *),
        xerbla_(char *, integer *, ftnlen);
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    bx_dim1 = *ldbx;
    bx_offset = 1 + bx_dim1;
    bx -= bx_offset;
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
    } else if (*n < *smlsiz) {
        *info = -3;
    } else if (*nrhs < 1) {
        *info = -4;
    } else if (*ldb < *n) {
        *info = -6;
    } else if (*ldbx < *n) {
        *info = -8;
    } else if (*ldu < *n) {
        *info = -10;
    } else if (*ldgcol < *n) {
        *info = -19;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLALSA", &i__1, (ftnlen)6);
        return 0;
    }
    inode = 1;
    ndiml = inode + *n;
    ndimr = ndiml + *n;
    dlasdt_(n, &nlvl, &nd, &iwork[inode], &iwork[ndiml], &iwork[ndimr], smlsiz);
    if (*icompq == 1) {
        goto L50;
    }
    ndb1 = (nd + 1) / 2;
    i__1 = nd;
    for (i__ = ndb1; i__ <= i__1; ++i__) {
        i1 = i__ - 1;
        ic = iwork[inode + i1];
        nl = iwork[ndiml + i1];
        nr = iwork[ndimr + i1];
        nlf = ic - nl;
        nrf = ic + 1;
        dgemm_((char *)"T", (char *)"N", &nl, nrhs, &nl, &c_b7, &u[nlf + u_dim1], ldu, &b[nlf + b_dim1], ldb, &c_b8,
               &bx[nlf + bx_dim1], ldbx, (ftnlen)1, (ftnlen)1);
        dgemm_((char *)"T", (char *)"N", &nr, nrhs, &nr, &c_b7, &u[nrf + u_dim1], ldu, &b[nrf + b_dim1], ldb, &c_b8,
               &bx[nrf + bx_dim1], ldbx, (ftnlen)1, (ftnlen)1);
    }
    i__1 = nd;
    for (i__ = 1; i__ <= i__1; ++i__) {
        ic = iwork[inode + i__ - 1];
        dcopy_(nrhs, &b[ic + b_dim1], ldb, &bx[ic + bx_dim1], ldbx);
    }
    j = pow_lmp_ii(&c__2, &nlvl);
    sqre = 0;
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
            --j;
            dlals0_(icompq, &nl, &nr, &sqre, nrhs, &bx[nlf + bx_dim1], ldbx, &b[nlf + b_dim1], ldb,
                    &perm[nlf + lvl * perm_dim1], &givptr[j], &givcol[nlf + lvl2 * givcol_dim1],
                    ldgcol, &givnum[nlf + lvl2 * givnum_dim1], ldu, &poles[nlf + lvl2 * poles_dim1],
                    &difl[nlf + lvl * difl_dim1], &difr[nlf + lvl2 * difr_dim1],
                    &z__[nlf + lvl * z_dim1], &k[j], &c__[j], &s[j], &work[1], info);
        }
    }
    goto L90;
L50:
    j = 0;
    i__1 = nlvl;
    for (lvl = 1; lvl <= i__1; ++lvl) {
        lvl2 = (lvl << 1) - 1;
        if (lvl == 1) {
            lf = 1;
            ll = 1;
        } else {
            i__2 = lvl - 1;
            lf = pow_lmp_ii(&c__2, &i__2);
            ll = (lf << 1) - 1;
        }
        i__2 = lf;
        for (i__ = ll; i__ >= i__2; --i__) {
            im1 = i__ - 1;
            ic = iwork[inode + im1];
            nl = iwork[ndiml + im1];
            nr = iwork[ndimr + im1];
            nlf = ic - nl;
            nrf = ic + 1;
            if (i__ == ll) {
                sqre = 0;
            } else {
                sqre = 1;
            }
            ++j;
            dlals0_(icompq, &nl, &nr, &sqre, nrhs, &b[nlf + b_dim1], ldb, &bx[nlf + bx_dim1], ldbx,
                    &perm[nlf + lvl * perm_dim1], &givptr[j], &givcol[nlf + lvl2 * givcol_dim1],
                    ldgcol, &givnum[nlf + lvl2 * givnum_dim1], ldu, &poles[nlf + lvl2 * poles_dim1],
                    &difl[nlf + lvl * difl_dim1], &difr[nlf + lvl2 * difr_dim1],
                    &z__[nlf + lvl * z_dim1], &k[j], &c__[j], &s[j], &work[1], info);
        }
    }
    ndb1 = (nd + 1) / 2;
    i__1 = nd;
    for (i__ = ndb1; i__ <= i__1; ++i__) {
        i1 = i__ - 1;
        ic = iwork[inode + i1];
        nl = iwork[ndiml + i1];
        nr = iwork[ndimr + i1];
        nlp1 = nl + 1;
        if (i__ == nd) {
            nrp1 = nr;
        } else {
            nrp1 = nr + 1;
        }
        nlf = ic - nl;
        nrf = ic + 1;
        dgemm_((char *)"T", (char *)"N", &nlp1, nrhs, &nlp1, &c_b7, &vt[nlf + vt_dim1], ldu, &b[nlf + b_dim1], ldb,
               &c_b8, &bx[nlf + bx_dim1], ldbx, (ftnlen)1, (ftnlen)1);
        dgemm_((char *)"T", (char *)"N", &nrp1, nrhs, &nrp1, &c_b7, &vt[nrf + vt_dim1], ldu, &b[nrf + b_dim1], ldb,
               &c_b8, &bx[nrf + bx_dim1], ldbx, (ftnlen)1, (ftnlen)1);
    }
L90:
    return 0;
}
#ifdef __cplusplus
}
#endif
