#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__0 = 0;
static integer c__2 = 2;
int dlasd0_(integer *n, integer *sqre, doublereal *d__, doublereal *e, doublereal *u, integer *ldu,
            doublereal *vt, integer *ldvt, integer *smlsiz, integer *iwork, doublereal *work,
            integer *info)
{
    integer u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2;
    integer pow_lmp_ii(integer *, integer *);
    integer i__, j, m, i1, ic, lf, nd, ll, nl, nr, im1, ncc, nlf, nrf, iwk, lvl, ndb1, nlp1, nrp1;
    doublereal beta;
    integer idxq, nlvl;
    doublereal alpha;
    integer inode, ndiml, idxqc, ndimr, itemp, sqrei;
    extern int dlasd1_(integer *, integer *, integer *, doublereal *, doublereal *, doublereal *,
                       doublereal *, integer *, doublereal *, integer *, integer *, integer *,
                       doublereal *, integer *),
        dlasdq_(char *, integer *, integer *, integer *, integer *, integer *, doublereal *,
                doublereal *, doublereal *, integer *, doublereal *, integer *, doublereal *,
                integer *, doublereal *, integer *, ftnlen),
        dlasdt_(integer *, integer *, integer *, integer *, integer *, integer *, integer *),
        xerbla_(char *, integer *, ftnlen);
    --d__;
    --e;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --iwork;
    --work;
    *info = 0;
    if (*n < 0) {
        *info = -1;
    } else if (*sqre < 0 || *sqre > 1) {
        *info = -2;
    }
    m = *n + *sqre;
    if (*ldu < *n) {
        *info = -6;
    } else if (*ldvt < m) {
        *info = -8;
    } else if (*smlsiz < 3) {
        *info = -9;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLASD0", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n <= *smlsiz) {
        dlasdq_((char *)"U", sqre, n, &m, n, &c__0, &d__[1], &e[1], &vt[vt_offset], ldvt, &u[u_offset], ldu,
                &u[u_offset], ldu, &work[1], info, (ftnlen)1);
        return 0;
    }
    inode = 1;
    ndiml = inode + *n;
    ndimr = ndiml + *n;
    idxq = ndimr + *n;
    iwk = idxq + *n;
    dlasdt_(n, &nlvl, &nd, &iwork[inode], &iwork[ndiml], &iwork[ndimr], smlsiz);
    ndb1 = (nd + 1) / 2;
    ncc = 0;
    i__1 = nd;
    for (i__ = ndb1; i__ <= i__1; ++i__) {
        i1 = i__ - 1;
        ic = iwork[inode + i1];
        nl = iwork[ndiml + i1];
        nlp1 = nl + 1;
        nr = iwork[ndimr + i1];
        nrp1 = nr + 1;
        nlf = ic - nl;
        nrf = ic + 1;
        sqrei = 1;
        dlasdq_((char *)"U", &sqrei, &nl, &nlp1, &nl, &ncc, &d__[nlf], &e[nlf], &vt[nlf + nlf * vt_dim1],
                ldvt, &u[nlf + nlf * u_dim1], ldu, &u[nlf + nlf * u_dim1], ldu, &work[1], info,
                (ftnlen)1);
        if (*info != 0) {
            return 0;
        }
        itemp = idxq + nlf - 2;
        i__2 = nl;
        for (j = 1; j <= i__2; ++j) {
            iwork[itemp + j] = j;
        }
        if (i__ == nd) {
            sqrei = *sqre;
        } else {
            sqrei = 1;
        }
        nrp1 = nr + sqrei;
        dlasdq_((char *)"U", &sqrei, &nr, &nrp1, &nr, &ncc, &d__[nrf], &e[nrf], &vt[nrf + nrf * vt_dim1],
                ldvt, &u[nrf + nrf * u_dim1], ldu, &u[nrf + nrf * u_dim1], ldu, &work[1], info,
                (ftnlen)1);
        if (*info != 0) {
            return 0;
        }
        itemp = idxq + ic;
        i__2 = nr;
        for (j = 1; j <= i__2; ++j) {
            iwork[itemp + j - 1] = j;
        }
    }
    for (lvl = nlvl; lvl >= 1; --lvl) {
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
            if (*sqre == 0 && i__ == ll) {
                sqrei = *sqre;
            } else {
                sqrei = 1;
            }
            idxqc = idxq + nlf - 1;
            alpha = d__[ic];
            beta = e[ic];
            dlasd1_(&nl, &nr, &sqrei, &d__[nlf], &alpha, &beta, &u[nlf + nlf * u_dim1], ldu,
                    &vt[nlf + nlf * vt_dim1], ldvt, &iwork[idxqc], &iwork[iwk], &work[1], info);
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
