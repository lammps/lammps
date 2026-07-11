#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
static integer c__2 = 2;
int dtrexc_(char *compq, integer *n, doublereal *t, integer *ldt, doublereal *q, integer *ldq,
            integer *ifst, integer *ilst, doublereal *work, integer *info, ftnlen compq_len)
{
    integer q_dim1, q_offset, t_dim1, t_offset, i__1;
    integer nbf, nbl, here;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    logical wantq;
    extern int dlaexc_(logical *, integer *, doublereal *, integer *, doublereal *, integer *,
                       integer *, integer *, integer *, doublereal *, integer *),
        xerbla_(char *, integer *, ftnlen);
    integer nbnext;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --work;
    *info = 0;
    wantq = lsame_(compq, (char *)"V", (ftnlen)1, (ftnlen)1);
    if (!wantq && !lsame_(compq, (char *)"N", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*ldt < max(1, *n)) {
        *info = -4;
    } else if (*ldq < 1 || wantq && *ldq < max(1, *n)) {
        *info = -6;
    } else if ((*ifst < 1 || *ifst > *n) && *n > 0) {
        *info = -7;
    } else if ((*ilst < 1 || *ilst > *n) && *n > 0) {
        *info = -8;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DTREXC", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n <= 1) {
        return 0;
    }
    if (*ifst > 1) {
        if (t[*ifst + (*ifst - 1) * t_dim1] != 0.) {
            --(*ifst);
        }
    }
    nbf = 1;
    if (*ifst < *n) {
        if (t[*ifst + 1 + *ifst * t_dim1] != 0.) {
            nbf = 2;
        }
    }
    if (*ilst > 1) {
        if (t[*ilst + (*ilst - 1) * t_dim1] != 0.) {
            --(*ilst);
        }
    }
    nbl = 1;
    if (*ilst < *n) {
        if (t[*ilst + 1 + *ilst * t_dim1] != 0.) {
            nbl = 2;
        }
    }
    if (*ifst == *ilst) {
        return 0;
    }
    if (*ifst < *ilst) {
        if (nbf == 2 && nbl == 1) {
            --(*ilst);
        }
        if (nbf == 1 && nbl == 2) {
            ++(*ilst);
        }
        here = *ifst;
    L10:
        if (nbf == 1 || nbf == 2) {
            nbnext = 1;
            if (here + nbf + 1 <= *n) {
                if (t[here + nbf + 1 + (here + nbf) * t_dim1] != 0.) {
                    nbnext = 2;
                }
            }
            dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &here, &nbf, &nbnext, &work[1],
                    info);
            if (*info != 0) {
                *ilst = here;
                return 0;
            }
            here += nbnext;
            if (nbf == 2) {
                if (t[here + 1 + here * t_dim1] == 0.) {
                    nbf = 3;
                }
            }
        } else {
            nbnext = 1;
            if (here + 3 <= *n) {
                if (t[here + 3 + (here + 2) * t_dim1] != 0.) {
                    nbnext = 2;
                }
            }
            i__1 = here + 1;
            dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &i__1, &c__1, &nbnext,
                    &work[1], info);
            if (*info != 0) {
                *ilst = here;
                return 0;
            }
            if (nbnext == 1) {
                dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &here, &c__1, &nbnext,
                        &work[1], info);
                ++here;
            } else {
                if (t[here + 2 + (here + 1) * t_dim1] == 0.) {
                    nbnext = 1;
                }
                if (nbnext == 2) {
                    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &here, &c__1, &nbnext,
                            &work[1], info);
                    if (*info != 0) {
                        *ilst = here;
                        return 0;
                    }
                    here += 2;
                } else {
                    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &here, &c__1, &c__1,
                            &work[1], info);
                    i__1 = here + 1;
                    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &i__1, &c__1, &c__1,
                            &work[1], info);
                    here += 2;
                }
            }
        }
        if (here < *ilst) {
            goto L10;
        }
    } else {
        here = *ifst;
    L20:
        if (nbf == 1 || nbf == 2) {
            nbnext = 1;
            if (here >= 3) {
                if (t[here - 1 + (here - 2) * t_dim1] != 0.) {
                    nbnext = 2;
                }
            }
            i__1 = here - nbnext;
            dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &i__1, &nbnext, &nbf, &work[1],
                    info);
            if (*info != 0) {
                *ilst = here;
                return 0;
            }
            here -= nbnext;
            if (nbf == 2) {
                if (t[here + 1 + here * t_dim1] == 0.) {
                    nbf = 3;
                }
            }
        } else {
            nbnext = 1;
            if (here >= 3) {
                if (t[here - 1 + (here - 2) * t_dim1] != 0.) {
                    nbnext = 2;
                }
            }
            i__1 = here - nbnext;
            dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &i__1, &nbnext, &c__1,
                    &work[1], info);
            if (*info != 0) {
                *ilst = here;
                return 0;
            }
            if (nbnext == 1) {
                dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &here, &nbnext, &c__1,
                        &work[1], info);
                --here;
            } else {
                if (t[here + (here - 1) * t_dim1] == 0.) {
                    nbnext = 1;
                }
                if (nbnext == 2) {
                    i__1 = here - 1;
                    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &i__1, &c__2, &c__1,
                            &work[1], info);
                    if (*info != 0) {
                        *ilst = here;
                        return 0;
                    }
                    here += -2;
                } else {
                    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &here, &c__1, &c__1,
                            &work[1], info);
                    i__1 = here - 1;
                    dlaexc_(&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &i__1, &c__1, &c__1,
                            &work[1], info);
                    here += -2;
                }
            }
        }
        if (here > *ilst) {
            goto L20;
        }
    }
    *ilst = here;
    return 0;
}
#ifdef __cplusplus
}
#endif
