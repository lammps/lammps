#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlaebz_(integer *ijob, integer *nitmax, integer *n, integer *mmax, integer *minp,
            integer *nbmin, doublereal *abstol, doublereal *reltol, doublereal *pivmin,
            doublereal *d__, doublereal *e, doublereal *e2, integer *nval, doublereal *ab,
            doublereal *c__, integer *mout, integer *nab, doublereal *work, integer *iwork,
            integer *info)
{
    integer nab_dim1, nab_offset, ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2, d__3, d__4;
    integer j, kf, ji, kl, jp, jit;
    doublereal tmp1, tmp2;
    integer itmp1, itmp2, kfnew, klnew;
    nab_dim1 = *mmax;
    nab_offset = 1 + nab_dim1;
    nab -= nab_offset;
    ab_dim1 = *mmax;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --d__;
    --e;
    --e2;
    --nval;
    --c__;
    --work;
    --iwork;
    *info = 0;
    if (*ijob < 1 || *ijob > 3) {
        *info = -1;
        return 0;
    }
    if (*ijob == 1) {
        *mout = 0;
        i__1 = *minp;
        for (ji = 1; ji <= i__1; ++ji) {
            for (jp = 1; jp <= 2; ++jp) {
                tmp1 = d__[1] - ab[ji + jp * ab_dim1];
                if (abs(tmp1) < *pivmin) {
                    tmp1 = -(*pivmin);
                }
                nab[ji + jp * nab_dim1] = 0;
                if (tmp1 <= 0.) {
                    nab[ji + jp * nab_dim1] = 1;
                }
                i__2 = *n;
                for (j = 2; j <= i__2; ++j) {
                    tmp1 = d__[j] - e2[j - 1] / tmp1 - ab[ji + jp * ab_dim1];
                    if (abs(tmp1) < *pivmin) {
                        tmp1 = -(*pivmin);
                    }
                    if (tmp1 <= 0.) {
                        ++nab[ji + jp * nab_dim1];
                    }
                }
            }
            *mout = *mout + nab[ji + (nab_dim1 << 1)] - nab[ji + nab_dim1];
        }
        return 0;
    }
    kf = 1;
    kl = *minp;
    if (*ijob == 2) {
        i__1 = *minp;
        for (ji = 1; ji <= i__1; ++ji) {
            c__[ji] = (ab[ji + ab_dim1] + ab[ji + (ab_dim1 << 1)]) * .5;
        }
    }
    i__1 = *nitmax;
    for (jit = 1; jit <= i__1; ++jit) {
        if (kl - kf + 1 >= *nbmin && *nbmin > 0) {
            i__2 = kl;
            for (ji = kf; ji <= i__2; ++ji) {
                work[ji] = d__[1] - c__[ji];
                iwork[ji] = 0;
                if (work[ji] <= *pivmin) {
                    iwork[ji] = 1;
                    d__1 = work[ji], d__2 = -(*pivmin);
                    work[ji] = min(d__1, d__2);
                }
                i__3 = *n;
                for (j = 2; j <= i__3; ++j) {
                    work[ji] = d__[j] - e2[j - 1] / work[ji] - c__[ji];
                    if (work[ji] <= *pivmin) {
                        ++iwork[ji];
                        d__1 = work[ji], d__2 = -(*pivmin);
                        work[ji] = min(d__1, d__2);
                    }
                }
            }
            if (*ijob <= 2) {
                klnew = kl;
                i__2 = kl;
                for (ji = kf; ji <= i__2; ++ji) {
                    i__5 = nab[ji + nab_dim1], i__6 = iwork[ji];
                    i__3 = nab[ji + (nab_dim1 << 1)], i__4 = max(i__5, i__6);
                    iwork[ji] = min(i__3, i__4);
                    if (iwork[ji] == nab[ji + (nab_dim1 << 1)]) {
                        ab[ji + (ab_dim1 << 1)] = c__[ji];
                    } else if (iwork[ji] == nab[ji + nab_dim1]) {
                        ab[ji + ab_dim1] = c__[ji];
                    } else {
                        ++klnew;
                        if (klnew <= *mmax) {
                            ab[klnew + (ab_dim1 << 1)] = ab[ji + (ab_dim1 << 1)];
                            nab[klnew + (nab_dim1 << 1)] = nab[ji + (nab_dim1 << 1)];
                            ab[klnew + ab_dim1] = c__[ji];
                            nab[klnew + nab_dim1] = iwork[ji];
                            ab[ji + (ab_dim1 << 1)] = c__[ji];
                            nab[ji + (nab_dim1 << 1)] = iwork[ji];
                        } else {
                            *info = *mmax + 1;
                        }
                    }
                }
                if (*info != 0) {
                    return 0;
                }
                kl = klnew;
            } else {
                i__2 = kl;
                for (ji = kf; ji <= i__2; ++ji) {
                    if (iwork[ji] <= nval[ji]) {
                        ab[ji + ab_dim1] = c__[ji];
                        nab[ji + nab_dim1] = iwork[ji];
                    }
                    if (iwork[ji] >= nval[ji]) {
                        ab[ji + (ab_dim1 << 1)] = c__[ji];
                        nab[ji + (nab_dim1 << 1)] = iwork[ji];
                    }
                }
            }
        } else {
            klnew = kl;
            i__2 = kl;
            for (ji = kf; ji <= i__2; ++ji) {
                tmp1 = c__[ji];
                tmp2 = d__[1] - tmp1;
                itmp1 = 0;
                if (tmp2 <= *pivmin) {
                    itmp1 = 1;
                    d__1 = tmp2, d__2 = -(*pivmin);
                    tmp2 = min(d__1, d__2);
                }
                i__3 = *n;
                for (j = 2; j <= i__3; ++j) {
                    tmp2 = d__[j] - e2[j - 1] / tmp2 - tmp1;
                    if (tmp2 <= *pivmin) {
                        ++itmp1;
                        d__1 = tmp2, d__2 = -(*pivmin);
                        tmp2 = min(d__1, d__2);
                    }
                }
                if (*ijob <= 2) {
                    i__5 = nab[ji + nab_dim1];
                    i__3 = nab[ji + (nab_dim1 << 1)], i__4 = max(i__5, itmp1);
                    itmp1 = min(i__3, i__4);
                    if (itmp1 == nab[ji + (nab_dim1 << 1)]) {
                        ab[ji + (ab_dim1 << 1)] = tmp1;
                    } else if (itmp1 == nab[ji + nab_dim1]) {
                        ab[ji + ab_dim1] = tmp1;
                    } else if (klnew < *mmax) {
                        ++klnew;
                        ab[klnew + (ab_dim1 << 1)] = ab[ji + (ab_dim1 << 1)];
                        nab[klnew + (nab_dim1 << 1)] = nab[ji + (nab_dim1 << 1)];
                        ab[klnew + ab_dim1] = tmp1;
                        nab[klnew + nab_dim1] = itmp1;
                        ab[ji + (ab_dim1 << 1)] = tmp1;
                        nab[ji + (nab_dim1 << 1)] = itmp1;
                    } else {
                        *info = *mmax + 1;
                        return 0;
                    }
                } else {
                    if (itmp1 <= nval[ji]) {
                        ab[ji + ab_dim1] = tmp1;
                        nab[ji + nab_dim1] = itmp1;
                    }
                    if (itmp1 >= nval[ji]) {
                        ab[ji + (ab_dim1 << 1)] = tmp1;
                        nab[ji + (nab_dim1 << 1)] = itmp1;
                    }
                }
            }
            kl = klnew;
        }
        kfnew = kf;
        i__2 = kl;
        for (ji = kf; ji <= i__2; ++ji) {
            tmp1 = (d__1 = ab[ji + (ab_dim1 << 1)] - ab[ji + ab_dim1], abs(d__1));
            d__3 = (d__1 = ab[ji + (ab_dim1 << 1)], abs(d__1)),
            d__4 = (d__2 = ab[ji + ab_dim1], abs(d__2));
            tmp2 = max(d__3, d__4);
            d__1 = max(*abstol, *pivmin), d__2 = *reltol * tmp2;
            if (tmp1 < max(d__1, d__2) || nab[ji + nab_dim1] >= nab[ji + (nab_dim1 << 1)]) {
                if (ji > kfnew) {
                    tmp1 = ab[ji + ab_dim1];
                    tmp2 = ab[ji + (ab_dim1 << 1)];
                    itmp1 = nab[ji + nab_dim1];
                    itmp2 = nab[ji + (nab_dim1 << 1)];
                    ab[ji + ab_dim1] = ab[kfnew + ab_dim1];
                    ab[ji + (ab_dim1 << 1)] = ab[kfnew + (ab_dim1 << 1)];
                    nab[ji + nab_dim1] = nab[kfnew + nab_dim1];
                    nab[ji + (nab_dim1 << 1)] = nab[kfnew + (nab_dim1 << 1)];
                    ab[kfnew + ab_dim1] = tmp1;
                    ab[kfnew + (ab_dim1 << 1)] = tmp2;
                    nab[kfnew + nab_dim1] = itmp1;
                    nab[kfnew + (nab_dim1 << 1)] = itmp2;
                    if (*ijob == 3) {
                        itmp1 = nval[ji];
                        nval[ji] = nval[kfnew];
                        nval[kfnew] = itmp1;
                    }
                }
                ++kfnew;
            }
        }
        kf = kfnew;
        i__2 = kl;
        for (ji = kf; ji <= i__2; ++ji) {
            c__[ji] = (ab[ji + ab_dim1] + ab[ji + (ab_dim1 << 1)]) * .5;
        }
        if (kf > kl) {
            goto L140;
        }
    }
L140:
    i__1 = kl + 1 - kf;
    *info = max(i__1, 0);
    *mout = kl;
    return 0;
}
#ifdef __cplusplus
}
#endif
