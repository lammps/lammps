#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
integer dlaneg_(integer *n, doublereal *d__, doublereal *lld, doublereal *sigma, doublereal *pivmin,
                integer *r__)
{
    integer ret_val, i__1, i__2, i__3, i__4;
    integer j;
    doublereal p, t;
    integer bj;
    doublereal tmp;
    integer neg1, neg2;
    doublereal bsav, gamma, dplus;
    extern logical disnan_(doublereal *);
    integer negcnt;
    logical sawnan;
    doublereal dminus;
    --lld;
    --d__;
    negcnt = 0;
    t = -(*sigma);
    i__1 = *r__ - 1;
    for (bj = 1; bj <= i__1; bj += 128) {
        neg1 = 0;
        bsav = t;
        i__3 = bj + 127, i__4 = *r__ - 1;
        i__2 = min(i__3, i__4);
        for (j = bj; j <= i__2; ++j) {
            dplus = d__[j] + t;
            if (dplus < 0.) {
                ++neg1;
            }
            tmp = t / dplus;
            t = tmp * lld[j] - *sigma;
        }
        sawnan = disnan_(&t);
        if (sawnan) {
            neg1 = 0;
            t = bsav;
            i__3 = bj + 127, i__4 = *r__ - 1;
            i__2 = min(i__3, i__4);
            for (j = bj; j <= i__2; ++j) {
                dplus = d__[j] + t;
                if (dplus < 0.) {
                    ++neg1;
                }
                tmp = t / dplus;
                if (disnan_(&tmp)) {
                    tmp = 1.;
                }
                t = tmp * lld[j] - *sigma;
            }
        }
        negcnt += neg1;
    }
    p = d__[*n] - *sigma;
    i__1 = *r__;
    for (bj = *n - 1; bj >= i__1; bj += -128) {
        neg2 = 0;
        bsav = p;
        i__3 = bj - 127;
        i__2 = max(i__3, *r__);
        for (j = bj; j >= i__2; --j) {
            dminus = lld[j] + p;
            if (dminus < 0.) {
                ++neg2;
            }
            tmp = p / dminus;
            p = tmp * d__[j] - *sigma;
        }
        sawnan = disnan_(&p);
        if (sawnan) {
            neg2 = 0;
            p = bsav;
            i__3 = bj - 127;
            i__2 = max(i__3, *r__);
            for (j = bj; j >= i__2; --j) {
                dminus = lld[j] + p;
                if (dminus < 0.) {
                    ++neg2;
                }
                tmp = p / dminus;
                if (disnan_(&tmp)) {
                    tmp = 1.;
                }
                p = tmp * d__[j] - *sigma;
            }
        }
        negcnt += neg2;
    }
    gamma = t + *sigma + p;
    if (gamma < 0.) {
        ++negcnt;
    }
    ret_val = negcnt;
    return ret_val;
}
#ifdef __cplusplus
}
#endif
