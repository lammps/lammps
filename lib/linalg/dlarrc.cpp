#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlarrc_(char *jobt, integer *n, doublereal *vl, doublereal *vu, doublereal *d__, doublereal *e,
            doublereal *pivmin, integer *eigcnt, integer *lcnt, integer *rcnt, integer *info,
            ftnlen jobt_len)
{
    integer i__1;
    doublereal d__1;
    integer i__;
    doublereal sl, su, tmp, tmp2;
    logical matt;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    doublereal lpivot, rpivot;
    --e;
    --d__;
    *info = 0;
    *lcnt = 0;
    *rcnt = 0;
    *eigcnt = 0;
    if (*n <= 0) {
        return 0;
    }
    matt = lsame_(jobt, (char *)"T", (ftnlen)1, (ftnlen)1);
    if (matt) {
        lpivot = d__[1] - *vl;
        rpivot = d__[1] - *vu;
        if (lpivot <= 0.) {
            ++(*lcnt);
        }
        if (rpivot <= 0.) {
            ++(*rcnt);
        }
        i__1 = *n - 1;
        for (i__ = 1; i__ <= i__1; ++i__) {
            d__1 = e[i__];
            tmp = d__1 * d__1;
            lpivot = d__[i__ + 1] - *vl - tmp / lpivot;
            rpivot = d__[i__ + 1] - *vu - tmp / rpivot;
            if (lpivot <= 0.) {
                ++(*lcnt);
            }
            if (rpivot <= 0.) {
                ++(*rcnt);
            }
        }
    } else {
        sl = -(*vl);
        su = -(*vu);
        i__1 = *n - 1;
        for (i__ = 1; i__ <= i__1; ++i__) {
            lpivot = d__[i__] + sl;
            rpivot = d__[i__] + su;
            if (lpivot <= 0.) {
                ++(*lcnt);
            }
            if (rpivot <= 0.) {
                ++(*rcnt);
            }
            tmp = e[i__] * d__[i__] * e[i__];
            tmp2 = tmp / lpivot;
            if (tmp2 == 0.) {
                sl = tmp - *vl;
            } else {
                sl = sl * tmp2 - *vl;
            }
            tmp2 = tmp / rpivot;
            if (tmp2 == 0.) {
                su = tmp - *vu;
            } else {
                su = su * tmp2 - *vu;
            }
        }
        lpivot = d__[*n] + sl;
        rpivot = d__[*n] + su;
        if (lpivot <= 0.) {
            ++(*lcnt);
        }
        if (rpivot <= 0.) {
            ++(*rcnt);
        }
    }
    *eigcnt = *rcnt - *lcnt;
    return 0;
}
#ifdef __cplusplus
}
#endif
