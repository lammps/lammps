#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlasrt_(char *id, integer *n, doublereal *d__, integer *info, ftnlen id_len)
{
    integer i__1, i__2;
    integer i__, j;
    doublereal d1, d2, d3;
    integer dir;
    doublereal tmp;
    integer endd;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer stack[64];
    doublereal dmnmx;
    integer start;
    extern int xerbla_(char *, integer *, ftnlen);
    integer stkpnt;
    --d__;
    *info = 0;
    dir = -1;
    if (lsame_(id, (char *)"D", (ftnlen)1, (ftnlen)1)) {
        dir = 0;
    } else if (lsame_(id, (char *)"I", (ftnlen)1, (ftnlen)1)) {
        dir = 1;
    }
    if (dir == -1) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DLASRT", &i__1, (ftnlen)6);
        return 0;
    }
    if (*n <= 1) {
        return 0;
    }
    stkpnt = 1;
    stack[0] = 1;
    stack[1] = *n;
L10:
    start = stack[(stkpnt << 1) - 2];
    endd = stack[(stkpnt << 1) - 1];
    --stkpnt;
    if (endd - start <= 20 && endd - start > 0) {
        if (dir == 0) {
            i__1 = endd;
            for (i__ = start + 1; i__ <= i__1; ++i__) {
                i__2 = start + 1;
                for (j = i__; j >= i__2; --j) {
                    if (d__[j] > d__[j - 1]) {
                        dmnmx = d__[j];
                        d__[j] = d__[j - 1];
                        d__[j - 1] = dmnmx;
                    } else {
                        goto L30;
                    }
                }
            L30:;
            }
        } else {
            i__1 = endd;
            for (i__ = start + 1; i__ <= i__1; ++i__) {
                i__2 = start + 1;
                for (j = i__; j >= i__2; --j) {
                    if (d__[j] < d__[j - 1]) {
                        dmnmx = d__[j];
                        d__[j] = d__[j - 1];
                        d__[j - 1] = dmnmx;
                    } else {
                        goto L50;
                    }
                }
            L50:;
            }
        }
    } else if (endd - start > 20) {
        d1 = d__[start];
        d2 = d__[endd];
        i__ = (start + endd) / 2;
        d3 = d__[i__];
        if (d1 < d2) {
            if (d3 < d1) {
                dmnmx = d1;
            } else if (d3 < d2) {
                dmnmx = d3;
            } else {
                dmnmx = d2;
            }
        } else {
            if (d3 < d2) {
                dmnmx = d2;
            } else if (d3 < d1) {
                dmnmx = d3;
            } else {
                dmnmx = d1;
            }
        }
        if (dir == 0) {
            i__ = start - 1;
            j = endd + 1;
        L60:
        L70:
            --j;
            if (d__[j] < dmnmx) {
                goto L70;
            }
        L80:
            ++i__;
            if (d__[i__] > dmnmx) {
                goto L80;
            }
            if (i__ < j) {
                tmp = d__[i__];
                d__[i__] = d__[j];
                d__[j] = tmp;
                goto L60;
            }
            if (j - start > endd - j - 1) {
                ++stkpnt;
                stack[(stkpnt << 1) - 2] = start;
                stack[(stkpnt << 1) - 1] = j;
                ++stkpnt;
                stack[(stkpnt << 1) - 2] = j + 1;
                stack[(stkpnt << 1) - 1] = endd;
            } else {
                ++stkpnt;
                stack[(stkpnt << 1) - 2] = j + 1;
                stack[(stkpnt << 1) - 1] = endd;
                ++stkpnt;
                stack[(stkpnt << 1) - 2] = start;
                stack[(stkpnt << 1) - 1] = j;
            }
        } else {
            i__ = start - 1;
            j = endd + 1;
        L90:
        L100:
            --j;
            if (d__[j] > dmnmx) {
                goto L100;
            }
        L110:
            ++i__;
            if (d__[i__] < dmnmx) {
                goto L110;
            }
            if (i__ < j) {
                tmp = d__[i__];
                d__[i__] = d__[j];
                d__[j] = tmp;
                goto L90;
            }
            if (j - start > endd - j - 1) {
                ++stkpnt;
                stack[(stkpnt << 1) - 2] = start;
                stack[(stkpnt << 1) - 1] = j;
                ++stkpnt;
                stack[(stkpnt << 1) - 2] = j + 1;
                stack[(stkpnt << 1) - 1] = endd;
            } else {
                ++stkpnt;
                stack[(stkpnt << 1) - 2] = j + 1;
                stack[(stkpnt << 1) - 1] = endd;
                ++stkpnt;
                stack[(stkpnt << 1) - 2] = start;
                stack[(stkpnt << 1) - 1] = j;
            }
        }
    }
    if (stkpnt > 0) {
        goto L10;
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
