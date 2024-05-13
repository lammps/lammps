#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlasdt_(integer *n, integer *lvl, integer *nd, integer *inode, integer *ndiml, integer *ndimr,
            integer *msub)
{
    integer i__1, i__2;
    double log(doublereal);
    integer i__, il, ir, maxn;
    doublereal temp;
    integer nlvl, llst, ncrnt;
    --ndimr;
    --ndiml;
    --inode;
    maxn = max(1, *n);
    temp = log((doublereal)maxn / (doublereal)(*msub + 1)) / log(2.);
    *lvl = (integer)temp + 1;
    i__ = *n / 2;
    inode[1] = i__ + 1;
    ndiml[1] = i__;
    ndimr[1] = *n - i__ - 1;
    il = 0;
    ir = 1;
    llst = 1;
    i__1 = *lvl - 1;
    for (nlvl = 1; nlvl <= i__1; ++nlvl) {
        i__2 = llst - 1;
        for (i__ = 0; i__ <= i__2; ++i__) {
            il += 2;
            ir += 2;
            ncrnt = llst + i__;
            ndiml[il] = ndiml[ncrnt] / 2;
            ndimr[il] = ndiml[ncrnt] - ndiml[il] - 1;
            inode[il] = inode[ncrnt] - ndimr[il] - 1;
            ndiml[ir] = ndimr[ncrnt] / 2;
            ndimr[ir] = ndimr[ncrnt] - ndiml[ir] - 1;
            inode[ir] = inode[ncrnt] + ndiml[ir] + 1;
        }
        llst <<= 1;
    }
    *nd = (llst << 1) - 1;
    return 0;
}
#ifdef __cplusplus
}
#endif
