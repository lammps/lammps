#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlarnv_(integer *idist, integer *iseed, integer *n, doublereal *x)
{
    integer i__1, i__2, i__3;
    double log(doublereal), sqrt(doublereal), cos(doublereal);
    integer i__;
    doublereal u[128];
    integer il, iv, il2;
    extern int dlaruv_(integer *, integer *, doublereal *);
    --x;
    --iseed;
    i__1 = *n;
    for (iv = 1; iv <= i__1; iv += 64) {
        i__2 = 64, i__3 = *n - iv + 1;
        il = min(i__2, i__3);
        if (*idist == 3) {
            il2 = il << 1;
        } else {
            il2 = il;
        }
        dlaruv_(&iseed[1], &il2, u);
        if (*idist == 1) {
            i__2 = il;
            for (i__ = 1; i__ <= i__2; ++i__) {
                x[iv + i__ - 1] = u[i__ - 1];
            }
        } else if (*idist == 2) {
            i__2 = il;
            for (i__ = 1; i__ <= i__2; ++i__) {
                x[iv + i__ - 1] = u[i__ - 1] * 2. - 1.;
            }
        } else if (*idist == 3) {
            i__2 = il;
            for (i__ = 1; i__ <= i__2; ++i__) {
                x[iv + i__ - 1] = sqrt(log(u[(i__ << 1) - 2]) * -2.) *
                                  cos(u[(i__ << 1) - 1] * 6.28318530717958647692528676655900576839);
            }
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
