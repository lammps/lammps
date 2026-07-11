#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dcombssq_(doublereal *v1, doublereal *v2)
{
    doublereal d__1;
    --v2;
    --v1;
    if (v1[1] >= v2[1]) {
        if (v1[1] != 0.) {
            d__1 = v2[1] / v1[1];
            v1[2] += d__1 * d__1 * v2[2];
        } else {
            v1[2] += v2[2];
        }
    } else {
        d__1 = v1[1] / v2[1];
        v1[2] = v2[2] + d__1 * d__1 * v1[2];
        v1[1] = v2[1];
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
