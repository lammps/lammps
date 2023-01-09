#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlabad_(doublereal *small, doublereal *large)
{
    double d_lmp_lg10(doublereal *), sqrt(doublereal);
    if (d_lmp_lg10(large) > 2e3) {
        *small = sqrt(*small);
        *large = sqrt(*large);
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
