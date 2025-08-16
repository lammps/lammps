#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
int dlaqr1_(integer *n, doublereal *h__, integer *ldh, doublereal *sr1, doublereal *si1,
            doublereal *sr2, doublereal *si2, doublereal *v)
{
    integer h_dim1, h_offset;
    doublereal d__1, d__2, d__3;
    doublereal s, h21s, h31s;
    h_dim1 = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --v;
    if (*n != 2 && *n != 3) {
        return 0;
    }
    if (*n == 2) {
        s = (d__1 = h__[h_dim1 + 1] - *sr2, abs(d__1)) + abs(*si2) +
            (d__2 = h__[h_dim1 + 2], abs(d__2));
        if (s == 0.) {
            v[1] = 0.;
            v[2] = 0.;
        } else {
            h21s = h__[h_dim1 + 2] / s;
            v[1] = h21s * h__[(h_dim1 << 1) + 1] +
                   (h__[h_dim1 + 1] - *sr1) * ((h__[h_dim1 + 1] - *sr2) / s) - *si1 * (*si2 / s);
            v[2] = h21s * (h__[h_dim1 + 1] + h__[(h_dim1 << 1) + 2] - *sr1 - *sr2);
        }
    } else {
        s = (d__1 = h__[h_dim1 + 1] - *sr2, abs(d__1)) + abs(*si2) +
            (d__2 = h__[h_dim1 + 2], abs(d__2)) + (d__3 = h__[h_dim1 + 3], abs(d__3));
        if (s == 0.) {
            v[1] = 0.;
            v[2] = 0.;
            v[3] = 0.;
        } else {
            h21s = h__[h_dim1 + 2] / s;
            h31s = h__[h_dim1 + 3] / s;
            v[1] = (h__[h_dim1 + 1] - *sr1) * ((h__[h_dim1 + 1] - *sr2) / s) - *si1 * (*si2 / s) +
                   h__[(h_dim1 << 1) + 1] * h21s + h__[h_dim1 * 3 + 1] * h31s;
            v[2] = h21s * (h__[h_dim1 + 1] + h__[(h_dim1 << 1) + 2] - *sr1 - *sr2) +
                   h__[h_dim1 * 3 + 2] * h31s;
            v[3] = h31s * (h__[h_dim1 + 1] + h__[h_dim1 * 3 + 3] - *sr1 - *sr2) +
                   h21s * h__[(h_dim1 << 1) + 3];
        }
    }
    return 0;
}
#ifdef __cplusplus
}
#endif
