#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
doublereal dlanst_(char *norm, integer *n, doublereal *d__, doublereal *e, ftnlen norm_len)
{
    integer i__1;
    doublereal ret_val, d__1, d__2, d__3;
    double sqrt(doublereal);
    integer i__;
    doublereal sum, scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    doublereal anorm;
    extern logical disnan_(doublereal *);
    extern int dlassq_(integer *, doublereal *, integer *, doublereal *, doublereal *);
    --e;
    --d__;
    if (*n <= 0) {
        anorm = 0.;
    } else if (lsame_(norm, (char *)"M", (ftnlen)1, (ftnlen)1)) {
        anorm = (d__1 = d__[*n], abs(d__1));
        i__1 = *n - 1;
        for (i__ = 1; i__ <= i__1; ++i__) {
            sum = (d__1 = d__[i__], abs(d__1));
            if (anorm < sum || disnan_(&sum)) {
                anorm = sum;
            }
            sum = (d__1 = e[i__], abs(d__1));
            if (anorm < sum || disnan_(&sum)) {
                anorm = sum;
            }
        }
    } else if (lsame_(norm, (char *)"O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1' ||
               lsame_(norm, (char *)"I", (ftnlen)1, (ftnlen)1)) {
        if (*n == 1) {
            anorm = abs(d__[1]);
        } else {
            anorm = abs(d__[1]) + abs(e[1]);
            sum = (d__1 = e[*n - 1], abs(d__1)) + (d__2 = d__[*n], abs(d__2));
            if (anorm < sum || disnan_(&sum)) {
                anorm = sum;
            }
            i__1 = *n - 1;
            for (i__ = 2; i__ <= i__1; ++i__) {
                sum = (d__1 = d__[i__], abs(d__1)) + (d__2 = e[i__], abs(d__2)) +
                      (d__3 = e[i__ - 1], abs(d__3));
                if (anorm < sum || disnan_(&sum)) {
                    anorm = sum;
                }
            }
        }
    } else if (lsame_(norm, (char *)"F", (ftnlen)1, (ftnlen)1) || lsame_(norm, (char *)"E", (ftnlen)1, (ftnlen)1)) {
        scale = 0.;
        sum = 1.;
        if (*n > 1) {
            i__1 = *n - 1;
            dlassq_(&i__1, &e[1], &c__1, &scale, &sum);
            sum *= 2;
        }
        dlassq_(n, &d__[1], &c__1, &scale, &sum);
        anorm = scale * sqrt(sum);
    }
    ret_val = anorm;
    return ret_val;
}
#ifdef __cplusplus
}
#endif
