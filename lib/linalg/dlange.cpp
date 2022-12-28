#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
doublereal dlange_(char *norm, integer *m, integer *n, doublereal *a, integer *lda,
                   doublereal *work, ftnlen norm_len)
{
    integer a_dim1, a_offset, i__1, i__2;
    doublereal ret_val, d__1;
    double sqrt(doublereal);
    integer i__, j;
    doublereal sum, temp, scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    doublereal value;
    extern logical disnan_(doublereal *);
    extern int dlassq_(integer *, doublereal *, integer *, doublereal *, doublereal *);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --work;
    if (min(*m, *n) == 0) {
        value = 0.;
    } else if (lsame_(norm, (char *)"M", (ftnlen)1, (ftnlen)1)) {
        value = 0.;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__2 = *m;
            for (i__ = 1; i__ <= i__2; ++i__) {
                temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
                if (value < temp || disnan_(&temp)) {
                    value = temp;
                }
            }
        }
    } else if (lsame_(norm, (char *)"O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {
        value = 0.;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            sum = 0.;
            i__2 = *m;
            for (i__ = 1; i__ <= i__2; ++i__) {
                sum += (d__1 = a[i__ + j * a_dim1], abs(d__1));
            }
            if (value < sum || disnan_(&sum)) {
                value = sum;
            }
        }
    } else if (lsame_(norm, (char *)"I", (ftnlen)1, (ftnlen)1)) {
        i__1 = *m;
        for (i__ = 1; i__ <= i__1; ++i__) {
            work[i__] = 0.;
        }
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__2 = *m;
            for (i__ = 1; i__ <= i__2; ++i__) {
                work[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1));
            }
        }
        value = 0.;
        i__1 = *m;
        for (i__ = 1; i__ <= i__1; ++i__) {
            temp = work[i__];
            if (value < temp || disnan_(&temp)) {
                value = temp;
            }
        }
    } else if (lsame_(norm, (char *)"F", (ftnlen)1, (ftnlen)1) || lsame_(norm, (char *)"E", (ftnlen)1, (ftnlen)1)) {
        scale = 0.;
        sum = 1.;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            dlassq_(m, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
        }
        value = scale * sqrt(sum);
    }
    ret_val = value;
    return ret_val;
}
#ifdef __cplusplus
}
#endif
