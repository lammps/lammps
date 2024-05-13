#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
doublereal zlanhe_(char *norm, char *uplo, integer *n, doublecomplex *a, integer *lda,
                   doublereal *work, ftnlen norm_len, ftnlen uplo_len)
{
    integer a_dim1, a_offset, i__1, i__2;
    doublereal ret_val, d__1;
    double z_lmp_abs(doublecomplex *), sqrt(doublereal);
    integer i__, j;
    doublereal sum, absa, scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    doublereal value;
    extern logical disnan_(doublereal *);
    extern int zlassq_(integer *, doublecomplex *, integer *, doublereal *, doublereal *);
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --work;
    if (*n == 0) {
        value = 0.;
    } else if (lsame_(norm, (char *)"M", (ftnlen)1, (ftnlen)1)) {
        value = 0.;
        if (lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1)) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = j - 1;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    sum = z_lmp_abs(&a[i__ + j * a_dim1]);
                    if (value < sum || disnan_(&sum)) {
                        value = sum;
                    }
                }
                i__2 = j + j * a_dim1;
                sum = (d__1 = a[i__2].r, abs(d__1));
                if (value < sum || disnan_(&sum)) {
                    value = sum;
                }
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = j + j * a_dim1;
                sum = (d__1 = a[i__2].r, abs(d__1));
                if (value < sum || disnan_(&sum)) {
                    value = sum;
                }
                i__2 = *n;
                for (i__ = j + 1; i__ <= i__2; ++i__) {
                    sum = z_lmp_abs(&a[i__ + j * a_dim1]);
                    if (value < sum || disnan_(&sum)) {
                        value = sum;
                    }
                }
            }
        }
    } else if (lsame_(norm, (char *)"I", (ftnlen)1, (ftnlen)1) || lsame_(norm, (char *)"O", (ftnlen)1, (ftnlen)1) ||
               *(unsigned char *)norm == '1') {
        value = 0.;
        if (lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1)) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                sum = 0.;
                i__2 = j - 1;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    absa = z_lmp_abs(&a[i__ + j * a_dim1]);
                    sum += absa;
                    work[i__] += absa;
                }
                i__2 = j + j * a_dim1;
                work[j] = sum + (d__1 = a[i__2].r, abs(d__1));
            }
            i__1 = *n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                sum = work[i__];
                if (value < sum || disnan_(&sum)) {
                    value = sum;
                }
            }
        } else {
            i__1 = *n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                work[i__] = 0.;
            }
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = j + j * a_dim1;
                sum = work[j] + (d__1 = a[i__2].r, abs(d__1));
                i__2 = *n;
                for (i__ = j + 1; i__ <= i__2; ++i__) {
                    absa = z_lmp_abs(&a[i__ + j * a_dim1]);
                    sum += absa;
                    work[i__] += absa;
                }
                if (value < sum || disnan_(&sum)) {
                    value = sum;
                }
            }
        }
    } else if (lsame_(norm, (char *)"F", (ftnlen)1, (ftnlen)1) || lsame_(norm, (char *)"E", (ftnlen)1, (ftnlen)1)) {
        scale = 0.;
        sum = 1.;
        if (lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1)) {
            i__1 = *n;
            for (j = 2; j <= i__1; ++j) {
                i__2 = j - 1;
                zlassq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
            }
        } else {
            i__1 = *n - 1;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *n - j;
                zlassq_(&i__2, &a[j + 1 + j * a_dim1], &c__1, &scale, &sum);
            }
        }
        sum *= 2;
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            i__2 = i__ + i__ * a_dim1;
            if (a[i__2].r != 0.) {
                i__2 = i__ + i__ * a_dim1;
                absa = (d__1 = a[i__2].r, abs(d__1));
                if (scale < absa) {
                    d__1 = scale / absa;
                    sum = sum * (d__1 * d__1) + 1.;
                    scale = absa;
                } else {
                    d__1 = absa / scale;
                    sum += d__1 * d__1;
                }
            }
        }
        value = scale * sqrt(sum);
    }
    ret_val = value;
    return ret_val;
}
#ifdef __cplusplus
}
#endif
