#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
static integer c__1 = 1;
int dgecon_(char *norm, integer *n, doublereal *a, integer *lda, doublereal *anorm,
            doublereal *rcond, doublereal *work, integer *iwork, integer *info, ftnlen norm_len)
{
    integer a_dim1, a_offset, i__1;
    doublereal d__1;
    doublereal sl;
    integer ix;
    doublereal su;
    integer kase, kase1;
    doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer isave[3];
    extern int drscl_(integer *, doublereal *, doublereal *, integer *),
        dlacn2_(integer *, doublereal *, doublereal *, integer *, doublereal *, integer *,
                integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern int xerbla_(char *, integer *, ftnlen);
    doublereal ainvnm;
    extern int dlatrs_(char *, char *, char *, char *, integer *, doublereal *, integer *,
                       doublereal *, doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen,
                       ftnlen);
    logical onenrm;
    char normin[1];
    doublereal smlnum;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --work;
    --iwork;
    *info = 0;
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, (char *)"O", (ftnlen)1, (ftnlen)1);
    if (!onenrm && !lsame_(norm, (char *)"I", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*lda < max(1, *n)) {
        *info = -4;
    } else if (*anorm < 0.) {
        *info = -5;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DGECON", &i__1, (ftnlen)6);
        return 0;
    }
    *rcond = 0.;
    if (*n == 0) {
        *rcond = 1.;
        return 0;
    } else if (*anorm == 0.) {
        return 0;
    }
    smlnum = dlamch_((char *)"Safe minimum", (ftnlen)12);
    ainvnm = 0.;
    *(unsigned char *)normin = 'N';
    if (onenrm) {
        kase1 = 1;
    } else {
        kase1 = 2;
    }
    kase = 0;
L10:
    dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
    if (kase != 0) {
        if (kase == kase1) {
            dlatrs_((char *)"Lower", (char *)"No transpose", (char *)"Unit", normin, n, &a[a_offset], lda, &work[1], &sl,
                    &work[(*n << 1) + 1], info, (ftnlen)5, (ftnlen)12, (ftnlen)4, (ftnlen)1);
            dlatrs_((char *)"Upper", (char *)"No transpose", (char *)"Non-unit", normin, n, &a[a_offset], lda, &work[1],
                    &su, &work[*n * 3 + 1], info, (ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
        } else {
            dlatrs_((char *)"Upper", (char *)"Transpose", (char *)"Non-unit", normin, n, &a[a_offset], lda, &work[1], &su,
                    &work[*n * 3 + 1], info, (ftnlen)5, (ftnlen)9, (ftnlen)8, (ftnlen)1);
            dlatrs_((char *)"Lower", (char *)"Transpose", (char *)"Unit", normin, n, &a[a_offset], lda, &work[1], &sl,
                    &work[(*n << 1) + 1], info, (ftnlen)5, (ftnlen)9, (ftnlen)4, (ftnlen)1);
        }
        scale = sl * su;
        *(unsigned char *)normin = 'Y';
        if (scale != 1.) {
            ix = idamax_(n, &work[1], &c__1);
            if (scale < (d__1 = work[ix], abs(d__1)) * smlnum || scale == 0.) {
                goto L20;
            }
            drscl_(n, &scale, &work[1], &c__1);
        }
        goto L10;
    }
    if (ainvnm != 0.) {
        *rcond = 1. / ainvnm / *anorm;
    }
L20:
    return 0;
}
#ifdef __cplusplus
}
#endif
