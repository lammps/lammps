/* fortran/ztptri.f -- translated by f2c (version 20200916).
   You must link the resulting object file with libf2c:
        on Microsoft Windows system, link with libf2c.lib;
        on Linux or Unix systems, link with .../path/to/libf2c.a -lm
        or, if you install libf2c.a in a standard place, with -lf2c -lm
        -- in that order, at the end of the command line, as in
                cc *.o -lf2c -lm
        Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

                http://www.netlib.org/f2c/libf2c.zip
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"

/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZTPTRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZTPTRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztptri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztptri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztptri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTPTRI( UPLO, DIAG, N, AP, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, UPLO */
/*       INTEGER            INFO, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTPTRI computes the inverse of a complex upper or lower triangular */
/* > matrix A stored in packed format. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  A is upper triangular; */
/* >          = 'L':  A is lower triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >          = 'N':  A is non-unit triangular; */
/* >          = 'U':  A is unit triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangular matrix A, stored */
/* >          columnwise in a linear array.  The j-th column of A is stored */
/* >          in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*((2*n-j)/2) = A(i,j) for j<=i<=n. */
/* >          See below for further details. */
/* >          On exit, the (triangular) inverse of the original matrix, in */
/* >          the same packed storage format. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, A(i,i) is exactly zero.  The triangular */
/* >                matrix is singular and its inverse can not be computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complex16OTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  A triangular matrix A can be transferred to packed storage using one */
/* >  of the following program segments: */
/* > */
/* >  UPLO = 'U':                      UPLO = 'L': */
/* > */
/* >        JC = 1                           JC = 1 */
/* >        DO 2 J = 1, N                    DO 2 J = 1, N */
/* >           DO 1 I = 1, J                    DO 1 I = J, N */
/* >              AP(JC+I-1) = A(I,J)              AP(JC+I-J) = A(I,J) */
/* >      1    CONTINUE                    1    CONTINUE */
/* >           JC = JC + J                      JC = JC + N - J + 1 */
/* >      2 CONTINUE                       2 CONTINUE */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ztptri_(char *uplo, char *diag, integer *n,
        doublecomplex *ap, integer *info, ftnlen uplo_len, ftnlen diag_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublecomplex z__1;

    /* Builtin functions */
    void z_lmp_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    integer j, jc, jj;
    doublecomplex ajj;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *,
            doublecomplex *, integer *);
    logical upper;
    extern /* Subroutine */ int ztpmv_(char *, char *, char *, integer *,
            doublecomplex *, doublecomplex *, integer *, ftnlen, ftnlen,
            ftnlen), xerbla_(char *, integer *, ftnlen);
    integer jclast;
    logical nounit;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --ap;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    nounit = lsame_(diag, (char *)"N", (ftnlen)1, (ftnlen)1);
    if (! upper && ! lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (! nounit && ! lsame_(diag, (char *)"U", (ftnlen)1, (ftnlen)1)) {
        *info = -2;
    } else if (*n < 0) {
        *info = -3;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"ZTPTRI", &i__1, (ftnlen)6);
        return 0;
    }

/*     Check for singularity if non-unit. */

    if (nounit) {
        if (upper) {
            jj = 0;
            i__1 = *n;
            for (*info = 1; *info <= i__1; ++(*info)) {
                jj += *info;
                i__2 = jj;
                if (ap[i__2].r == 0. && ap[i__2].i == 0.) {
                    return 0;
                }
/* L10: */
            }
        } else {
            jj = 1;
            i__1 = *n;
            for (*info = 1; *info <= i__1; ++(*info)) {
                i__2 = jj;
                if (ap[i__2].r == 0. && ap[i__2].i == 0.) {
                    return 0;
                }
                jj = jj + *n - *info + 1;
/* L20: */
            }
        }
        *info = 0;
    }

    if (upper) {

/*        Compute inverse of upper triangular matrix. */

        jc = 1;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            if (nounit) {
                i__2 = jc + j - 1;
                z_lmp_div(&z__1, &c_b1, &ap[jc + j - 1]);
                ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
                i__2 = jc + j - 1;
                z__1.r = -ap[i__2].r, z__1.i = -ap[i__2].i;
                ajj.r = z__1.r, ajj.i = z__1.i;
            } else {
                z__1.r = -1., z__1.i = -0.;
                ajj.r = z__1.r, ajj.i = z__1.i;
            }

/*           Compute elements 1:j-1 of j-th column. */

            i__2 = j - 1;
            ztpmv_((char *)"Upper", (char *)"No transpose", diag, &i__2, &ap[1], &ap[jc], &
                    c__1, (ftnlen)5, (ftnlen)12, (ftnlen)1);
            i__2 = j - 1;
            zscal_(&i__2, &ajj, &ap[jc], &c__1);
            jc += j;
/* L30: */
        }

    } else {

/*        Compute inverse of lower triangular matrix. */

        jc = *n * (*n + 1) / 2;
        for (j = *n; j >= 1; --j) {
            if (nounit) {
                i__1 = jc;
                z_lmp_div(&z__1, &c_b1, &ap[jc]);
                ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
                i__1 = jc;
                z__1.r = -ap[i__1].r, z__1.i = -ap[i__1].i;
                ajj.r = z__1.r, ajj.i = z__1.i;
            } else {
                z__1.r = -1., z__1.i = -0.;
                ajj.r = z__1.r, ajj.i = z__1.i;
            }
            if (j < *n) {

/*              Compute elements j+1:n of j-th column. */

                i__1 = *n - j;
                ztpmv_((char *)"Lower", (char *)"No transpose", diag, &i__1, &ap[jclast], &ap[
                        jc + 1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)1);
                i__1 = *n - j;
                zscal_(&i__1, &ajj, &ap[jc + 1], &c__1);
            }
            jclast = jc;
            jc = jc - *n + j - 2;
/* L40: */
        }
    }

    return 0;

/*     End of ZTPTRI */

} /* ztptri_ */

#ifdef __cplusplus
        }
#endif
