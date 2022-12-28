/* fortran/dpotrs.f -- translated by f2c (version 20200916).
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

static doublereal c_b9 = 1.;

/* > \brief \b DPOTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPOTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpotrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpotrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpotrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDB, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPOTRS solves a system of linear equations A*X = B with a symmetric */
/* > positive definite matrix A using the Cholesky factorization */
/* > A = U**T*U or A = L*L**T computed by DPOTRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored; */
/* >          = 'L':  Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrix B.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          The triangular factor U or L from the Cholesky factorization */
/* >          A = U**T*U or A = L*L**T, as computed by DPOTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
/* >          On entry, the right hand side matrix B. */
/* >          On exit, the solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doublePOcomputational */

/*  ===================================================================== */
/* Subroutine */ int dpotrs_(char *uplo, integer *n, integer *nrhs,
        doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
        info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *,
            integer *, integer *, doublereal *, doublereal *, integer *,
            doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    if (! upper && ! lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*nrhs < 0) {
        *info = -3;
    } else if (*lda < max(1,*n)) {
        *info = -5;
    } else if (*ldb < max(1,*n)) {
        *info = -7;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DPOTRS", &i__1, (ftnlen)6);
        return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *nrhs == 0) {
        return 0;
    }

    if (upper) {

/*        Solve A*X = B where A = U**T *U. */

/*        Solve U**T *X = B, overwriting B with X. */

        dtrsm_((char *)"Left", (char *)"Upper", (char *)"Transpose", (char *)"Non-unit", n, nrhs, &c_b9, &a[
                a_offset], lda, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (
                ftnlen)9, (ftnlen)8);

/*        Solve U*X = B, overwriting B with X. */

        dtrsm_((char *)"Left", (char *)"Upper", (char *)"No transpose", (char *)"Non-unit", n, nrhs, &c_b9, &
                a[a_offset], lda, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (
                ftnlen)12, (ftnlen)8);
    } else {

/*        Solve A*X = B where A = L*L**T. */

/*        Solve L*X = B, overwriting B with X. */

        dtrsm_((char *)"Left", (char *)"Lower", (char *)"No transpose", (char *)"Non-unit", n, nrhs, &c_b9, &
                a[a_offset], lda, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (
                ftnlen)12, (ftnlen)8);

/*        Solve L**T *X = B, overwriting B with X. */

        dtrsm_((char *)"Left", (char *)"Lower", (char *)"Transpose", (char *)"Non-unit", n, nrhs, &c_b9, &a[
                a_offset], lda, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (
                ftnlen)9, (ftnlen)8);
    }

    return 0;

/*     End of DPOTRS */

} /* dpotrs_ */

#ifdef __cplusplus
        }
#endif
