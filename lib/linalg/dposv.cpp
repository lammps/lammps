/* fortran/dposv.f -- translated by f2c (version 20200916).
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

/* > \brief <b> DPOSV computes the solution to system of linear equations A * X = B for PO matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPOSV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dposv.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dposv.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dposv.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO ) */

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
/* > DPOSV computes the solution to a real system of linear equations */
/* >    A * X = B, */
/* > where A is an N-by-N symmetric positive definite matrix and X and B */
/* > are N-by-NRHS matrices. */
/* > */
/* > The Cholesky decomposition is used to factor A as */
/* >    A = U**T* U,  if UPLO = 'U', or */
/* >    A = L * L**T,  if UPLO = 'L', */
/* > where U is an upper triangular matrix and L is a lower triangular */
/* > matrix.  The factored form of A is then used to solve the system of */
/* > equations A * X = B. */
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
/* >          The number of linear equations, i.e., the order of the */
/* >          matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrix B.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
/* >          N-by-N upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading N-by-N lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* > */
/* >          On exit, if INFO = 0, the factor U or L from the Cholesky */
/* >          factorization A = U**T*U or A = L*L**T. */
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
/* >          On entry, the N-by-NRHS right hand side matrix B. */
/* >          On exit, if INFO = 0, the N-by-NRHS solution matrix X. */
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
/* >          > 0:  if INFO = i, the leading minor of order i of A is not */
/* >                positive definite, so the factorization could not be */
/* >                completed, and the solution has not been computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doublePOsolve */

/*  ===================================================================== */
/* Subroutine */ int dposv_(char *uplo, integer *n, integer *nrhs, doublereal 
	*a, integer *lda, doublereal *b, integer *ldb, integer *info, ftnlen 
	uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dpotrf_(
	    char *, integer *, doublereal *, integer *, integer *, ftnlen), 
	    dpotrs_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen);


/*  -- LAPACK driver routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

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
    if (! lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, (char *)"L", (
	    ftnlen)1, (ftnlen)1)) {
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
	xerbla_((char *)"DPOSV ", &i__1, (ftnlen)6);
	return 0;
    }

/*     Compute the Cholesky factorization A = U**T*U or A = L*L**T. */

    dpotrf_(uplo, n, &a[a_offset], lda, info, (ftnlen)1);
    if (*info == 0) {

/*        Solve the system A*X = B, overwriting B with X. */

	dpotrs_(uplo, n, nrhs, &a[a_offset], lda, &b[b_offset], ldb, info, (
		ftnlen)1);

    }
    return 0;

/*     End of DPOSV */

} /* dposv_ */

#ifdef __cplusplus
	}
#endif
