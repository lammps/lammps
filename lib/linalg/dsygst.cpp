/* fortran/dsygst.f -- translated by f2c (version 20200916).
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

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b14 = 1.;
static doublereal c_b16 = -.5;
static doublereal c_b19 = -1.;
static doublereal c_b52 = .5;

/* > \brief \b DSYGST */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYGST + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsygst.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsygst.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsygst.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, ITYPE, LDA, LDB, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYGST reduces a real symmetric-definite generalized eigenproblem */
/* > to standard form. */
/* > */
/* > If ITYPE = 1, the problem is A*x = lambda*B*x, */
/* > and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T) */
/* > */
/* > If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or */
/* > B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L. */
/* > */
/* > B must have been previously factorized as U**T*U or L*L**T by DPOTRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ITYPE */
/* > \verbatim */
/* >          ITYPE is INTEGER */
/* >          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T); */
/* >          = 2 or 3: compute U*A*U**T or L**T*A*L. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored and B is factored as */
/* >                  U**T*U; */
/* >          = 'L':  Lower triangle of A is stored and B is factored as */
/* >                  L*L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices A and B.  N >= 0. */
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
/* >          On exit, if INFO = 0, the transformed matrix, stored in the */
/* >          same format as A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,N) */
/* >          The triangular factor from the Cholesky factorization of B, */
/* >          as returned by DPOTRF. */
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

/* > \ingroup doubleSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int dsygst_(integer *itype, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    integer k, kb, nb;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dsymm_(
	    char *, char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen, ftnlen);
    logical upper;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dsygs2_(
	    integer *, char *, integer *, doublereal *, integer *, doublereal 
	    *, integer *, integer *, ftnlen), dsyr2k_(char *, char *, integer 
	    *, integer *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen)
	    , xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
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
    if (*itype < 1 || *itype > 3) {
	*info = -1;
    } else if (! upper && ! lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*ldb < max(1,*n)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_((char *)"DSYGST", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Determine the block size for this environment. */

    nb = ilaenv_(&c__1, (char *)"DSYGST", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);

    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code */

	dsygs2_(itype, uplo, n, &a[a_offset], lda, &b[b_offset], ldb, info, (
		ftnlen)1);
    } else {

/*        Use blocked code */

	if (*itype == 1) {
	    if (upper) {

/*              Compute inv(U**T)*A*inv(U) */

		i__1 = *n;
		i__2 = nb;
		for (k = 1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2) {
/* Computing MIN */
		    i__3 = *n - k + 1;
		    kb = min(i__3,nb);

/*                 Update the upper triangle of A(k:n,k:n) */

		    dsygs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + 
			    k * b_dim1], ldb, info, (ftnlen)1);
		    if (k + kb <= *n) {
			i__3 = *n - k - kb + 1;
			dtrsm_((char *)"Left", uplo, (char *)"Transpose", (char *)"Non-unit", &kb, &
				i__3, &c_b14, &b[k + k * b_dim1], ldb, &a[k + 
				(k + kb) * a_dim1], lda, (ftnlen)4, (ftnlen)1,
				 (ftnlen)9, (ftnlen)8);
			i__3 = *n - k - kb + 1;
			dsymm_((char *)"Left", uplo, &kb, &i__3, &c_b16, &a[k + k * 
				a_dim1], lda, &b[k + (k + kb) * b_dim1], ldb, 
				&c_b14, &a[k + (k + kb) * a_dim1], lda, (
				ftnlen)4, (ftnlen)1);
			i__3 = *n - k - kb + 1;
			dsyr2k_(uplo, (char *)"Transpose", &i__3, &kb, &c_b19, &a[k + 
				(k + kb) * a_dim1], lda, &b[k + (k + kb) * 
				b_dim1], ldb, &c_b14, &a[k + kb + (k + kb) * 
				a_dim1], lda, (ftnlen)1, (ftnlen)9);
			i__3 = *n - k - kb + 1;
			dsymm_((char *)"Left", uplo, &kb, &i__3, &c_b16, &a[k + k * 
				a_dim1], lda, &b[k + (k + kb) * b_dim1], ldb, 
				&c_b14, &a[k + (k + kb) * a_dim1], lda, (
				ftnlen)4, (ftnlen)1);
			i__3 = *n - k - kb + 1;
			dtrsm_((char *)"Right", uplo, (char *)"No transpose", (char *)"Non-unit", &kb,
				 &i__3, &c_b14, &b[k + kb + (k + kb) * b_dim1]
				, ldb, &a[k + (k + kb) * a_dim1], lda, (
				ftnlen)5, (ftnlen)1, (ftnlen)12, (ftnlen)8);
		    }
/* L10: */
		}
	    } else {

/*              Compute inv(L)*A*inv(L**T) */

		i__2 = *n;
		i__1 = nb;
		for (k = 1; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1) {
/* Computing MIN */
		    i__3 = *n - k + 1;
		    kb = min(i__3,nb);

/*                 Update the lower triangle of A(k:n,k:n) */

		    dsygs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + 
			    k * b_dim1], ldb, info, (ftnlen)1);
		    if (k + kb <= *n) {
			i__3 = *n - k - kb + 1;
			dtrsm_((char *)"Right", uplo, (char *)"Transpose", (char *)"Non-unit", &i__3, 
				&kb, &c_b14, &b[k + k * b_dim1], ldb, &a[k + 
				kb + k * a_dim1], lda, (ftnlen)5, (ftnlen)1, (
				ftnlen)9, (ftnlen)8);
			i__3 = *n - k - kb + 1;
			dsymm_((char *)"Right", uplo, &i__3, &kb, &c_b16, &a[k + k * 
				a_dim1], lda, &b[k + kb + k * b_dim1], ldb, &
				c_b14, &a[k + kb + k * a_dim1], lda, (ftnlen)
				5, (ftnlen)1);
			i__3 = *n - k - kb + 1;
			dsyr2k_(uplo, (char *)"No transpose", &i__3, &kb, &c_b19, &a[
				k + kb + k * a_dim1], lda, &b[k + kb + k * 
				b_dim1], ldb, &c_b14, &a[k + kb + (k + kb) * 
				a_dim1], lda, (ftnlen)1, (ftnlen)12);
			i__3 = *n - k - kb + 1;
			dsymm_((char *)"Right", uplo, &i__3, &kb, &c_b16, &a[k + k * 
				a_dim1], lda, &b[k + kb + k * b_dim1], ldb, &
				c_b14, &a[k + kb + k * a_dim1], lda, (ftnlen)
				5, (ftnlen)1);
			i__3 = *n - k - kb + 1;
			dtrsm_((char *)"Left", uplo, (char *)"No transpose", (char *)"Non-unit", &
				i__3, &kb, &c_b14, &b[k + kb + (k + kb) * 
				b_dim1], ldb, &a[k + kb + k * a_dim1], lda, (
				ftnlen)4, (ftnlen)1, (ftnlen)12, (ftnlen)8);
		    }
/* L20: */
		}
	    }
	} else {
	    if (upper) {

/*              Compute U*A*U**T */

		i__1 = *n;
		i__2 = nb;
		for (k = 1; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2) {
/* Computing MIN */
		    i__3 = *n - k + 1;
		    kb = min(i__3,nb);

/*                 Update the upper triangle of A(1:k+kb-1,1:k+kb-1) */

		    i__3 = k - 1;
		    dtrmm_((char *)"Left", uplo, (char *)"No transpose", (char *)"Non-unit", &i__3, &
			    kb, &c_b14, &b[b_offset], ldb, &a[k * a_dim1 + 1],
			     lda, (ftnlen)4, (ftnlen)1, (ftnlen)12, (ftnlen)8)
			    ;
		    i__3 = k - 1;
		    dsymm_((char *)"Right", uplo, &i__3, &kb, &c_b52, &a[k + k * 
			    a_dim1], lda, &b[k * b_dim1 + 1], ldb, &c_b14, &a[
			    k * a_dim1 + 1], lda, (ftnlen)5, (ftnlen)1);
		    i__3 = k - 1;
		    dsyr2k_(uplo, (char *)"No transpose", &i__3, &kb, &c_b14, &a[k * 
			    a_dim1 + 1], lda, &b[k * b_dim1 + 1], ldb, &c_b14,
			     &a[a_offset], lda, (ftnlen)1, (ftnlen)12);
		    i__3 = k - 1;
		    dsymm_((char *)"Right", uplo, &i__3, &kb, &c_b52, &a[k + k * 
			    a_dim1], lda, &b[k * b_dim1 + 1], ldb, &c_b14, &a[
			    k * a_dim1 + 1], lda, (ftnlen)5, (ftnlen)1);
		    i__3 = k - 1;
		    dtrmm_((char *)"Right", uplo, (char *)"Transpose", (char *)"Non-unit", &i__3, &kb,
			     &c_b14, &b[k + k * b_dim1], ldb, &a[k * a_dim1 + 
			    1], lda, (ftnlen)5, (ftnlen)1, (ftnlen)9, (ftnlen)
			    8);
		    dsygs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + 
			    k * b_dim1], ldb, info, (ftnlen)1);
/* L30: */
		}
	    } else {

/*              Compute L**T*A*L */

		i__2 = *n;
		i__1 = nb;
		for (k = 1; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1) {
/* Computing MIN */
		    i__3 = *n - k + 1;
		    kb = min(i__3,nb);

/*                 Update the lower triangle of A(1:k+kb-1,1:k+kb-1) */

		    i__3 = k - 1;
		    dtrmm_((char *)"Right", uplo, (char *)"No transpose", (char *)"Non-unit", &kb, &
			    i__3, &c_b14, &b[b_offset], ldb, &a[k + a_dim1], 
			    lda, (ftnlen)5, (ftnlen)1, (ftnlen)12, (ftnlen)8);
		    i__3 = k - 1;
		    dsymm_((char *)"Left", uplo, &kb, &i__3, &c_b52, &a[k + k * 
			    a_dim1], lda, &b[k + b_dim1], ldb, &c_b14, &a[k + 
			    a_dim1], lda, (ftnlen)4, (ftnlen)1);
		    i__3 = k - 1;
		    dsyr2k_(uplo, (char *)"Transpose", &i__3, &kb, &c_b14, &a[k + 
			    a_dim1], lda, &b[k + b_dim1], ldb, &c_b14, &a[
			    a_offset], lda, (ftnlen)1, (ftnlen)9);
		    i__3 = k - 1;
		    dsymm_((char *)"Left", uplo, &kb, &i__3, &c_b52, &a[k + k * 
			    a_dim1], lda, &b[k + b_dim1], ldb, &c_b14, &a[k + 
			    a_dim1], lda, (ftnlen)4, (ftnlen)1);
		    i__3 = k - 1;
		    dtrmm_((char *)"Left", uplo, (char *)"Transpose", (char *)"Non-unit", &kb, &i__3, 
			    &c_b14, &b[k + k * b_dim1], ldb, &a[k + a_dim1], 
			    lda, (ftnlen)4, (ftnlen)1, (ftnlen)9, (ftnlen)8);
		    dsygs2_(itype, uplo, &kb, &a[k + k * a_dim1], lda, &b[k + 
			    k * b_dim1], ldb, info, (ftnlen)1);
/* L40: */
		}
	    }
	}
    }
    return 0;

/*     End of DSYGST */

} /* dsygst_ */

#ifdef __cplusplus
	}
#endif
