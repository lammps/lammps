/* fortran/dsygs2.f -- translated by f2c (version 20200916).
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

static doublereal c_b6 = -1.;
static integer c__1 = 1;
static doublereal c_b27 = 1.;

/* > \brief \b DSYGS2 reduces a symmetric definite generalized eigenproblem to standard form, using the factor
ization results obtained from spotrf (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYGS2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsygs2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsygs2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsygs2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO ) */

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
/* > DSYGS2 reduces a real symmetric-definite generalized eigenproblem */
/* > to standard form. */
/* > */
/* > If ITYPE = 1, the problem is A*x = lambda*B*x, */
/* > and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T) */
/* > */
/* > If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or */
/* > B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T *A*L. */
/* > */
/* > B must have been previously factorized as U**T *U or L*L**T by DPOTRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ITYPE */
/* > \verbatim */
/* >          ITYPE is INTEGER */
/* >          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T); */
/* >          = 2 or 3: compute U*A*U**T or L**T *A*L. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the upper or lower triangular part of the */
/* >          symmetric matrix A is stored, and how B has been factorized. */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
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
/* >          n by n upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading n by n lower triangular part of A contains the lower */
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
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doubleSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int dsygs2_(integer *itype, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    integer k;
    doublereal ct, akk, bkk;
    extern /* Subroutine */ int dsyr2_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen), dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    logical upper;
    extern /* Subroutine */ int dtrmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen), dtrsv_(char *, char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);


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
	xerbla_((char *)"DSYGS2", &i__1, (ftnlen)6);
	return 0;
    }

    if (*itype == 1) {
	if (upper) {

/*           Compute inv(U**T)*A*inv(U) */

	    i__1 = *n;
	    for (k = 1; k <= i__1; ++k) {

/*              Update the upper triangle of A(k:n,k:n) */

		akk = a[k + k * a_dim1];
		bkk = b[k + k * b_dim1];
/* Computing 2nd power */
		d__1 = bkk;
		akk /= d__1 * d__1;
		a[k + k * a_dim1] = akk;
		if (k < *n) {
		    i__2 = *n - k;
		    d__1 = 1. / bkk;
		    dscal_(&i__2, &d__1, &a[k + (k + 1) * a_dim1], lda);
		    ct = akk * -.5;
		    i__2 = *n - k;
		    daxpy_(&i__2, &ct, &b[k + (k + 1) * b_dim1], ldb, &a[k + (
			    k + 1) * a_dim1], lda);
		    i__2 = *n - k;
		    dsyr2_(uplo, &i__2, &c_b6, &a[k + (k + 1) * a_dim1], lda, 
			    &b[k + (k + 1) * b_dim1], ldb, &a[k + 1 + (k + 1) 
			    * a_dim1], lda, (ftnlen)1);
		    i__2 = *n - k;
		    daxpy_(&i__2, &ct, &b[k + (k + 1) * b_dim1], ldb, &a[k + (
			    k + 1) * a_dim1], lda);
		    i__2 = *n - k;
		    dtrsv_(uplo, (char *)"Transpose", (char *)"Non-unit", &i__2, &b[k + 1 + (
			    k + 1) * b_dim1], ldb, &a[k + (k + 1) * a_dim1], 
			    lda, (ftnlen)1, (ftnlen)9, (ftnlen)8);
		}
/* L10: */
	    }
	} else {

/*           Compute inv(L)*A*inv(L**T) */

	    i__1 = *n;
	    for (k = 1; k <= i__1; ++k) {

/*              Update the lower triangle of A(k:n,k:n) */

		akk = a[k + k * a_dim1];
		bkk = b[k + k * b_dim1];
/* Computing 2nd power */
		d__1 = bkk;
		akk /= d__1 * d__1;
		a[k + k * a_dim1] = akk;
		if (k < *n) {
		    i__2 = *n - k;
		    d__1 = 1. / bkk;
		    dscal_(&i__2, &d__1, &a[k + 1 + k * a_dim1], &c__1);
		    ct = akk * -.5;
		    i__2 = *n - k;
		    daxpy_(&i__2, &ct, &b[k + 1 + k * b_dim1], &c__1, &a[k + 
			    1 + k * a_dim1], &c__1);
		    i__2 = *n - k;
		    dsyr2_(uplo, &i__2, &c_b6, &a[k + 1 + k * a_dim1], &c__1, 
			    &b[k + 1 + k * b_dim1], &c__1, &a[k + 1 + (k + 1) 
			    * a_dim1], lda, (ftnlen)1);
		    i__2 = *n - k;
		    daxpy_(&i__2, &ct, &b[k + 1 + k * b_dim1], &c__1, &a[k + 
			    1 + k * a_dim1], &c__1);
		    i__2 = *n - k;
		    dtrsv_(uplo, (char *)"No transpose", (char *)"Non-unit", &i__2, &b[k + 1 
			    + (k + 1) * b_dim1], ldb, &a[k + 1 + k * a_dim1], 
			    &c__1, (ftnlen)1, (ftnlen)12, (ftnlen)8);
		}
/* L20: */
	    }
	}
    } else {
	if (upper) {

/*           Compute U*A*U**T */

	    i__1 = *n;
	    for (k = 1; k <= i__1; ++k) {

/*              Update the upper triangle of A(1:k,1:k) */

		akk = a[k + k * a_dim1];
		bkk = b[k + k * b_dim1];
		i__2 = k - 1;
		dtrmv_(uplo, (char *)"No transpose", (char *)"Non-unit", &i__2, &b[b_offset], 
			ldb, &a[k * a_dim1 + 1], &c__1, (ftnlen)1, (ftnlen)12,
			 (ftnlen)8);
		ct = akk * .5;
		i__2 = k - 1;
		daxpy_(&i__2, &ct, &b[k * b_dim1 + 1], &c__1, &a[k * a_dim1 + 
			1], &c__1);
		i__2 = k - 1;
		dsyr2_(uplo, &i__2, &c_b27, &a[k * a_dim1 + 1], &c__1, &b[k * 
			b_dim1 + 1], &c__1, &a[a_offset], lda, (ftnlen)1);
		i__2 = k - 1;
		daxpy_(&i__2, &ct, &b[k * b_dim1 + 1], &c__1, &a[k * a_dim1 + 
			1], &c__1);
		i__2 = k - 1;
		dscal_(&i__2, &bkk, &a[k * a_dim1 + 1], &c__1);
/* Computing 2nd power */
		d__1 = bkk;
		a[k + k * a_dim1] = akk * (d__1 * d__1);
/* L30: */
	    }
	} else {

/*           Compute L**T *A*L */

	    i__1 = *n;
	    for (k = 1; k <= i__1; ++k) {

/*              Update the lower triangle of A(1:k,1:k) */

		akk = a[k + k * a_dim1];
		bkk = b[k + k * b_dim1];
		i__2 = k - 1;
		dtrmv_(uplo, (char *)"Transpose", (char *)"Non-unit", &i__2, &b[b_offset], 
			ldb, &a[k + a_dim1], lda, (ftnlen)1, (ftnlen)9, (
			ftnlen)8);
		ct = akk * .5;
		i__2 = k - 1;
		daxpy_(&i__2, &ct, &b[k + b_dim1], ldb, &a[k + a_dim1], lda);
		i__2 = k - 1;
		dsyr2_(uplo, &i__2, &c_b27, &a[k + a_dim1], lda, &b[k + 
			b_dim1], ldb, &a[a_offset], lda, (ftnlen)1);
		i__2 = k - 1;
		daxpy_(&i__2, &ct, &b[k + b_dim1], ldb, &a[k + a_dim1], lda);
		i__2 = k - 1;
		dscal_(&i__2, &bkk, &a[k + a_dim1], lda);
/* Computing 2nd power */
		d__1 = bkk;
		a[k + k * a_dim1] = akk * (d__1 * d__1);
/* L40: */
	    }
	}
    }
    return 0;

/*     End of DSYGS2 */

} /* dsygs2_ */

#ifdef __cplusplus
	}
#endif
