/* fortran/dlansy.f -- translated by f2c (version 20200916).
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

/* > \brief \b DLANSY returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a real symmetric matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLANSY + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlansy.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlansy.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlansy.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DLANSY( NORM, UPLO, N, A, LDA, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM, UPLO */
/*       INTEGER            LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLANSY  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > real symmetric matrix A. */
/* > \endverbatim */
/* > */
/* > \return DLANSY */
/* > \verbatim */
/* > */
/* >    DLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
/* >             ( */
/* >             ( norm1(A),         NORM = '1', 'O' or 'o' */
/* >             ( */
/* >             ( normI(A),         NORM = 'I' or 'i' */
/* >             ( */
/* >             ( normF(A),         NORM = 'F', 'f', 'E' or 'e' */
/* > */
/* > where  norm1  denotes the  one norm of a matrix (maximum column sum), */
/* > normI  denotes the  infinity norm  of a matrix  (maximum row sum) and */
/* > normF  denotes the  Frobenius norm of a matrix (square root of sum of */
/* > squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] NORM */
/* > \verbatim */
/* >          NORM is CHARACTER*1 */
/* >          Specifies the value to be returned in DLANSY as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the upper or lower triangular part of the */
/* >          symmetric matrix A is to be referenced. */
/* >          = 'U':  Upper triangular part of A is referenced */
/* >          = 'L':  Lower triangular part of A is referenced */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0.  When N = 0, DLANSY is */
/* >          set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          The symmetric matrix A.  If UPLO = 'U', the leading n by n */
/* >          upper triangular part of A contains the upper triangular part */
/* >          of the matrix A, and the strictly lower triangular part of A */
/* >          is not referenced.  If UPLO = 'L', the leading n by n lower */
/* >          triangular part of A contains the lower triangular part of */
/* >          the matrix A, and the strictly upper triangular part of A is */
/* >          not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(N,1). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)), */
/* >          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise, */
/* >          WORK is not referenced. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doubleSYauxiliary */

/*  ===================================================================== */
doublereal dlansy_(char *norm, char *uplo, integer *n, doublereal *a, integer 
	*lda, doublereal *work, ftnlen norm_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__, j;
    doublereal sum, absa, scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    doublereal value;
    extern logical disnan_(doublereal *);
    extern /* Subroutine */ int dlassq_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);


/*  -- LAPACK auxiliary routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --work;

    /* Function Body */
    if (*n == 0) {
	value = 0.;
    } else if (lsame_(norm, (char *)"M", (ftnlen)1, (ftnlen)1)) {

/*        Find max(abs(A(i,j))). */

	value = 0.;
	if (lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1)) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    sum = (d__1 = a[i__ + j * a_dim1], abs(d__1));
		    if (value < sum || disnan_(&sum)) {
			value = sum;
		    }
/* L10: */
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = j; i__ <= i__2; ++i__) {
		    sum = (d__1 = a[i__ + j * a_dim1], abs(d__1));
		    if (value < sum || disnan_(&sum)) {
			value = sum;
		    }
/* L30: */
		}
/* L40: */
	    }
	}
    } else if (lsame_(norm, (char *)"I", (ftnlen)1, (ftnlen)1) || lsame_(norm, (char *)"O", (
	    ftnlen)1, (ftnlen)1) || *(unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is symmetric). */

	value = 0.;
	if (lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1)) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		sum = 0.;
		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    absa = (d__1 = a[i__ + j * a_dim1], abs(d__1));
		    sum += absa;
		    work[i__] += absa;
/* L50: */
		}
		work[j] = sum + (d__1 = a[j + j * a_dim1], abs(d__1));
/* L60: */
	    }
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		sum = work[i__];
		if (value < sum || disnan_(&sum)) {
		    value = sum;
		}
/* L70: */
	    }
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		work[i__] = 0.;
/* L80: */
	    }
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		sum = work[j] + (d__1 = a[j + j * a_dim1], abs(d__1));
		i__2 = *n;
		for (i__ = j + 1; i__ <= i__2; ++i__) {
		    absa = (d__1 = a[i__ + j * a_dim1], abs(d__1));
		    sum += absa;
		    work[i__] += absa;
/* L90: */
		}
		if (value < sum || disnan_(&sum)) {
		    value = sum;
		}
/* L100: */
	    }
	}
    } else if (lsame_(norm, (char *)"F", (ftnlen)1, (ftnlen)1) || lsame_(norm, (char *)"E", (
	    ftnlen)1, (ftnlen)1)) {

/*        Find normF(A). */

	scale = 0.;
	sum = 1.;
	if (lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1)) {
	    i__1 = *n;
	    for (j = 2; j <= i__1; ++j) {
		i__2 = j - 1;
		dlassq_(&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
/* L110: */
	    }
	} else {
	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n - j;
		dlassq_(&i__2, &a[j + 1 + j * a_dim1], &c__1, &scale, &sum);
/* L120: */
	    }
	}
	sum *= 2;
	i__1 = *lda + 1;
	dlassq_(n, &a[a_offset], &i__1, &scale, &sum);
	value = scale * sqrt(sum);
    }

    ret_val = value;
    return ret_val;

/*     End of DLANSY */

} /* dlansy_ */

#ifdef __cplusplus
	}
#endif
