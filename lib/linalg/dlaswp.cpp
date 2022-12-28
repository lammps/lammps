/* fortran/dlaswp.f -- translated by f2c (version 20200916).
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

/* > \brief \b DLASWP performs a series of row interchanges on a general rectangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASWP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaswp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaswp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaswp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, K1, K2, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASWP performs a series of row interchanges on the matrix A. */
/* > One row interchange is initiated for each of rows K1 through K2 of A. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the matrix of column dimension N to which the row */
/* >          interchanges will be applied. */
/* >          On exit, the permuted matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. */
/* > \endverbatim */
/* > */
/* > \param[in] K1 */
/* > \verbatim */
/* >          K1 is INTEGER */
/* >          The first element of IPIV for which a row interchange will */
/* >          be done. */
/* > \endverbatim */
/* > */
/* > \param[in] K2 */
/* > \verbatim */
/* >          K2 is INTEGER */
/* >          (K2-K1+1) is the number of elements of IPIV for which a row */
/* >          interchange will be done. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (K1+(K2-K1)*abs(INCX)) */
/* >          The vector of pivot indices. Only the elements in positions */
/* >          K1 through K1+(K2-K1)*abs(INCX) of IPIV are accessed. */
/* >          IPIV(K1+(K-K1)*abs(INCX)) = L implies rows K and L are to be */
/* >          interchanged. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >          The increment between successive values of IPIV. If INCX */
/* >          is negative, the pivots are applied in reverse order. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doubleOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Modified by */
/* >   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlaswp_(integer *n, doublereal *a, integer *lda, integer 
	*k1, integer *k2, integer *ipiv, integer *incx)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    integer i__, j, k, i1, i2, n32, ip, ix, ix0, inc;
    doublereal temp;


/*  -- LAPACK auxiliary routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Interchange row I with row IPIV(K1+(I-K1)*abs(INCX)) for each of rows */
/*     K1 through K2. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;

    /* Function Body */
    if (*incx > 0) {
	ix0 = *k1;
	i1 = *k1;
	i2 = *k2;
	inc = 1;
    } else if (*incx < 0) {
	ix0 = *k1 + (*k1 - *k2) * *incx;
	i1 = *k2;
	i2 = *k1;
	inc = -1;
    } else {
	return 0;
    }

    n32 = *n / 32 << 5;
    if (n32 != 0) {
	i__1 = n32;
	for (j = 1; j <= i__1; j += 32) {
	    ix = ix0;
	    i__2 = i2;
	    i__3 = inc;
	    for (i__ = i1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) 
		    {
		ip = ipiv[ix];
		if (ip != i__) {
		    i__4 = j + 31;
		    for (k = j; k <= i__4; ++k) {
			temp = a[i__ + k * a_dim1];
			a[i__ + k * a_dim1] = a[ip + k * a_dim1];
			a[ip + k * a_dim1] = temp;
/* L10: */
		    }
		}
		ix += *incx;
/* L20: */
	    }
/* L30: */
	}
    }
    if (n32 != *n) {
	++n32;
	ix = ix0;
	i__1 = i2;
	i__3 = inc;
	for (i__ = i1; i__3 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__3) {
	    ip = ipiv[ix];
	    if (ip != i__) {
		i__2 = *n;
		for (k = n32; k <= i__2; ++k) {
		    temp = a[i__ + k * a_dim1];
		    a[i__ + k * a_dim1] = a[ip + k * a_dim1];
		    a[ip + k * a_dim1] = temp;
/* L40: */
		}
	    }
	    ix += *incx;
/* L50: */
	}
    }

    return 0;

/*     End of DLASWP */

} /* dlaswp_ */

#ifdef __cplusplus
	}
#endif
