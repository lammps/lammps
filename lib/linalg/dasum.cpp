/* fortran/dasum.f -- translated by f2c (version 20200916).
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

/* > \brief \b DASUM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION DX(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    DASUM takes the sum of the absolute values. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         number of elements in input vector(s) */
/* > \endverbatim */
/* > */
/* > \param[in] DX */
/* > \verbatim */
/* >          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >         storage spacing between elements of DX */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup double_blas_level1 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >     jack dongarra, linpack, 3/11/78. */
/* >     modified 3/93 to return if incx .le. 0. */
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
doublereal dasum_(integer *n, doublereal *dx, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    integer i__, m, mp1;
    doublereal dtemp;
    integer nincx;


/*  -- Reference BLAS level1 routine -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    --dx;

    /* Function Body */
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0 || *incx <= 0) {
	return ret_val;
    }
    if (*incx == 1) {
/*        code for increment equal to 1 */


/*        clean-up loop */

	m = *n % 6;
	if (m != 0) {
	    i__1 = m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dtemp += (d__1 = dx[i__], abs(d__1));
	    }
	    if (*n < 6) {
		ret_val = dtemp;
		return ret_val;
	    }
	}
	mp1 = m + 1;
	i__1 = *n;
	for (i__ = mp1; i__ <= i__1; i__ += 6) {
	    dtemp = dtemp + (d__1 = dx[i__], abs(d__1)) + (d__2 = dx[i__ + 1],
		     abs(d__2)) + (d__3 = dx[i__ + 2], abs(d__3)) + (d__4 = 
		    dx[i__ + 3], abs(d__4)) + (d__5 = dx[i__ + 4], abs(d__5)) 
		    + (d__6 = dx[i__ + 5], abs(d__6));
	}
    } else {

/*        code for increment not equal to 1 */

	nincx = *n * *incx;
	i__1 = nincx;
	i__2 = *incx;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	    dtemp += (d__1 = dx[i__], abs(d__1));
	}
    }
    ret_val = dtemp;
    return ret_val;

/*     End of DASUM */

} /* dasum_ */

#ifdef __cplusplus
	}
#endif
