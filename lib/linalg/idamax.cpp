/* fortran/idamax.f -- translated by f2c (version 20200916).
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

/* > \brief \b IDAMAX */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       INTEGER FUNCTION IDAMAX(N,DX,INCX) */

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
/* >    IDAMAX finds the index of the first element having maximum absolute value. */
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

/* > \ingroup aux_blas */

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
integer idamax_(integer *n, doublereal *dx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;
    doublereal d__1;

    /* Local variables */
    integer i__, ix;
    doublereal dmax__;


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
    ret_val = 0;
    if (*n < 1 || *incx <= 0) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {

/*        code for increment equal to 1 */

	dmax__ = abs(dx[1]);
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    if ((d__1 = dx[i__], abs(d__1)) > dmax__) {
		ret_val = i__;
		dmax__ = (d__1 = dx[i__], abs(d__1));
	    }
	}
    } else {

/*        code for increment not equal to 1 */

	ix = 1;
	dmax__ = abs(dx[1]);
	ix += *incx;
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    if ((d__1 = dx[ix], abs(d__1)) > dmax__) {
		ret_val = i__;
		dmax__ = (d__1 = dx[ix], abs(d__1));
	    }
	    ix += *incx;
	}
    }
    return ret_val;

/*     End of IDAMAX */

} /* idamax_ */

#ifdef __cplusplus
	}
#endif
