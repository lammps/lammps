/* fortran/zswap.f -- translated by f2c (version 20200916).
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

/* > \brief \b ZSWAP */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSWAP(N,ZX,INCX,ZY,INCY) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,INCY,N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 ZX(*),ZY(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    ZSWAP interchanges two vectors. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         number of elements in input vector(s) */
/* > \endverbatim */
/* > */
/* > \param[in,out] ZX */
/* > \verbatim */
/* >          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >         storage spacing between elements of ZX */
/* > \endverbatim */
/* > */
/* > \param[in,out] ZY */
/* > \verbatim */
/* >          ZY is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCY ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >         storage spacing between elements of ZY */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complex16_blas_level1 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >     jack dongarra, 3/11/78. */
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zswap_(integer *n, doublecomplex *zx, integer *incx, 
	doublecomplex *zy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    integer i__, ix, iy;
    doublecomplex ztemp;


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
    /* Parameter adjustments */
    --zy;
    --zx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {

/*       code for both increments equal to 1 */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = i__;
	    ztemp.r = zx[i__2].r, ztemp.i = zx[i__2].i;
	    i__2 = i__;
	    i__3 = i__;
	    zx[i__2].r = zy[i__3].r, zx[i__2].i = zy[i__3].i;
	    i__2 = i__;
	    zy[i__2].r = ztemp.r, zy[i__2].i = ztemp.i;
	}
    } else {

/*       code for unequal increments or equal increments not equal */
/*         to 1 */

	ix = 1;
	iy = 1;
	if (*incx < 0) {
	    ix = (-(*n) + 1) * *incx + 1;
	}
	if (*incy < 0) {
	    iy = (-(*n) + 1) * *incy + 1;
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = ix;
	    ztemp.r = zx[i__2].r, ztemp.i = zx[i__2].i;
	    i__2 = ix;
	    i__3 = iy;
	    zx[i__2].r = zy[i__3].r, zx[i__2].i = zy[i__3].i;
	    i__2 = iy;
	    zy[i__2].r = ztemp.r, zy[i__2].i = ztemp.i;
	    ix += *incx;
	    iy += *incy;
	}
    }
    return 0;

/*     End of ZSWAP */

} /* zswap_ */

#ifdef __cplusplus
	}
#endif
