/* fortran/zdotc.f -- translated by f2c (version 20200916).
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

/* > \brief \b ZDOTC */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       COMPLEX*16 FUNCTION ZDOTC(N,ZX,INCX,ZY,INCY) */

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
/* > ZDOTC forms the dot product of two complex vectors */
/* >      ZDOTC = X^H * Y */
/* > */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         number of elements in input vector(s) */
/* > \endverbatim */
/* > */
/* > \param[in] ZX */
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
/* > \param[in] ZY */
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
/* Double Complex */ VOID zdotc_(doublecomplex * ret_val, integer *n, 
	doublecomplex *zx, integer *incx, doublecomplex *zy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

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
/*     .. Intrinsic Functions .. */
/*     .. */
    /* Parameter adjustments */
    --zy;
    --zx;

    /* Function Body */
    ztemp.r = 0., ztemp.i = 0.;
     ret_val->r = 0.,  ret_val->i = 0.;
    if (*n <= 0) {
	return ;
    }
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    d_cnjg(&z__3, &zx[i__]);
	    i__2 = i__;
	    z__2.r = z__3.r * zy[i__2].r - z__3.i * zy[i__2].i, z__2.i = 
		    z__3.r * zy[i__2].i + z__3.i * zy[i__2].r;
	    z__1.r = ztemp.r + z__2.r, z__1.i = ztemp.i + z__2.i;
	    ztemp.r = z__1.r, ztemp.i = z__1.i;
	}
    } else {

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

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
	    d_cnjg(&z__3, &zx[ix]);
	    i__2 = iy;
	    z__2.r = z__3.r * zy[i__2].r - z__3.i * zy[i__2].i, z__2.i = 
		    z__3.r * zy[i__2].i + z__3.i * zy[i__2].r;
	    z__1.r = ztemp.r + z__2.r, z__1.i = ztemp.i + z__2.i;
	    ztemp.r = z__1.r, ztemp.i = z__1.i;
	    ix += *incx;
	    iy += *incy;
	}
    }
     ret_val->r = ztemp.r,  ret_val->i = ztemp.i;
    return ;

/*     End of ZDOTC */

} /* zdotc_ */

#ifdef __cplusplus
	}
#endif
