/* fortran/zdrot.f -- translated by f2c (version 20200916).
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

/* > \brief \b ZDROT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZDROT( N, ZX, INCX, ZY, INCY, C, S ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INCX, INCY, N */
/*       DOUBLE PRECISION   C, S */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         ZX( * ), ZY( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > Applies a plane rotation, where the cos and sin (c and s) are real */
/* > and the vectors cx and cy are complex. */
/* > jack dongarra, linpack, 3/11/78. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the order of the vectors cx and cy. */
/* >           N must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in,out] ZX */
/* > \verbatim */
/* >          ZX is COMPLEX*16 array, dimension at least */
/* >           ( 1 + ( N - 1 )*abs( INCX ) ). */
/* >           Before entry, the incremented array ZX must contain the n */
/* >           element vector cx. On exit, ZX is overwritten by the updated */
/* >           vector cx. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >           On entry, INCX specifies the increment for the elements of */
/* >           ZX. INCX must not be zero. */
/* > \endverbatim */
/* > */
/* > \param[in,out] ZY */
/* > \verbatim */
/* >          ZY is COMPLEX*16 array, dimension at least */
/* >           ( 1 + ( N - 1 )*abs( INCY ) ). */
/* >           Before entry, the incremented array ZY must contain the n */
/* >           element vector cy. On exit, ZY is overwritten by the updated */
/* >           vector cy. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >           On entry, INCY specifies the increment for the elements of */
/* >           ZY. INCY must not be zero. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION */
/* >           On entry, C specifies the cosine, cos. */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* >          S is DOUBLE PRECISION */
/* >           On entry, S specifies the sine, sin. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complex16_blas_level1 */

/*  ===================================================================== */
/* Subroutine */ int zdrot_(integer *n, doublecomplex *zx, integer *incx, 
	doublecomplex *zy, integer *incy, doublereal *c__, doublereal *s)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2, z__3;

    /* Local variables */
    integer i__, ix, iy;
    doublecomplex ctemp;


/*  -- Reference BLAS level1 routine -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --zy;
    --zx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {

/*        code for both increments equal to 1 */

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = i__;
	    z__2.r = *c__ * zx[i__2].r, z__2.i = *c__ * zx[i__2].i;
	    i__3 = i__;
	    z__3.r = *s * zy[i__3].r, z__3.i = *s * zy[i__3].i;
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	    ctemp.r = z__1.r, ctemp.i = z__1.i;
	    i__2 = i__;
	    i__3 = i__;
	    z__2.r = *c__ * zy[i__3].r, z__2.i = *c__ * zy[i__3].i;
	    i__4 = i__;
	    z__3.r = *s * zx[i__4].r, z__3.i = *s * zx[i__4].i;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    zy[i__2].r = z__1.r, zy[i__2].i = z__1.i;
	    i__2 = i__;
	    zx[i__2].r = ctemp.r, zx[i__2].i = ctemp.i;
	}
    } else {

/*        code for unequal increments or equal increments not equal */
/*          to 1 */

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
	    z__2.r = *c__ * zx[i__2].r, z__2.i = *c__ * zx[i__2].i;
	    i__3 = iy;
	    z__3.r = *s * zy[i__3].r, z__3.i = *s * zy[i__3].i;
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	    ctemp.r = z__1.r, ctemp.i = z__1.i;
	    i__2 = iy;
	    i__3 = iy;
	    z__2.r = *c__ * zy[i__3].r, z__2.i = *c__ * zy[i__3].i;
	    i__4 = ix;
	    z__3.r = *s * zx[i__4].r, z__3.i = *s * zx[i__4].i;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    zy[i__2].r = z__1.r, zy[i__2].i = z__1.i;
	    i__2 = ix;
	    zx[i__2].r = ctemp.r, zx[i__2].i = ctemp.i;
	    ix += *incx;
	    iy += *incy;
	}
    }
    return 0;

/*     End of ZDROT */

} /* zdrot_ */

#ifdef __cplusplus
	}
#endif
