/* fortran/dcabs1.f -- translated by f2c (version 20200916).
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

/* > \brief \b DCABS1 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DCABS1(Z) */

/*       .. Scalar Arguments .. */
/*       COMPLEX*16 Z */
/*       .. */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DCABS1 computes |Re(.)| + |Im(.)| of a double complex number */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] Z */
/* > \verbatim */
/* >          Z is COMPLEX*16 */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup double_blas_level1 */

/*  ===================================================================== */
doublereal dcabs1_(doublecomplex *z__)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);


/*  -- Reference BLAS level1 routine -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. */
/*  ===================================================================== */

/*     .. Intrinsic Functions .. */

    ret_val = (d__1 = z__->r, abs(d__1)) + (d__2 = d_imag(z__), abs(d__2));
    return ret_val;

/*     End of DCABS1 */

} /* dcabs1_ */

#ifdef __cplusplus
	}
#endif
