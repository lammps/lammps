/* fortran/dlas2.f -- translated by f2c (version 20200916).
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

/* > \brief \b DLAS2 computes singular values of a 2-by-2 triangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAS2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlas2.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlas2.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlas2.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAS2( F, G, H, SSMIN, SSMAX ) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION   F, G, H, SSMAX, SSMIN */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAS2  computes the singular values of the 2-by-2 matrix */
/* >    [  F   G  ] */
/* >    [  0   H  ]. */
/* > On return, SSMIN is the smaller singular value and SSMAX is the */
/* > larger singular value. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] F */
/* > \verbatim */
/* >          F is DOUBLE PRECISION */
/* >          The (1,1) element of the 2-by-2 matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] G */
/* > \verbatim */
/* >          G is DOUBLE PRECISION */
/* >          The (1,2) element of the 2-by-2 matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] H */
/* > \verbatim */
/* >          H is DOUBLE PRECISION */
/* >          The (2,2) element of the 2-by-2 matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] SSMIN */
/* > \verbatim */
/* >          SSMIN is DOUBLE PRECISION */
/* >          The smaller singular value. */
/* > \endverbatim */
/* > */
/* > \param[out] SSMAX */
/* > \verbatim */
/* >          SSMAX is DOUBLE PRECISION */
/* >          The larger singular value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup OTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Barring over/underflow, all output quantities are correct to within */
/* >  a few units in the last place (ulps), even in the absence of a guard */
/* >  digit in addition/subtraction. */
/* > */
/* >  In IEEE arithmetic, the code works correctly if one matrix element is */
/* >  infinite. */
/* > */
/* >  Overflow will not occur unless the largest singular value itself */
/* >  overflows, or is within a few ulps of overflow. (On machines with */
/* >  partial overflow, like the Cray, overflow may occur if the largest */
/* >  singular value is within a factor of 2 of overflow.) */
/* > */
/* >  Underflow is harmless if underflow is gradual. Otherwise, results */
/* >  may correspond to a matrix modified by perturbations of size near */
/* >  the underflow threshold. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlas2_(doublereal *f, doublereal *g, doublereal *h__, 
	doublereal *ssmin, doublereal *ssmax)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal c__, fa, ga, ha, as, at, au, fhmn, fhmx;


/*  -- LAPACK auxiliary routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */

/*  ==================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    fa = abs(*f);
    ga = abs(*g);
    ha = abs(*h__);
    fhmn = min(fa,ha);
    fhmx = max(fa,ha);
    if (fhmn == 0.) {
	*ssmin = 0.;
	if (fhmx == 0.) {
	    *ssmax = ga;
	} else {
/* Computing 2nd power */
	    d__1 = min(fhmx,ga) / max(fhmx,ga);
	    *ssmax = max(fhmx,ga) * sqrt(d__1 * d__1 + 1.);
	}
    } else {
	if (ga < fhmx) {
	    as = fhmn / fhmx + 1.;
	    at = (fhmx - fhmn) / fhmx;
/* Computing 2nd power */
	    d__1 = ga / fhmx;
	    au = d__1 * d__1;
	    c__ = 2. / (sqrt(as * as + au) + sqrt(at * at + au));
	    *ssmin = fhmn * c__;
	    *ssmax = fhmx / c__;
	} else {
	    au = fhmx / ga;
	    if (au == 0.) {

/*              Avoid possible harmful underflow if exponent range */
/*              asymmetric (true SSMIN may not underflow even if */
/*              AU underflows) */

		*ssmin = fhmn * fhmx / ga;
		*ssmax = ga;
	    } else {
		as = fhmn / fhmx + 1.;
		at = (fhmx - fhmn) / fhmx;
/* Computing 2nd power */
		d__1 = as * au;
/* Computing 2nd power */
		d__2 = at * au;
		c__ = 1. / (sqrt(d__1 * d__1 + 1.) + sqrt(d__2 * d__2 + 1.));
		*ssmin = fhmn * c__ * au;
		*ssmin += *ssmin;
		*ssmax = ga / (c__ + c__);
	    }
	}
    }
    return 0;

/*     End of DLAS2 */

} /* dlas2_ */

#ifdef __cplusplus
	}
#endif
