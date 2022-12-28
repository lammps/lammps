/* fortran/dlapy3.f -- translated by f2c (version 20200916).
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

/* > \brief \b DLAPY3 returns sqrt(x2+y2+z2). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAPY3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlapy3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlapy3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlapy3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       DOUBLE PRECISION FUNCTION DLAPY3( X, Y, Z ) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION   X, Y, Z */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause */
/* > unnecessary overflow and unnecessary underflow. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] Y */
/* > \verbatim */
/* >          Y is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION */
/* >          X, Y and Z specify the values x, y and z. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup OTHERauxiliary */

/*  ===================================================================== */
doublereal dlapy3_(doublereal *x, doublereal *y, doublereal *z__)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal w, xabs, yabs, zabs;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal hugeval;


/*  -- LAPACK auxiliary routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
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
/*     .. Executable Statements .. */

    hugeval = dlamch_((char *)"Overflow", (ftnlen)8);
    xabs = abs(*x);
    yabs = abs(*y);
    zabs = abs(*z__);
/* Computing MAX */
    d__1 = max(xabs,yabs);
    w = max(d__1,zabs);
    if (w == 0. || w > hugeval) {
/*     W can be zero for max(0,nan,0) */
/*     adding all three entries together will make sure */
/*     NaN will not disappear. */
	ret_val = xabs + yabs + zabs;
    } else {
/* Computing 2nd power */
	d__1 = xabs / w;
/* Computing 2nd power */
	d__2 = yabs / w;
/* Computing 2nd power */
	d__3 = zabs / w;
	ret_val = w * sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    }
    return ret_val;

/*     End of DLAPY3 */

} /* dlapy3_ */

#ifdef __cplusplus
	}
#endif
