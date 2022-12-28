/* fortran/dlaev2.f -- translated by f2c (version 20200916).
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

/* > \brief \b DLAEV2 computes the eigenvalues and eigenvectors of a 2-by-2 symmetric/Hermitian matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAEV2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaev2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaev2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaev2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 ) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION   A, B, C, CS1, RT1, RT2, SN1 */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix */
/* >    [  A   B  ] */
/* >    [  B   C  ]. */
/* > On return, RT1 is the eigenvalue of larger absolute value, RT2 is the */
/* > eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right */
/* > eigenvector for RT1, giving the decomposition */
/* > */
/* >    [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ] */
/* >    [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ]. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION */
/* >          The (1,1) element of the 2-by-2 matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION */
/* >          The (1,2) element and the conjugate of the (2,1) element of */
/* >          the 2-by-2 matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION */
/* >          The (2,2) element of the 2-by-2 matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] RT1 */
/* > \verbatim */
/* >          RT1 is DOUBLE PRECISION */
/* >          The eigenvalue of larger absolute value. */
/* > \endverbatim */
/* > */
/* > \param[out] RT2 */
/* > \verbatim */
/* >          RT2 is DOUBLE PRECISION */
/* >          The eigenvalue of smaller absolute value. */
/* > \endverbatim */
/* > */
/* > \param[out] CS1 */
/* > \verbatim */
/* >          CS1 is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[out] SN1 */
/* > \verbatim */
/* >          SN1 is DOUBLE PRECISION */
/* >          The vector (CS1, SN1) is a unit right eigenvector for RT1. */
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
/* >  RT1 is accurate to a few ulps barring over/underflow. */
/* > */
/* >  RT2 may be inaccurate if there is massive cancellation in the */
/* >  determinant A*C-B*B; higher precision or correctly rounded or */
/* >  correctly truncated arithmetic would be needed to compute RT2 */
/* >  accurately in all cases. */
/* > */
/* >  CS1 and SN1 are accurate to a few ulps barring over/underflow. */
/* > */
/* >  Overflow is possible only if RT1 is within a factor of 5 of overflow. */
/* >  Underflow is harmless if the input data is 0 or exceeds */
/* >     underflow_threshold / macheps. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlaev2_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *rt1, doublereal *rt2, doublereal *cs1, doublereal *sn1)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal ab, df, cs, ct, tb, sm, tn, rt, adf, acs;
    integer sgn1, sgn2;
    doublereal acmn, acmx;


/*  -- LAPACK auxiliary routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Compute the eigenvalues */

    sm = *a + *c__;
    df = *a - *c__;
    adf = abs(df);
    tb = *b + *b;
    ab = abs(tb);
    if (abs(*a) > abs(*c__)) {
	acmx = *a;
	acmn = *c__;
    } else {
	acmx = *c__;
	acmn = *a;
    }
    if (adf > ab) {
/* Computing 2nd power */
	d__1 = ab / adf;
	rt = adf * sqrt(d__1 * d__1 + 1.);
    } else if (adf < ab) {
/* Computing 2nd power */
	d__1 = adf / ab;
	rt = ab * sqrt(d__1 * d__1 + 1.);
    } else {

/*        Includes case AB=ADF=0 */

	rt = ab * sqrt(2.);
    }
    if (sm < 0.) {
	*rt1 = (sm - rt) * .5;
	sgn1 = -1;

/*        Order of execution important. */
/*        To get fully accurate smaller eigenvalue, */
/*        next line needs to be executed in higher precision. */

	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else if (sm > 0.) {
	*rt1 = (sm + rt) * .5;
	sgn1 = 1;

/*        Order of execution important. */
/*        To get fully accurate smaller eigenvalue, */
/*        next line needs to be executed in higher precision. */

	*rt2 = acmx / *rt1 * acmn - *b / *rt1 * *b;
    } else {

/*        Includes case RT1 = RT2 = 0 */

	*rt1 = rt * .5;
	*rt2 = rt * -.5;
	sgn1 = 1;
    }

/*     Compute the eigenvector */

    if (df >= 0.) {
	cs = df + rt;
	sgn2 = 1;
    } else {
	cs = df - rt;
	sgn2 = -1;
    }
    acs = abs(cs);
    if (acs > ab) {
	ct = -tb / cs;
	*sn1 = 1. / sqrt(ct * ct + 1.);
	*cs1 = ct * *sn1;
    } else {
	if (ab == 0.) {
	    *cs1 = 1.;
	    *sn1 = 0.;
	} else {
	    tn = -cs / tb;
	    *cs1 = 1. / sqrt(tn * tn + 1.);
	    *sn1 = tn * *cs1;
	}
    }
    if (sgn1 == sgn2) {
	tn = *cs1;
	*cs1 = -(*sn1);
	*sn1 = tn;
    }
    return 0;

/*     End of DLAEV2 */

} /* dlaev2_ */

#ifdef __cplusplus
	}
#endif
