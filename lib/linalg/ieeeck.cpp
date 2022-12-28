/* fortran/ieeeck.f -- translated by f2c (version 20200916).
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

/* > \brief \b IEEECK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download IEEECK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ieeeck.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ieeeck.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ieeeck.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            ISPEC */
/*       REAL               ONE, ZERO */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > IEEECK is called from the ILAENV to verify that Infinity and */
/* > possibly NaN arithmetic is safe (i.e. will not trap). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ISPEC */
/* > \verbatim */
/* >          ISPEC is INTEGER */
/* >          Specifies whether to test just for infinity arithmetic */
/* >          or whether to test for infinity and NaN arithmetic. */
/* >          = 0: Verify infinity arithmetic only. */
/* >          = 1: Verify infinity and NaN arithmetic. */
/* > \endverbatim */
/* > */
/* > \param[in] ZERO */
/* > \verbatim */
/* >          ZERO is REAL */
/* >          Must contain the value 0.0 */
/* >          This is passed to prevent the compiler from optimizing */
/* >          away this code. */
/* > \endverbatim */
/* > */
/* > \param[in] ONE */
/* > \verbatim */
/* >          ONE is REAL */
/* >          Must contain the value 1.0 */
/* >          This is passed to prevent the compiler from optimizing */
/* >          away this code. */
/* > */
/* >  RETURN VALUE:  INTEGER */
/* >          = 0:  Arithmetic failed to produce the correct answers */
/* >          = 1:  Arithmetic produced the correct answers */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup OTHERauxiliary */

/*  ===================================================================== */
integer ieeeck_(integer *ispec, real *zero, real *one)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    real nan1, nan2, nan3, nan4, nan5, nan6, neginf, posinf, negzro, newzro;


/*  -- LAPACK auxiliary routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */
    ret_val = 1;

    posinf = *one / *zero;
    if (posinf <= *one) {
        ret_val = 0;
        return ret_val;
    }

    neginf = -(*one) / *zero;
    if (neginf >= *zero) {
        ret_val = 0;
        return ret_val;
    }

    negzro = *one / (neginf + *one);
    if (negzro != *zero) {
        ret_val = 0;
        return ret_val;
    }

    neginf = *one / negzro;
    if (neginf >= *zero) {
        ret_val = 0;
        return ret_val;
    }

    newzro = negzro + *zero;
    if (newzro != *zero) {
        ret_val = 0;
        return ret_val;
    }

    posinf = *one / newzro;
    if (posinf <= *one) {
        ret_val = 0;
        return ret_val;
    }

    neginf *= posinf;
    if (neginf >= *zero) {
        ret_val = 0;
        return ret_val;
    }

    posinf *= posinf;
    if (posinf <= *one) {
        ret_val = 0;
        return ret_val;
    }




/*     Return if we were only asked to check infinity arithmetic */

    if (*ispec == 0) {
        return ret_val;
    }

    nan1 = posinf + neginf;

    nan2 = posinf / neginf;

    nan3 = posinf / posinf;

    nan4 = posinf * *zero;

    nan5 = neginf * negzro;

    nan6 = nan5 * *zero;

    if (nan1 == nan1) {
        ret_val = 0;
        return ret_val;
    }

    if (nan2 == nan2) {
        ret_val = 0;
        return ret_val;
    }

    if (nan3 == nan3) {
        ret_val = 0;
        return ret_val;
    }

    if (nan4 == nan4) {
        ret_val = 0;
        return ret_val;
    }

    if (nan5 == nan5) {
        ret_val = 0;
        return ret_val;
    }

    if (nan6 == nan6) {
        ret_val = 0;
        return ret_val;
    }

    return ret_val;
} /* ieeeck_ */

#ifdef __cplusplus
        }
#endif
