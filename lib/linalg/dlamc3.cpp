/* static/dlamc3.f -- translated by f2c (version 20200916).
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

/* > \brief \b DLAMC3 */
/* > \details */
/* > \b Purpose: */
/* > \verbatim */
/* > DLAMC3  is intended to force  A  and  B  to be stored prior to doing */
/* > the addition of  A  and  B ,  for use in situations where optimizers */
/* > might hold one of these in a register. */
/* > \endverbatim */
/* > \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ.
of Colorado Denver and NAG Ltd.. */
/* > \date December 2016 */
/* > \ingroup auxOTHERauxiliary */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is a DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is a DOUBLE PRECISION */
/* >          The values A and B. */
/* > \endverbatim */
/* > */
doublereal dlamc3_(doublereal *a, doublereal *b)
{
    /* System generated locals */
    doublereal ret_val;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2010 */

/*     .. Scalar Arguments .. */
/*     .. */
/* ===================================================================== */

/*     .. Executable Statements .. */

    ret_val = *a + *b;

    return ret_val;

/*     End of DLAMC3 */

} /* dlamc3_ */

#ifdef __cplusplus
        }
#endif
