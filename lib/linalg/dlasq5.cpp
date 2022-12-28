/* fortran/dlasq5.f -- translated by f2c (version 20200916).
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

/* > \brief \b DLASQ5 computes one dqds transform in ping-pong form. Used by sbdsqr and sstegr. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASQ5 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasq5.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasq5.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasq5.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASQ5( I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2, DN, */
/*                          DNM1, DNM2, IEEE, EPS ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            IEEE */
/*       INTEGER            I0, N0, PP */
/*       DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DNM1, DNM2, TAU, SIGMA, EPS */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASQ5 computes one dqds transform in ping-pong form, one */
/* > version for IEEE machines another for non IEEE machines. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] I0 */
/* > \verbatim */
/* >          I0 is INTEGER */
/* >        First index. */
/* > \endverbatim */
/* > */
/* > \param[in] N0 */
/* > \verbatim */
/* >          N0 is INTEGER */
/* >        Last index. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension ( 4*N ) */
/* >        Z holds the qd array. EMIN is stored in Z(4*N0) to avoid */
/* >        an extra argument. */
/* > \endverbatim */
/* > */
/* > \param[in] PP */
/* > \verbatim */
/* >          PP is INTEGER */
/* >        PP=0 for ping, PP=1 for pong. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is DOUBLE PRECISION */
/* >        This is the shift. */
/* > \endverbatim */
/* > */
/* > \param[in] SIGMA */
/* > \verbatim */
/* >          SIGMA is DOUBLE PRECISION */
/* >        This is the accumulated shift up to this step. */
/* > \endverbatim */
/* > */
/* > \param[out] DMIN */
/* > \verbatim */
/* >          DMIN is DOUBLE PRECISION */
/* >        Minimum value of d. */
/* > \endverbatim */
/* > */
/* > \param[out] DMIN1 */
/* > \verbatim */
/* >          DMIN1 is DOUBLE PRECISION */
/* >        Minimum value of d, excluding D( N0 ). */
/* > \endverbatim */
/* > */
/* > \param[out] DMIN2 */
/* > \verbatim */
/* >          DMIN2 is DOUBLE PRECISION */
/* >        Minimum value of d, excluding D( N0 ) and D( N0-1 ). */
/* > \endverbatim */
/* > */
/* > \param[out] DN */
/* > \verbatim */
/* >          DN is DOUBLE PRECISION */
/* >        d(N0), the last value of d. */
/* > \endverbatim */
/* > */
/* > \param[out] DNM1 */
/* > \verbatim */
/* >          DNM1 is DOUBLE PRECISION */
/* >        d(N0-1). */
/* > \endverbatim */
/* > */
/* > \param[out] DNM2 */
/* > \verbatim */
/* >          DNM2 is DOUBLE PRECISION */
/* >        d(N0-2). */
/* > \endverbatim */
/* > */
/* > \param[in] IEEE */
/* > \verbatim */
/* >          IEEE is LOGICAL */
/* >        Flag for IEEE or non IEEE arithmetic. */
/* > \endverbatim */
/* > */
/* > \param[in] EPS */
/* > \verbatim */
/* >          EPS is DOUBLE PRECISION */
/* >        This is the value of epsilon used. */
/* > \endverbatim */
/* > */
/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dlasq5_(integer *i0, integer *n0, doublereal *z__,
        integer *pp, doublereal *tau, doublereal *sigma, doublereal *dmin__,
        doublereal *dmin1, doublereal *dmin2, doublereal *dn, doublereal *
        dnm1, doublereal *dnm2, logical *ieee, doublereal *eps)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    doublereal d__;
    integer j4, j4p2;
    doublereal emin, temp, dthresh;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameter .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    /* Parameter adjustments */
    --z__;

    /* Function Body */
    if (*n0 - *i0 - 1 <= 0) {
        return 0;
    }

    dthresh = *eps * (*sigma + *tau);
    if (*tau < dthresh * .5) {
        *tau = 0.;
    }
    if (*tau != 0.) {
        j4 = (*i0 << 2) + *pp - 3;
        emin = z__[j4 + 4];
        d__ = z__[j4] - *tau;
        *dmin__ = d__;
        *dmin1 = -z__[j4];

        if (*ieee) {

/*        Code for IEEE arithmetic. */

            if (*pp == 0) {
                i__1 = *n0 - 3 << 2;
                for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
                    z__[j4 - 2] = d__ + z__[j4 - 1];
                    temp = z__[j4 + 1] / z__[j4 - 2];
                    d__ = d__ * temp - *tau;
                    *dmin__ = min(*dmin__,d__);
                    z__[j4] = z__[j4 - 1] * temp;
/* Computing MIN */
                    d__1 = z__[j4];
                    emin = min(d__1,emin);
/* L10: */
                }
            } else {
                i__1 = *n0 - 3 << 2;
                for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
                    z__[j4 - 3] = d__ + z__[j4];
                    temp = z__[j4 + 2] / z__[j4 - 3];
                    d__ = d__ * temp - *tau;
                    *dmin__ = min(*dmin__,d__);
                    z__[j4 - 1] = z__[j4] * temp;
/* Computing MIN */
                    d__1 = z__[j4 - 1];
                    emin = min(d__1,emin);
/* L20: */
                }
            }

/*        Unroll last two steps. */

            *dnm2 = d__;
            *dmin2 = *dmin__;
            j4 = (*n0 - 2 << 2) - *pp;
            j4p2 = j4 + (*pp << 1) - 1;
            z__[j4 - 2] = *dnm2 + z__[j4p2];
            z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
            *dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
            *dmin__ = min(*dmin__,*dnm1);

            *dmin1 = *dmin__;
            j4 += 4;
            j4p2 = j4 + (*pp << 1) - 1;
            z__[j4 - 2] = *dnm1 + z__[j4p2];
            z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
            *dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
            *dmin__ = min(*dmin__,*dn);

        } else {

/*        Code for non IEEE arithmetic. */

            if (*pp == 0) {
                i__1 = *n0 - 3 << 2;
                for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
                    z__[j4 - 2] = d__ + z__[j4 - 1];
                    if (d__ < 0.) {
                        return 0;
                    } else {
                        z__[j4] = z__[j4 + 1] * (z__[j4 - 1] / z__[j4 - 2]);
                        d__ = z__[j4 + 1] * (d__ / z__[j4 - 2]) - *tau;
                    }
                    *dmin__ = min(*dmin__,d__);
/* Computing MIN */
                    d__1 = emin, d__2 = z__[j4];
                    emin = min(d__1,d__2);
/* L30: */
                }
            } else {
                i__1 = *n0 - 3 << 2;
                for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
                    z__[j4 - 3] = d__ + z__[j4];
                    if (d__ < 0.) {
                        return 0;
                    } else {
                        z__[j4 - 1] = z__[j4 + 2] * (z__[j4] / z__[j4 - 3]);
                        d__ = z__[j4 + 2] * (d__ / z__[j4 - 3]) - *tau;
                    }
                    *dmin__ = min(*dmin__,d__);
/* Computing MIN */
                    d__1 = emin, d__2 = z__[j4 - 1];
                    emin = min(d__1,d__2);
/* L40: */
                }
            }

/*        Unroll last two steps. */

            *dnm2 = d__;
            *dmin2 = *dmin__;
            j4 = (*n0 - 2 << 2) - *pp;
            j4p2 = j4 + (*pp << 1) - 1;
            z__[j4 - 2] = *dnm2 + z__[j4p2];
            if (*dnm2 < 0.) {
                return 0;
            } else {
                z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
                *dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
            }
            *dmin__ = min(*dmin__,*dnm1);

            *dmin1 = *dmin__;
            j4 += 4;
            j4p2 = j4 + (*pp << 1) - 1;
            z__[j4 - 2] = *dnm1 + z__[j4p2];
            if (*dnm1 < 0.) {
                return 0;
            } else {
                z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
                *dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
            }
            *dmin__ = min(*dmin__,*dn);

        }
    } else {
/*     This is the version that sets d's to zero if they are small enough */
        j4 = (*i0 << 2) + *pp - 3;
        emin = z__[j4 + 4];
        d__ = z__[j4] - *tau;
        *dmin__ = d__;
        *dmin1 = -z__[j4];
        if (*ieee) {

/*     Code for IEEE arithmetic. */

            if (*pp == 0) {
                i__1 = *n0 - 3 << 2;
                for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
                    z__[j4 - 2] = d__ + z__[j4 - 1];
                    temp = z__[j4 + 1] / z__[j4 - 2];
                    d__ = d__ * temp - *tau;
                    if (d__ < dthresh) {
                        d__ = 0.;
                    }
                    *dmin__ = min(*dmin__,d__);
                    z__[j4] = z__[j4 - 1] * temp;
/* Computing MIN */
                    d__1 = z__[j4];
                    emin = min(d__1,emin);
/* L50: */
                }
            } else {
                i__1 = *n0 - 3 << 2;
                for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
                    z__[j4 - 3] = d__ + z__[j4];
                    temp = z__[j4 + 2] / z__[j4 - 3];
                    d__ = d__ * temp - *tau;
                    if (d__ < dthresh) {
                        d__ = 0.;
                    }
                    *dmin__ = min(*dmin__,d__);
                    z__[j4 - 1] = z__[j4] * temp;
/* Computing MIN */
                    d__1 = z__[j4 - 1];
                    emin = min(d__1,emin);
/* L60: */
                }
            }

/*     Unroll last two steps. */

            *dnm2 = d__;
            *dmin2 = *dmin__;
            j4 = (*n0 - 2 << 2) - *pp;
            j4p2 = j4 + (*pp << 1) - 1;
            z__[j4 - 2] = *dnm2 + z__[j4p2];
            z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
            *dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
            *dmin__ = min(*dmin__,*dnm1);

            *dmin1 = *dmin__;
            j4 += 4;
            j4p2 = j4 + (*pp << 1) - 1;
            z__[j4 - 2] = *dnm1 + z__[j4p2];
            z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
            *dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
            *dmin__ = min(*dmin__,*dn);

        } else {

/*     Code for non IEEE arithmetic. */

            if (*pp == 0) {
                i__1 = *n0 - 3 << 2;
                for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
                    z__[j4 - 2] = d__ + z__[j4 - 1];
                    if (d__ < 0.) {
                        return 0;
                    } else {
                        z__[j4] = z__[j4 + 1] * (z__[j4 - 1] / z__[j4 - 2]);
                        d__ = z__[j4 + 1] * (d__ / z__[j4 - 2]) - *tau;
                    }
                    if (d__ < dthresh) {
                        d__ = 0.;
                    }
                    *dmin__ = min(*dmin__,d__);
/* Computing MIN */
                    d__1 = emin, d__2 = z__[j4];
                    emin = min(d__1,d__2);
/* L70: */
                }
            } else {
                i__1 = *n0 - 3 << 2;
                for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
                    z__[j4 - 3] = d__ + z__[j4];
                    if (d__ < 0.) {
                        return 0;
                    } else {
                        z__[j4 - 1] = z__[j4 + 2] * (z__[j4] / z__[j4 - 3]);
                        d__ = z__[j4 + 2] * (d__ / z__[j4 - 3]) - *tau;
                    }
                    if (d__ < dthresh) {
                        d__ = 0.;
                    }
                    *dmin__ = min(*dmin__,d__);
/* Computing MIN */
                    d__1 = emin, d__2 = z__[j4 - 1];
                    emin = min(d__1,d__2);
/* L80: */
                }
            }

/*     Unroll last two steps. */

            *dnm2 = d__;
            *dmin2 = *dmin__;
            j4 = (*n0 - 2 << 2) - *pp;
            j4p2 = j4 + (*pp << 1) - 1;
            z__[j4 - 2] = *dnm2 + z__[j4p2];
            if (*dnm2 < 0.) {
                return 0;
            } else {
                z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
                *dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
            }
            *dmin__ = min(*dmin__,*dnm1);

            *dmin1 = *dmin__;
            j4 += 4;
            j4p2 = j4 + (*pp << 1) - 1;
            z__[j4 - 2] = *dnm1 + z__[j4p2];
            if (*dnm1 < 0.) {
                return 0;
            } else {
                z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
                *dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
            }
            *dmin__ = min(*dmin__,*dn);

        }
    }

    z__[j4 + 2] = *dn;
    z__[(*n0 << 2) - *pp] = emin;
    return 0;

/*     End of DLASQ5 */

} /* dlasq5_ */

#ifdef __cplusplus
        }
#endif
