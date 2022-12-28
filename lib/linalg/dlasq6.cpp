/* fortran/dlasq6.f -- translated by f2c (version 20200916).
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

/* > \brief \b DLASQ6 computes one dqd transform in ping-pong form. Used by sbdsqr and sstegr. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASQ6 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasq6.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasq6.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasq6.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, */
/*                          DNM1, DNM2 ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            I0, N0, PP */
/*       DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DNM1, DNM2 */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASQ6 computes one dqd (shift equal to zero) transform in */
/* > ping-pong form, with protection against underflow and overflow. */
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

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dlasq6_(integer *i0, integer *n0, doublereal *z__,
        integer *pp, doublereal *dmin__, doublereal *dmin1, doublereal *dmin2,
         doublereal *dn, doublereal *dnm1, doublereal *dnm2)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    doublereal d__;
    integer j4, j4p2;
    doublereal emin, temp;
    extern doublereal dlamch_(char *, ftnlen);
    doublereal safmin;


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
/*     .. External Function .. */
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

    safmin = dlamch_((char *)"Safe minimum", (ftnlen)12);
    j4 = (*i0 << 2) + *pp - 3;
    emin = z__[j4 + 4];
    d__ = z__[j4];
    *dmin__ = d__;

    if (*pp == 0) {
        i__1 = *n0 - 3 << 2;
        for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
            z__[j4 - 2] = d__ + z__[j4 - 1];
            if (z__[j4 - 2] == 0.) {
                z__[j4] = 0.;
                d__ = z__[j4 + 1];
                *dmin__ = d__;
                emin = 0.;
            } else if (safmin * z__[j4 + 1] < z__[j4 - 2] && safmin * z__[j4
                    - 2] < z__[j4 + 1]) {
                temp = z__[j4 + 1] / z__[j4 - 2];
                z__[j4] = z__[j4 - 1] * temp;
                d__ *= temp;
            } else {
                z__[j4] = z__[j4 + 1] * (z__[j4 - 1] / z__[j4 - 2]);
                d__ = z__[j4 + 1] * (d__ / z__[j4 - 2]);
            }
            *dmin__ = min(*dmin__,d__);
/* Computing MIN */
            d__1 = emin, d__2 = z__[j4];
            emin = min(d__1,d__2);
/* L10: */
        }
    } else {
        i__1 = *n0 - 3 << 2;
        for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
            z__[j4 - 3] = d__ + z__[j4];
            if (z__[j4 - 3] == 0.) {
                z__[j4 - 1] = 0.;
                d__ = z__[j4 + 2];
                *dmin__ = d__;
                emin = 0.;
            } else if (safmin * z__[j4 + 2] < z__[j4 - 3] && safmin * z__[j4
                    - 3] < z__[j4 + 2]) {
                temp = z__[j4 + 2] / z__[j4 - 3];
                z__[j4 - 1] = z__[j4] * temp;
                d__ *= temp;
            } else {
                z__[j4 - 1] = z__[j4 + 2] * (z__[j4] / z__[j4 - 3]);
                d__ = z__[j4 + 2] * (d__ / z__[j4 - 3]);
            }
            *dmin__ = min(*dmin__,d__);
/* Computing MIN */
            d__1 = emin, d__2 = z__[j4 - 1];
            emin = min(d__1,d__2);
/* L20: */
        }
    }

/*     Unroll last two steps. */

    *dnm2 = d__;
    *dmin2 = *dmin__;
    j4 = (*n0 - 2 << 2) - *pp;
    j4p2 = j4 + (*pp << 1) - 1;
    z__[j4 - 2] = *dnm2 + z__[j4p2];
    if (z__[j4 - 2] == 0.) {
        z__[j4] = 0.;
        *dnm1 = z__[j4p2 + 2];
        *dmin__ = *dnm1;
        emin = 0.;
    } else if (safmin * z__[j4p2 + 2] < z__[j4 - 2] && safmin * z__[j4 - 2] <
            z__[j4p2 + 2]) {
        temp = z__[j4p2 + 2] / z__[j4 - 2];
        z__[j4] = z__[j4p2] * temp;
        *dnm1 = *dnm2 * temp;
    } else {
        z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
        *dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]);
    }
    *dmin__ = min(*dmin__,*dnm1);

    *dmin1 = *dmin__;
    j4 += 4;
    j4p2 = j4 + (*pp << 1) - 1;
    z__[j4 - 2] = *dnm1 + z__[j4p2];
    if (z__[j4 - 2] == 0.) {
        z__[j4] = 0.;
        *dn = z__[j4p2 + 2];
        *dmin__ = *dn;
        emin = 0.;
    } else if (safmin * z__[j4p2 + 2] < z__[j4 - 2] && safmin * z__[j4 - 2] <
            z__[j4p2 + 2]) {
        temp = z__[j4p2 + 2] / z__[j4 - 2];
        z__[j4] = z__[j4p2] * temp;
        *dn = *dnm1 * temp;
    } else {
        z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
        *dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]);
    }
    *dmin__ = min(*dmin__,*dn);

    z__[j4 + 2] = *dn;
    z__[(*n0 << 2) - *pp] = emin;
    return 0;

/*     End of DLASQ6 */

} /* dlasq6_ */

#ifdef __cplusplus
        }
#endif
