/* fortran/iparmq.f -- translated by f2c (version 20200916).
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

/* > \brief \b IPARMQ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download IPARMQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iparmq.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iparmq.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iparmq.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            IHI, ILO, ISPEC, LWORK, N */
/*       CHARACTER          NAME*( * ), OPTS*( * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >      This program sets problem and machine dependent parameters */
/* >      useful for xHSEQR and related subroutines for eigenvalue */
/* >      problems. It is called whenever */
/* >      IPARMQ is called with 12 <= ISPEC <= 16 */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ISPEC */
/* > \verbatim */
/* >          ISPEC is INTEGER */
/* >              ISPEC specifies which tunable parameter IPARMQ should */
/* >              return. */
/* > */
/* >              ISPEC=12: (INMIN)  Matrices of order nmin or less */
/* >                        are sent directly to xLAHQR, the implicit */
/* >                        double shift QR algorithm.  NMIN must be */
/* >                        at least 11. */
/* > */
/* >              ISPEC=13: (INWIN)  Size of the deflation window. */
/* >                        This is best set greater than or equal to */
/* >                        the number of simultaneous shifts NS. */
/* >                        Larger matrices benefit from larger deflation */
/* >                        windows. */
/* > */
/* >              ISPEC=14: (INIBL) Determines when to stop nibbling and */
/* >                        invest in an (expensive) multi-shift QR sweep. */
/* >                        If the aggressive early deflation subroutine */
/* >                        finds LD converged eigenvalues from an order */
/* >                        NW deflation window and LD > (NW*NIBBLE)/100, */
/* >                        then the next QR sweep is skipped and early */
/* >                        deflation is applied immediately to the */
/* >                        remaining active diagonal block.  Setting */
/* >                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a */
/* >                        multi-shift QR sweep whenever early deflation */
/* >                        finds a converged eigenvalue.  Setting */
/* >                        IPARMQ(ISPEC=14) greater than or equal to 100 */
/* >                        prevents TTQRE from skipping a multi-shift */
/* >                        QR sweep. */
/* > */
/* >              ISPEC=15: (NSHFTS) The number of simultaneous shifts in */
/* >                        a multi-shift QR iteration. */
/* > */
/* >              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the */
/* >                        following meanings. */
/* >                        0:  During the multi-shift QR/QZ sweep, */
/* >                            blocked eigenvalue reordering, blocked */
/* >                            Hessenberg-triangular reduction, */
/* >                            reflections and/or rotations are not */
/* >                            accumulated when updating the */
/* >                            far-from-diagonal matrix entries. */
/* >                        1:  During the multi-shift QR/QZ sweep, */
/* >                            blocked eigenvalue reordering, blocked */
/* >                            Hessenberg-triangular reduction, */
/* >                            reflections and/or rotations are */
/* >                            accumulated, and matrix-matrix */
/* >                            multiplication is used to update the */
/* >                            far-from-diagonal matrix entries. */
/* >                        2:  During the multi-shift QR/QZ sweep, */
/* >                            blocked eigenvalue reordering, blocked */
/* >                            Hessenberg-triangular reduction, */
/* >                            reflections and/or rotations are */
/* >                            accumulated, and 2-by-2 block structure */
/* >                            is exploited during matrix-matrix */
/* >                            multiplies. */
/* >                        (If xTRMM is slower than xGEMM, then */
/* >                        IPARMQ(ISPEC=16)=1 may be more efficient than */
/* >                        IPARMQ(ISPEC=16)=2 despite the greater level of */
/* >                        arithmetic work implied by the latter choice.) */
/* > */
/* >              ISPEC=17: (ICOST) An estimate of the relative cost of flops */
/* >                        within the near-the-diagonal shift chase compared */
/* >                        to flops within the BLAS calls of a QZ sweep. */
/* > \endverbatim */
/* > */
/* > \param[in] NAME */
/* > \verbatim */
/* >          NAME is CHARACTER string */
/* >               Name of the calling subroutine */
/* > \endverbatim */
/* > */
/* > \param[in] OPTS */
/* > \verbatim */
/* >          OPTS is CHARACTER string */
/* >               This is a concatenation of the string arguments to */
/* >               TTQRE. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >               N is the order of the Hessenberg matrix H. */
/* > \endverbatim */
/* > */
/* > \param[in] ILO */
/* > \verbatim */
/* >          ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* >          IHI is INTEGER */
/* >               It is assumed that H is already upper triangular */
/* >               in rows and columns 1:ILO-1 and IHI+1:N. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >               The amount of workspace available. */
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
/* >       Little is known about how best to choose these parameters. */
/* >       It is possible to use different values of the parameters */
/* >       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR. */
/* > */
/* >       It is probably best to choose different parameters for */
/* >       different matrices and different parameters at different */
/* >       times during the iteration, but this has not been */
/* >       implemented --- yet. */
/* > */
/* > */
/* >       The best choices of most of the parameters depend */
/* >       in an ill-understood way on the relative execution */
/* >       rate of xLAQR3 and xLAQR5 and on the nature of each */
/* >       particular eigenvalue problem.  Experiment may be the */
/* >       only practical way to determine which choices are most */
/* >       effective. */
/* > */
/* >       Following is a list of default values supplied by IPARMQ. */
/* >       These defaults may be adjusted in order to attain better */
/* >       performance in any particular computational environment. */
/* > */
/* >       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point. */
/* >                        Default: 75. (Must be at least 11.) */
/* > */
/* >       IPARMQ(ISPEC=13) Recommended deflation window size. */
/* >                        This depends on ILO, IHI and NS, the */
/* >                        number of simultaneous shifts returned */
/* >                        by IPARMQ(ISPEC=15).  The default for */
/* >                        (IHI-ILO+1) <= 500 is NS.  The default */
/* >                        for (IHI-ILO+1) > 500 is 3*NS/2. */
/* > */
/* >       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14. */
/* > */
/* >       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS. */
/* >                        a multi-shift QR iteration. */
/* > */
/* >                        If IHI-ILO+1 is ... */
/* > */
/* >                        greater than      ...but less    ... the */
/* >                        or equal to ...      than        default is */
/* > */
/* >                                0               30       NS =   2+ */
/* >                               30               60       NS =   4+ */
/* >                               60              150       NS =  10 */
/* >                              150              590       NS =  ** */
/* >                              590             3000       NS =  64 */
/* >                             3000             6000       NS = 128 */
/* >                             6000             infinity   NS = 256 */
/* > */
/* >                    (+)  By default matrices of this order are */
/* >                         passed to the implicit double shift routine */
/* >                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These */
/* >                         values of NS are used only in case of a rare */
/* >                         xLAHQR failure. */
/* > */
/* >                    (**) The asterisks (**) indicate an ad-hoc */
/* >                         function increasing from 10 to 64. */
/* > */
/* >       IPARMQ(ISPEC=16) Select structured matrix multiply. */
/* >                        (See ISPEC=16 above for details.) */
/* >                        Default: 3. */
/* > */
/* >       IPARMQ(ISPEC=17) Relative cost heuristic for blocksize selection. */
/* >                        Expressed as a percentage. */
/* >                        Default: 10. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
integer iparmq_(integer *ispec, char *name__, char *opts, integer *n, integer
        *ilo, integer *ihi, integer *lwork, ftnlen name_len, ftnlen opts_len)
{
    /* System generated locals */
    integer ret_val, i__1, i__2;
    real r__1;

    /* Builtin functions */
    double log(doublereal);
    integer i_lmp_nint(real *);
    /* Subroutine */ int s_lmp_copy(char *, char *, ftnlen, ftnlen);
    integer s_lmp_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    integer i__, ic, nh, ns, iz;
    char subnam[6];


/*  -- LAPACK auxiliary routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */

/*  ================================================================ */
/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */
    if (*ispec == 15 || *ispec == 13 || *ispec == 16) {

/*        ==== Set the number simultaneous shifts ==== */

        nh = *ihi - *ilo + 1;
        ns = 2;
        if (nh >= 30) {
            ns = 4;
        }
        if (nh >= 60) {
            ns = 10;
        }
        if (nh >= 150) {
/* Computing MAX */
            r__1 = log((real) nh) / log((float)2.);
            i__1 = 10, i__2 = nh / i_lmp_nint(&r__1);
            ns = max(i__1,i__2);
        }
        if (nh >= 590) {
            ns = 64;
        }
        if (nh >= 3000) {
            ns = 128;
        }
        if (nh >= 6000) {
            ns = 256;
        }
/* Computing MAX */
        i__1 = 2, i__2 = ns - ns % 2;
        ns = max(i__1,i__2);
    }

    if (*ispec == 12) {


/*        ===== Matrices of order smaller than NMIN get sent */
/*        .     to xLAHQR, the classic double shift algorithm. */
/*        .     This must be at least 11. ==== */

        ret_val = 75;

    } else if (*ispec == 14) {

/*        ==== INIBL: skip a multi-shift qr iteration and */
/*        .    whenever aggressive early deflation finds */
/*        .    at least (NIBBLE*(window size)/100) deflations. ==== */

        ret_val = 14;

    } else if (*ispec == 15) {

/*        ==== NSHFTS: The number of simultaneous shifts ===== */

        ret_val = ns;

    } else if (*ispec == 13) {

/*        ==== NW: deflation window size.  ==== */

        if (nh <= 500) {
            ret_val = ns;
        } else {
            ret_val = ns * 3 / 2;
        }

    } else if (*ispec == 16) {

/*        ==== IACC22: Whether to accumulate reflections */
/*        .     before updating the far-from-diagonal elements */
/*        .     and whether to use 2-by-2 block structure while */
/*        .     doing it.  A small amount of work could be saved */
/*        .     by making this choice dependent also upon the */
/*        .     NH=IHI-ILO+1. */


/*        Convert NAME to upper case if the first character is lower case. */

        ret_val = 0;
        s_lmp_copy(subnam, name__, (ftnlen)6, name_len);
        ic = *(unsigned char *)subnam;
        iz = 'Z';
        if (iz == 90 || iz == 122) {

/*           ASCII character set */

            if (ic >= 97 && ic <= 122) {
                *(unsigned char *)subnam = (char) (ic - 32);
                for (i__ = 2; i__ <= 6; ++i__) {
                    ic = *(unsigned char *)&subnam[i__ - 1];
                    if (ic >= 97 && ic <= 122) {
                        *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
                    }
                }
            }

        } else if (iz == 233 || iz == 169) {

/*           EBCDIC character set */

            if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162
                    && ic <= 169) {
                *(unsigned char *)subnam = (char) (ic + 64);
                for (i__ = 2; i__ <= 6; ++i__) {
                    ic = *(unsigned char *)&subnam[i__ - 1];
                    if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 ||
                            ic >= 162 && ic <= 169) {
                        *(unsigned char *)&subnam[i__ - 1] = (char) (ic + 64);
                    }
                }
            }

        } else if (iz == 218 || iz == 250) {

/*           Prime machines:  ASCII+128 */

            if (ic >= 225 && ic <= 250) {
                *(unsigned char *)subnam = (char) (ic - 32);
                for (i__ = 2; i__ <= 6; ++i__) {
                    ic = *(unsigned char *)&subnam[i__ - 1];
                    if (ic >= 225 && ic <= 250) {
                        *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
                    }
                }
            }
        }

        if (s_lmp_cmp(subnam + 1, (char *)"GGHRD", (ftnlen)5, (ftnlen)5) == 0 || s_lmp_cmp(
                subnam + 1, (char *)"GGHD3", (ftnlen)5, (ftnlen)5) == 0) {
            ret_val = 1;
            if (nh >= 14) {
                ret_val = 2;
            }
        } else if (s_lmp_cmp(subnam + 3, (char *)"EXC", (ftnlen)3, (ftnlen)3) == 0) {
            if (nh >= 14) {
                ret_val = 1;
            }
            if (nh >= 14) {
                ret_val = 2;
            }
        } else if (s_lmp_cmp(subnam + 1, (char *)"HSEQR", (ftnlen)5, (ftnlen)5) == 0 ||
                s_lmp_cmp(subnam + 1, (char *)"LAQR", (ftnlen)4, (ftnlen)4) == 0) {
            if (ns >= 14) {
                ret_val = 1;
            }
            if (ns >= 14) {
                ret_val = 2;
            }
        }

    } else if (*ispec == 17) {

/*        === Relative cost of near-the-diagonal chase vs */
/*            BLAS updates === */

        ret_val = 10;
    } else {
/*        ===== invalid value of ispec ===== */
        ret_val = -1;

    }

/*     ==== End of IPARMQ ==== */

    return ret_val;
} /* iparmq_ */

#ifdef __cplusplus
        }
#endif
