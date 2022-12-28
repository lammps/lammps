/* fortran/dlasd4.f -- translated by f2c (version 20200916).
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

/* > \brief \b DLASD4 computes the square root of the i-th updated eigenvalue of a positive symmetric rank-one
 modification to a positive diagonal matrix. Used by dbdsdc. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASD4 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasd4.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasd4.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasd4.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASD4( N, I, D, Z, DELTA, RHO, SIGMA, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            I, INFO, N */
/*       DOUBLE PRECISION   RHO, SIGMA */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), DELTA( * ), WORK( * ), Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > This subroutine computes the square root of the I-th updated */
/* > eigenvalue of a positive symmetric rank-one modification to */
/* > a positive diagonal matrix whose entries are given as the squares */
/* > of the corresponding entries in the array d, and that */
/* > */
/* >        0 <= D(i) < D(j)  for  i < j */
/* > */
/* > and that RHO > 0. This is arranged by the calling routine, and is */
/* > no loss in generality.  The rank-one modified system is thus */
/* > */
/* >        diag( D ) * diag( D ) +  RHO * Z * Z_transpose. */
/* > */
/* > where we assume the Euclidean norm of Z is 1. */
/* > */
/* > The method consists of approximating the rational functions in the */
/* > secular equation by simpler interpolating rational functions. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         The length of all arrays. */
/* > \endverbatim */
/* > */
/* > \param[in] I */
/* > \verbatim */
/* >          I is INTEGER */
/* >         The index of the eigenvalue to be computed.  1 <= I <= N. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension ( N ) */
/* >         The original eigenvalues.  It is assumed that they are in */
/* >         order, 0 <= D(I) < D(J)  for I < J. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension ( N ) */
/* >         The components of the updating vector. */
/* > \endverbatim */
/* > */
/* > \param[out] DELTA */
/* > \verbatim */
/* >          DELTA is DOUBLE PRECISION array, dimension ( N ) */
/* >         If N .ne. 1, DELTA contains (D(j) - sigma_I) in its  j-th */
/* >         component.  If N = 1, then DELTA(1) = 1.  The vector DELTA */
/* >         contains the information necessary to construct the */
/* >         (singular) eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] RHO */
/* > \verbatim */
/* >          RHO is DOUBLE PRECISION */
/* >         The scalar in the symmetric updating formula. */
/* > \endverbatim */
/* > */
/* > \param[out] SIGMA */
/* > \verbatim */
/* >          SIGMA is DOUBLE PRECISION */
/* >         The computed sigma_I, the I-th updated eigenvalue. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension ( N ) */
/* >         If N .ne. 1, WORK contains (D(j) + sigma_I) in its  j-th */
/* >         component.  If N = 1, then WORK( 1 ) = 1. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >         = 0:  successful exit */
/* >         > 0:  if INFO = 1, the updating process failed. */
/* > \endverbatim */

/* > \par Internal Parameters: */
/*  ========================= */
/* > */
/* > \verbatim */
/* >  Logical variable ORGATI (origin-at-i?) is used for distinguishing */
/* >  whether D(i) or D(i+1) is treated as the origin. */
/* > */
/* >            ORGATI = .true.    origin at i */
/* >            ORGATI = .false.   origin at i+1 */
/* > */
/* >  Logical variable SWTCH3 (switch-for-3-poles?) is for noting */
/* >  if we are working with THREE poles! */
/* > */
/* >  MAXIT is the maximum number of iterations allowed for each */
/* >  eigenvalue. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup OTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ren-Cang Li, Computer Science Division, University of California */
/* >     at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlasd4_(integer *n, integer *i__, doublereal *d__,
        doublereal *z__, doublereal *delta, doublereal *rho, doublereal *
        sigma, doublereal *work, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal a, b, c__;
    integer j;
    doublereal w, dd[3];
    integer ii;
    doublereal dw, zz[3];
    integer ip1;
    doublereal sq2, eta, phi, eps, tau, psi;
    integer iim1, iip1;
    doublereal tau2, dphi, sglb, dpsi, sgub;
    integer iter;
    doublereal temp, prew, temp1, temp2, dtiim, delsq, dtiip;
    integer niter;
    doublereal dtisq;
    logical swtch;
    doublereal dtnsq;
    extern /* Subroutine */ int dlaed6_(integer *, logical *, doublereal *,
            doublereal *, doublereal *, doublereal *, doublereal *, integer *)
            , dlasd5_(integer *, doublereal *, doublereal *, doublereal *,
            doublereal *, doublereal *, doublereal *);
    doublereal delsq2, dtnsq1;
    logical swtch3;
    extern doublereal dlamch_(char *, ftnlen);
    logical orgati;
    doublereal erretm, dtipsq, rhoinv;
    logical geomavg;


/*  -- LAPACK auxiliary routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Since this routine is called in an inner loop, we do no argument */
/*     checking. */

/*     Quick return for N=1 and 2. */

    /* Parameter adjustments */
    --work;
    --delta;
    --z__;
    --d__;

    /* Function Body */
    *info = 0;
    if (*n == 1) {

/*        Presumably, I=1 upon entry */

        *sigma = sqrt(d__[1] * d__[1] + *rho * z__[1] * z__[1]);
        delta[1] = 1.;
        work[1] = 1.;
        return 0;
    }
    if (*n == 2) {
        dlasd5_(i__, &d__[1], &z__[1], &delta[1], rho, sigma, &work[1]);
        return 0;
    }

/*     Compute machine epsilon */

    eps = dlamch_((char *)"Epsilon", (ftnlen)7);
    rhoinv = 1. / *rho;
    tau2 = 0.;

/*     The case I = N */

    if (*i__ == *n) {

/*        Initialize some basic variables */

        ii = *n - 1;
        niter = 1;

/*        Calculate initial guess */

        temp = *rho / 2.;

/*        If ||Z||_2 is not one, then TEMP should be set to */
/*        RHO * ||Z||_2^2 / TWO */

        temp1 = temp / (d__[*n] + sqrt(d__[*n] * d__[*n] + temp));
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            work[j] = d__[j] + d__[*n] + temp1;
            delta[j] = d__[j] - d__[*n] - temp1;
/* L10: */
        }

        psi = 0.;
        i__1 = *n - 2;
        for (j = 1; j <= i__1; ++j) {
            psi += z__[j] * z__[j] / (delta[j] * work[j]);
/* L20: */
        }

        c__ = rhoinv + psi;
        w = c__ + z__[ii] * z__[ii] / (delta[ii] * work[ii]) + z__[*n] * z__[*
                n] / (delta[*n] * work[*n]);

        if (w <= 0.) {
            temp1 = sqrt(d__[*n] * d__[*n] + *rho);
            temp = z__[*n - 1] * z__[*n - 1] / ((d__[*n - 1] + temp1) * (d__[*
                    n] - d__[*n - 1] + *rho / (d__[*n] + temp1))) + z__[*n] *
                    z__[*n] / *rho;

/*           The following TAU2 is to approximate */
/*           SIGMA_n^2 - D( N )*D( N ) */

            if (c__ <= temp) {
                tau = *rho;
            } else {
                delsq = (d__[*n] - d__[*n - 1]) * (d__[*n] + d__[*n - 1]);
                a = -c__ * delsq + z__[*n - 1] * z__[*n - 1] + z__[*n] * z__[*
                        n];
                b = z__[*n] * z__[*n] * delsq;
                if (a < 0.) {
                    tau2 = b * 2. / (sqrt(a * a + b * 4. * c__) - a);
                } else {
                    tau2 = (a + sqrt(a * a + b * 4. * c__)) / (c__ * 2.);
                }
                tau = tau2 / (d__[*n] + sqrt(d__[*n] * d__[*n] + tau2));
            }

/*           It can be proved that */
/*               D(N)^2+RHO/2 <= SIGMA_n^2 < D(N)^2+TAU2 <= D(N)^2+RHO */

        } else {
            delsq = (d__[*n] - d__[*n - 1]) * (d__[*n] + d__[*n - 1]);
            a = -c__ * delsq + z__[*n - 1] * z__[*n - 1] + z__[*n] * z__[*n];
            b = z__[*n] * z__[*n] * delsq;

/*           The following TAU2 is to approximate */
/*           SIGMA_n^2 - D( N )*D( N ) */

            if (a < 0.) {
                tau2 = b * 2. / (sqrt(a * a + b * 4. * c__) - a);
            } else {
                tau2 = (a + sqrt(a * a + b * 4. * c__)) / (c__ * 2.);
            }
            tau = tau2 / (d__[*n] + sqrt(d__[*n] * d__[*n] + tau2));

/*           It can be proved that */
/*           D(N)^2 < D(N)^2+TAU2 < SIGMA(N)^2 < D(N)^2+RHO/2 */

        }

/*        The following TAU is to approximate SIGMA_n - D( N ) */

/*         TAU = TAU2 / ( D( N )+SQRT( D( N )*D( N )+TAU2 ) ) */

        *sigma = d__[*n] + tau;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            delta[j] = d__[j] - d__[*n] - tau;
            work[j] = d__[j] + d__[*n] + tau;
/* L30: */
        }

/*        Evaluate PSI and the derivative DPSI */

        dpsi = 0.;
        psi = 0.;
        erretm = 0.;
        i__1 = ii;
        for (j = 1; j <= i__1; ++j) {
            temp = z__[j] / (delta[j] * work[j]);
            psi += z__[j] * temp;
            dpsi += temp * temp;
            erretm += psi;
/* L40: */
        }
        erretm = abs(erretm);

/*        Evaluate PHI and the derivative DPHI */

        temp = z__[*n] / (delta[*n] * work[*n]);
        phi = z__[*n] * temp;
        dphi = temp * temp;
        erretm = (-phi - psi) * 8. + erretm - phi + rhoinv;
/*    $          + ABS( TAU2 )*( DPSI+DPHI ) */

        w = rhoinv + phi + psi;

/*        Test for convergence */

        if (abs(w) <= eps * erretm) {
            goto L240;
        }

/*        Calculate the new step */

        ++niter;
        dtnsq1 = work[*n - 1] * delta[*n - 1];
        dtnsq = work[*n] * delta[*n];
        c__ = w - dtnsq1 * dpsi - dtnsq * dphi;
        a = (dtnsq + dtnsq1) * w - dtnsq * dtnsq1 * (dpsi + dphi);
        b = dtnsq * dtnsq1 * w;
        if (c__ < 0.) {
            c__ = abs(c__);
        }
        if (c__ == 0.) {
            eta = *rho - *sigma * *sigma;
        } else if (a >= 0.) {
            eta = (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (c__
                    * 2.);
        } else {
            eta = b * 2. / (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))
                    );
        }

/*        Note, eta should be positive if w is negative, and */
/*        eta should be negative otherwise. However, */
/*        if for some reason caused by roundoff, eta*w > 0, */
/*        we simply use one Newton step instead. This way */
/*        will guarantee eta*w < 0. */

        if (w * eta > 0.) {
            eta = -w / (dpsi + dphi);
        }
        temp = eta - dtnsq;
        if (temp > *rho) {
            eta = *rho + dtnsq;
        }

        eta /= *sigma + sqrt(eta + *sigma * *sigma);
        tau += eta;
        *sigma += eta;

        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            delta[j] -= eta;
            work[j] += eta;
/* L50: */
        }

/*        Evaluate PSI and the derivative DPSI */

        dpsi = 0.;
        psi = 0.;
        erretm = 0.;
        i__1 = ii;
        for (j = 1; j <= i__1; ++j) {
            temp = z__[j] / (work[j] * delta[j]);
            psi += z__[j] * temp;
            dpsi += temp * temp;
            erretm += psi;
/* L60: */
        }
        erretm = abs(erretm);

/*        Evaluate PHI and the derivative DPHI */

        tau2 = work[*n] * delta[*n];
        temp = z__[*n] / tau2;
        phi = z__[*n] * temp;
        dphi = temp * temp;
        erretm = (-phi - psi) * 8. + erretm - phi + rhoinv;
/*    $          + ABS( TAU2 )*( DPSI+DPHI ) */

        w = rhoinv + phi + psi;

/*        Main loop to update the values of the array   DELTA */

        iter = niter + 1;

        for (niter = iter; niter <= 400; ++niter) {

/*           Test for convergence */

            if (abs(w) <= eps * erretm) {
                goto L240;
            }

/*           Calculate the new step */

            dtnsq1 = work[*n - 1] * delta[*n - 1];
            dtnsq = work[*n] * delta[*n];
            c__ = w - dtnsq1 * dpsi - dtnsq * dphi;
            a = (dtnsq + dtnsq1) * w - dtnsq1 * dtnsq * (dpsi + dphi);
            b = dtnsq1 * dtnsq * w;
            if (a >= 0.) {
                eta = (a + sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (
                        c__ * 2.);
            } else {
                eta = b * 2. / (a - sqrt((d__1 = a * a - b * 4. * c__, abs(
                        d__1))));
            }

/*           Note, eta should be positive if w is negative, and */
/*           eta should be negative otherwise. However, */
/*           if for some reason caused by roundoff, eta*w > 0, */
/*           we simply use one Newton step instead. This way */
/*           will guarantee eta*w < 0. */

            if (w * eta > 0.) {
                eta = -w / (dpsi + dphi);
            }
            temp = eta - dtnsq;
            if (temp <= 0.) {
                eta /= 2.;
            }

            eta /= *sigma + sqrt(eta + *sigma * *sigma);
            tau += eta;
            *sigma += eta;

            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                delta[j] -= eta;
                work[j] += eta;
/* L70: */
            }

/*           Evaluate PSI and the derivative DPSI */

            dpsi = 0.;
            psi = 0.;
            erretm = 0.;
            i__1 = ii;
            for (j = 1; j <= i__1; ++j) {
                temp = z__[j] / (work[j] * delta[j]);
                psi += z__[j] * temp;
                dpsi += temp * temp;
                erretm += psi;
/* L80: */
            }
            erretm = abs(erretm);

/*           Evaluate PHI and the derivative DPHI */

            tau2 = work[*n] * delta[*n];
            temp = z__[*n] / tau2;
            phi = z__[*n] * temp;
            dphi = temp * temp;
            erretm = (-phi - psi) * 8. + erretm - phi + rhoinv;
/*    $             + ABS( TAU2 )*( DPSI+DPHI ) */

            w = rhoinv + phi + psi;
/* L90: */
        }

/*        Return with INFO = 1, NITER = MAXIT and not converged */

        *info = 1;
        goto L240;

/*        End for the case I = N */

    } else {

/*        The case for I < N */

        niter = 1;
        ip1 = *i__ + 1;

/*        Calculate initial guess */

        delsq = (d__[ip1] - d__[*i__]) * (d__[ip1] + d__[*i__]);
        delsq2 = delsq / 2.;
        sq2 = sqrt((d__[*i__] * d__[*i__] + d__[ip1] * d__[ip1]) / 2.);
        temp = delsq2 / (d__[*i__] + sq2);
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            work[j] = d__[j] + d__[*i__] + temp;
            delta[j] = d__[j] - d__[*i__] - temp;
/* L100: */
        }

        psi = 0.;
        i__1 = *i__ - 1;
        for (j = 1; j <= i__1; ++j) {
            psi += z__[j] * z__[j] / (work[j] * delta[j]);
/* L110: */
        }

        phi = 0.;
        i__1 = *i__ + 2;
        for (j = *n; j >= i__1; --j) {
            phi += z__[j] * z__[j] / (work[j] * delta[j]);
/* L120: */
        }
        c__ = rhoinv + psi + phi;
        w = c__ + z__[*i__] * z__[*i__] / (work[*i__] * delta[*i__]) + z__[
                ip1] * z__[ip1] / (work[ip1] * delta[ip1]);

        geomavg = FALSE_;
        if (w > 0.) {

/*           d(i)^2 < the ith sigma^2 < (d(i)^2+d(i+1)^2)/2 */

/*           We choose d(i) as origin. */

            orgati = TRUE_;
            ii = *i__;
            sglb = 0.;
            sgub = delsq2 / (d__[*i__] + sq2);
            a = c__ * delsq + z__[*i__] * z__[*i__] + z__[ip1] * z__[ip1];
            b = z__[*i__] * z__[*i__] * delsq;
            if (a > 0.) {
                tau2 = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(
                        d__1))));
            } else {
                tau2 = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) /
                        (c__ * 2.);
            }

/*           TAU2 now is an estimation of SIGMA^2 - D( I )^2. The */
/*           following, however, is the corresponding estimation of */
/*           SIGMA - D( I ). */

            tau = tau2 / (d__[*i__] + sqrt(d__[*i__] * d__[*i__] + tau2));
            temp = sqrt(eps);
            if (d__[*i__] <= temp * d__[ip1] && (d__1 = z__[*i__], abs(d__1))
                    <= temp && d__[*i__] > 0.) {
/* Computing MIN */
                d__1 = d__[*i__] * 10.;
                tau = min(d__1,sgub);
                geomavg = TRUE_;
            }
        } else {

/*           (d(i)^2+d(i+1)^2)/2 <= the ith sigma^2 < d(i+1)^2/2 */

/*           We choose d(i+1) as origin. */

            orgati = FALSE_;
            ii = ip1;
            sglb = -delsq2 / (d__[ii] + sq2);
            sgub = 0.;
            a = c__ * delsq - z__[*i__] * z__[*i__] - z__[ip1] * z__[ip1];
            b = z__[ip1] * z__[ip1] * delsq;
            if (a < 0.) {
                tau2 = b * 2. / (a - sqrt((d__1 = a * a + b * 4. * c__, abs(
                        d__1))));
            } else {
                tau2 = -(a + sqrt((d__1 = a * a + b * 4. * c__, abs(d__1)))) /
                         (c__ * 2.);
            }

/*           TAU2 now is an estimation of SIGMA^2 - D( IP1 )^2. The */
/*           following, however, is the corresponding estimation of */
/*           SIGMA - D( IP1 ). */

            tau = tau2 / (d__[ip1] + sqrt((d__1 = d__[ip1] * d__[ip1] + tau2,
                    abs(d__1))));
        }

        *sigma = d__[ii] + tau;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            work[j] = d__[j] + d__[ii] + tau;
            delta[j] = d__[j] - d__[ii] - tau;
/* L130: */
        }
        iim1 = ii - 1;
        iip1 = ii + 1;

/*        Evaluate PSI and the derivative DPSI */

        dpsi = 0.;
        psi = 0.;
        erretm = 0.;
        i__1 = iim1;
        for (j = 1; j <= i__1; ++j) {
            temp = z__[j] / (work[j] * delta[j]);
            psi += z__[j] * temp;
            dpsi += temp * temp;
            erretm += psi;
/* L150: */
        }
        erretm = abs(erretm);

/*        Evaluate PHI and the derivative DPHI */

        dphi = 0.;
        phi = 0.;
        i__1 = iip1;
        for (j = *n; j >= i__1; --j) {
            temp = z__[j] / (work[j] * delta[j]);
            phi += z__[j] * temp;
            dphi += temp * temp;
            erretm += phi;
/* L160: */
        }

        w = rhoinv + phi + psi;

/*        W is the value of the secular function with */
/*        its ii-th element removed. */

        swtch3 = FALSE_;
        if (orgati) {
            if (w < 0.) {
                swtch3 = TRUE_;
            }
        } else {
            if (w > 0.) {
                swtch3 = TRUE_;
            }
        }
        if (ii == 1 || ii == *n) {
            swtch3 = FALSE_;
        }

        temp = z__[ii] / (work[ii] * delta[ii]);
        dw = dpsi + dphi + temp * temp;
        temp = z__[ii] * temp;
        w += temp;
        erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + abs(temp) * 3.;
/*    $          + ABS( TAU2 )*DW */

/*        Test for convergence */

        if (abs(w) <= eps * erretm) {
            goto L240;
        }

        if (w <= 0.) {
            sglb = max(sglb,tau);
        } else {
            sgub = min(sgub,tau);
        }

/*        Calculate the new step */

        ++niter;
        if (! swtch3) {
            dtipsq = work[ip1] * delta[ip1];
            dtisq = work[*i__] * delta[*i__];
            if (orgati) {
/* Computing 2nd power */
                d__1 = z__[*i__] / dtisq;
                c__ = w - dtipsq * dw + delsq * (d__1 * d__1);
            } else {
/* Computing 2nd power */
                d__1 = z__[ip1] / dtipsq;
                c__ = w - dtisq * dw - delsq * (d__1 * d__1);
            }
            a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
            b = dtipsq * dtisq * w;
            if (c__ == 0.) {
                if (a == 0.) {
                    if (orgati) {
                        a = z__[*i__] * z__[*i__] + dtipsq * dtipsq * (dpsi +
                                dphi);
                    } else {
                        a = z__[ip1] * z__[ip1] + dtisq * dtisq * (dpsi +
                                dphi);
                    }
                }
                eta = b / a;
            } else if (a <= 0.) {
                eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1)))) / (
                        c__ * 2.);
            } else {
                eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__, abs(
                        d__1))));
            }
        } else {

/*           Interpolation using THREE most relevant poles */

            dtiim = work[iim1] * delta[iim1];
            dtiip = work[iip1] * delta[iip1];
            temp = rhoinv + psi + phi;
            if (orgati) {
                temp1 = z__[iim1] / dtiim;
                temp1 *= temp1;
                c__ = temp - dtiip * (dpsi + dphi) - (d__[iim1] - d__[iip1]) *
                         (d__[iim1] + d__[iip1]) * temp1;
                zz[0] = z__[iim1] * z__[iim1];
                if (dpsi < temp1) {
                    zz[2] = dtiip * dtiip * dphi;
                } else {
                    zz[2] = dtiip * dtiip * (dpsi - temp1 + dphi);
                }
            } else {
                temp1 = z__[iip1] / dtiip;
                temp1 *= temp1;
                c__ = temp - dtiim * (dpsi + dphi) - (d__[iip1] - d__[iim1]) *
                         (d__[iim1] + d__[iip1]) * temp1;
                if (dphi < temp1) {
                    zz[0] = dtiim * dtiim * dpsi;
                } else {
                    zz[0] = dtiim * dtiim * (dpsi + (dphi - temp1));
                }
                zz[2] = z__[iip1] * z__[iip1];
            }
            zz[1] = z__[ii] * z__[ii];
            dd[0] = dtiim;
            dd[1] = delta[ii] * work[ii];
            dd[2] = dtiip;
            dlaed6_(&niter, &orgati, &c__, dd, zz, &w, &eta, info);

            if (*info != 0) {

/*              If INFO is not 0, i.e., DLAED6 failed, switch back */
/*              to 2 pole interpolation. */

                swtch3 = FALSE_;
                *info = 0;
                dtipsq = work[ip1] * delta[ip1];
                dtisq = work[*i__] * delta[*i__];
                if (orgati) {
/* Computing 2nd power */
                    d__1 = z__[*i__] / dtisq;
                    c__ = w - dtipsq * dw + delsq * (d__1 * d__1);
                } else {
/* Computing 2nd power */
                    d__1 = z__[ip1] / dtipsq;
                    c__ = w - dtisq * dw - delsq * (d__1 * d__1);
                }
                a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
                b = dtipsq * dtisq * w;
                if (c__ == 0.) {
                    if (a == 0.) {
                        if (orgati) {
                            a = z__[*i__] * z__[*i__] + dtipsq * dtipsq * (
                                    dpsi + dphi);
                        } else {
                            a = z__[ip1] * z__[ip1] + dtisq * dtisq * (dpsi +
                                    dphi);
                        }
                    }
                    eta = b / a;
                } else if (a <= 0.) {
                    eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))))
                             / (c__ * 2.);
                } else {
                    eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__,
                            abs(d__1))));
                }
            }
        }

/*        Note, eta should be positive if w is negative, and */
/*        eta should be negative otherwise. However, */
/*        if for some reason caused by roundoff, eta*w > 0, */
/*        we simply use one Newton step instead. This way */
/*        will guarantee eta*w < 0. */

        if (w * eta >= 0.) {
            eta = -w / dw;
        }

        eta /= *sigma + sqrt(*sigma * *sigma + eta);
        temp = tau + eta;
        if (temp > sgub || temp < sglb) {
            if (w < 0.) {
                eta = (sgub - tau) / 2.;
            } else {
                eta = (sglb - tau) / 2.;
            }
            if (geomavg) {
                if (w < 0.) {
                    if (tau > 0.) {
                        eta = sqrt(sgub * tau) - tau;
                    }
                } else {
                    if (sglb > 0.) {
                        eta = sqrt(sglb * tau) - tau;
                    }
                }
            }
        }

        prew = w;

        tau += eta;
        *sigma += eta;

        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            work[j] += eta;
            delta[j] -= eta;
/* L170: */
        }

/*        Evaluate PSI and the derivative DPSI */

        dpsi = 0.;
        psi = 0.;
        erretm = 0.;
        i__1 = iim1;
        for (j = 1; j <= i__1; ++j) {
            temp = z__[j] / (work[j] * delta[j]);
            psi += z__[j] * temp;
            dpsi += temp * temp;
            erretm += psi;
/* L180: */
        }
        erretm = abs(erretm);

/*        Evaluate PHI and the derivative DPHI */

        dphi = 0.;
        phi = 0.;
        i__1 = iip1;
        for (j = *n; j >= i__1; --j) {
            temp = z__[j] / (work[j] * delta[j]);
            phi += z__[j] * temp;
            dphi += temp * temp;
            erretm += phi;
/* L190: */
        }

        tau2 = work[ii] * delta[ii];
        temp = z__[ii] / tau2;
        dw = dpsi + dphi + temp * temp;
        temp = z__[ii] * temp;
        w = rhoinv + phi + psi + temp;
        erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + abs(temp) * 3.;
/*    $          + ABS( TAU2 )*DW */

        swtch = FALSE_;
        if (orgati) {
            if (-w > abs(prew) / 10.) {
                swtch = TRUE_;
            }
        } else {
            if (w > abs(prew) / 10.) {
                swtch = TRUE_;
            }
        }

/*        Main loop to update the values of the array   DELTA and WORK */

        iter = niter + 1;

        for (niter = iter; niter <= 400; ++niter) {

/*           Test for convergence */

            if (abs(w) <= eps * erretm) {
/*     $          .OR. (SGUB-SGLB).LE.EIGHT*ABS(SGUB+SGLB) ) THEN */
                goto L240;
            }

            if (w <= 0.) {
                sglb = max(sglb,tau);
            } else {
                sgub = min(sgub,tau);
            }

/*           Calculate the new step */

            if (! swtch3) {
                dtipsq = work[ip1] * delta[ip1];
                dtisq = work[*i__] * delta[*i__];
                if (! swtch) {
                    if (orgati) {
/* Computing 2nd power */
                        d__1 = z__[*i__] / dtisq;
                        c__ = w - dtipsq * dw + delsq * (d__1 * d__1);
                    } else {
/* Computing 2nd power */
                        d__1 = z__[ip1] / dtipsq;
                        c__ = w - dtisq * dw - delsq * (d__1 * d__1);
                    }
                } else {
                    temp = z__[ii] / (work[ii] * delta[ii]);
                    if (orgati) {
                        dpsi += temp * temp;
                    } else {
                        dphi += temp * temp;
                    }
                    c__ = w - dtisq * dpsi - dtipsq * dphi;
                }
                a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
                b = dtipsq * dtisq * w;
                if (c__ == 0.) {
                    if (a == 0.) {
                        if (! swtch) {
                            if (orgati) {
                                a = z__[*i__] * z__[*i__] + dtipsq * dtipsq *
                                        (dpsi + dphi);
                            } else {
                                a = z__[ip1] * z__[ip1] + dtisq * dtisq * (
                                        dpsi + dphi);
                            }
                        } else {
                            a = dtisq * dtisq * dpsi + dtipsq * dtipsq * dphi;
                        }
                    }
                    eta = b / a;
                } else if (a <= 0.) {
                    eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(d__1))))
                             / (c__ * 2.);
                } else {
                    eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__,
                            abs(d__1))));
                }
            } else {

/*              Interpolation using THREE most relevant poles */

                dtiim = work[iim1] * delta[iim1];
                dtiip = work[iip1] * delta[iip1];
                temp = rhoinv + psi + phi;
                if (swtch) {
                    c__ = temp - dtiim * dpsi - dtiip * dphi;
                    zz[0] = dtiim * dtiim * dpsi;
                    zz[2] = dtiip * dtiip * dphi;
                } else {
                    if (orgati) {
                        temp1 = z__[iim1] / dtiim;
                        temp1 *= temp1;
                        temp2 = (d__[iim1] - d__[iip1]) * (d__[iim1] + d__[
                                iip1]) * temp1;
                        c__ = temp - dtiip * (dpsi + dphi) - temp2;
                        zz[0] = z__[iim1] * z__[iim1];
                        if (dpsi < temp1) {
                            zz[2] = dtiip * dtiip * dphi;
                        } else {
                            zz[2] = dtiip * dtiip * (dpsi - temp1 + dphi);
                        }
                    } else {
                        temp1 = z__[iip1] / dtiip;
                        temp1 *= temp1;
                        temp2 = (d__[iip1] - d__[iim1]) * (d__[iim1] + d__[
                                iip1]) * temp1;
                        c__ = temp - dtiim * (dpsi + dphi) - temp2;
                        if (dphi < temp1) {
                            zz[0] = dtiim * dtiim * dpsi;
                        } else {
                            zz[0] = dtiim * dtiim * (dpsi + (dphi - temp1));
                        }
                        zz[2] = z__[iip1] * z__[iip1];
                    }
                }
                dd[0] = dtiim;
                dd[1] = delta[ii] * work[ii];
                dd[2] = dtiip;
                dlaed6_(&niter, &orgati, &c__, dd, zz, &w, &eta, info);

                if (*info != 0) {

/*                 If INFO is not 0, i.e., DLAED6 failed, switch */
/*                 back to two pole interpolation */

                    swtch3 = FALSE_;
                    *info = 0;
                    dtipsq = work[ip1] * delta[ip1];
                    dtisq = work[*i__] * delta[*i__];
                    if (! swtch) {
                        if (orgati) {
/* Computing 2nd power */
                            d__1 = z__[*i__] / dtisq;
                            c__ = w - dtipsq * dw + delsq * (d__1 * d__1);
                        } else {
/* Computing 2nd power */
                            d__1 = z__[ip1] / dtipsq;
                            c__ = w - dtisq * dw - delsq * (d__1 * d__1);
                        }
                    } else {
                        temp = z__[ii] / (work[ii] * delta[ii]);
                        if (orgati) {
                            dpsi += temp * temp;
                        } else {
                            dphi += temp * temp;
                        }
                        c__ = w - dtisq * dpsi - dtipsq * dphi;
                    }
                    a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
                    b = dtipsq * dtisq * w;
                    if (c__ == 0.) {
                        if (a == 0.) {
                            if (! swtch) {
                                if (orgati) {
                                    a = z__[*i__] * z__[*i__] + dtipsq *
                                            dtipsq * (dpsi + dphi);
                                } else {
                                    a = z__[ip1] * z__[ip1] + dtisq * dtisq *
                                            (dpsi + dphi);
                                }
                            } else {
                                a = dtisq * dtisq * dpsi + dtipsq * dtipsq *
                                        dphi;
                            }
                        }
                        eta = b / a;
                    } else if (a <= 0.) {
                        eta = (a - sqrt((d__1 = a * a - b * 4. * c__, abs(
                                d__1)))) / (c__ * 2.);
                    } else {
                        eta = b * 2. / (a + sqrt((d__1 = a * a - b * 4. * c__,
                                 abs(d__1))));
                    }
                }
            }

/*           Note, eta should be positive if w is negative, and */
/*           eta should be negative otherwise. However, */
/*           if for some reason caused by roundoff, eta*w > 0, */
/*           we simply use one Newton step instead. This way */
/*           will guarantee eta*w < 0. */

            if (w * eta >= 0.) {
                eta = -w / dw;
            }

            eta /= *sigma + sqrt(*sigma * *sigma + eta);
            temp = tau + eta;
            if (temp > sgub || temp < sglb) {
                if (w < 0.) {
                    eta = (sgub - tau) / 2.;
                } else {
                    eta = (sglb - tau) / 2.;
                }
                if (geomavg) {
                    if (w < 0.) {
                        if (tau > 0.) {
                            eta = sqrt(sgub * tau) - tau;
                        }
                    } else {
                        if (sglb > 0.) {
                            eta = sqrt(sglb * tau) - tau;
                        }
                    }
                }
            }

            prew = w;

            tau += eta;
            *sigma += eta;

            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                work[j] += eta;
                delta[j] -= eta;
/* L200: */
            }

/*           Evaluate PSI and the derivative DPSI */

            dpsi = 0.;
            psi = 0.;
            erretm = 0.;
            i__1 = iim1;
            for (j = 1; j <= i__1; ++j) {
                temp = z__[j] / (work[j] * delta[j]);
                psi += z__[j] * temp;
                dpsi += temp * temp;
                erretm += psi;
/* L210: */
            }
            erretm = abs(erretm);

/*           Evaluate PHI and the derivative DPHI */

            dphi = 0.;
            phi = 0.;
            i__1 = iip1;
            for (j = *n; j >= i__1; --j) {
                temp = z__[j] / (work[j] * delta[j]);
                phi += z__[j] * temp;
                dphi += temp * temp;
                erretm += phi;
/* L220: */
            }

            tau2 = work[ii] * delta[ii];
            temp = z__[ii] / tau2;
            dw = dpsi + dphi + temp * temp;
            temp = z__[ii] * temp;
            w = rhoinv + phi + psi + temp;
            erretm = (phi - psi) * 8. + erretm + rhoinv * 2. + abs(temp) * 3.;
/*    $             + ABS( TAU2 )*DW */

            if (w * prew > 0. && abs(w) > abs(prew) / 10.) {
                swtch = ! swtch;
            }

/* L230: */
        }

/*        Return with INFO = 1, NITER = MAXIT and not converged */

        *info = 1;

    }

L240:
    return 0;

/*     End of DLASD4 */

} /* dlasd4_ */

#ifdef __cplusplus
        }
#endif
