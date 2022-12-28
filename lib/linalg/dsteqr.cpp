/* fortran/dsteqr.f -- translated by f2c (version 20200916).
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

/* Table of constant values */

static doublereal c_b9 = 0.;
static doublereal c_b10 = 1.;
static integer c__0 = 0;
static integer c__1 = 1;
static integer c__2 = 2;

/* > \brief \b DSTEQR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSTEQR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsteqr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsteqr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsteqr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          COMPZ */
/*       INTEGER            INFO, LDZ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSTEQR computes all eigenvalues and, optionally, eigenvectors of a */
/* > symmetric tridiagonal matrix using the implicit QL or QR method. */
/* > The eigenvectors of a full or band symmetric matrix can also be found */
/* > if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to */
/* > tridiagonal form. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] COMPZ */
/* > \verbatim */
/* >          COMPZ is CHARACTER*1 */
/* >          = 'N':  Compute eigenvalues only. */
/* >          = 'V':  Compute eigenvalues and eigenvectors of the original */
/* >                  symmetric matrix.  On entry, Z must contain the */
/* >                  orthogonal matrix used to reduce the original matrix */
/* >                  to tridiagonal form. */
/* >          = 'I':  Compute eigenvalues and eigenvectors of the */
/* >                  tridiagonal matrix.  Z is initialized to the identity */
/* >                  matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the diagonal elements of the tridiagonal matrix. */
/* >          On exit, if INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >          On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* >          matrix. */
/* >          On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (LDZ, N) */
/* >          On entry, if  COMPZ = 'V', then Z contains the orthogonal */
/* >          matrix used in the reduction to tridiagonal form. */
/* >          On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the */
/* >          orthonormal eigenvectors of the original symmetric matrix, */
/* >          and if COMPZ = 'I', Z contains the orthonormal eigenvectors */
/* >          of the symmetric tridiagonal matrix. */
/* >          If COMPZ = 'N', then Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z.  LDZ >= 1, and if */
/* >          eigenvectors are desired, then  LDZ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (max(1,2*N-2)) */
/* >          If COMPZ = 'N', then WORK is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  the algorithm has failed to find all the eigenvalues in */
/* >                a total of 30*N iterations; if INFO = i, then i */
/* >                elements of E have not converged to zero; on exit, D */
/* >                and E contain the elements of a symmetric tridiagonal */
/* >                matrix which is orthogonally similar to the original */
/* >                matrix. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dsteqr_(char *compz, integer *n, doublereal *d__,
        doublereal *e, doublereal *z__, integer *ldz, doublereal *work,
        integer *info, ftnlen compz_len)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    doublereal b, c__, f, g;
    integer i__, j, k, l, m;
    doublereal p, r__, s;
    integer l1, ii, mm, lm1, mm1, nm1;
    doublereal rt1, rt2, eps;
    integer lsv;
    doublereal tst, eps2;
    integer lend, jtot;
    extern /* Subroutine */ int dlae2_(doublereal *, doublereal *, doublereal
            *, doublereal *, doublereal *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlasr_(char *, char *, char *, integer *,
            integer *, doublereal *, doublereal *, doublereal *, integer *,
            ftnlen, ftnlen, ftnlen);
    doublereal anorm;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *,
            doublereal *, integer *), dlaev2_(doublereal *, doublereal *,
            doublereal *, doublereal *, doublereal *, doublereal *,
            doublereal *);
    integer lendm1, lendp1;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *,
            ftnlen);
    integer iscale;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *,
            doublereal *, doublereal *, integer *, integer *, doublereal *,
            integer *, integer *, ftnlen), dlaset_(char *, integer *, integer
            *, doublereal *, doublereal *, doublereal *, integer *, ftnlen);
    doublereal safmin;
    extern /* Subroutine */ int dlartg_(doublereal *, doublereal *,
            doublereal *, doublereal *, doublereal *);
    doublereal safmax;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *,
            ftnlen);
    extern /* Subroutine */ int dlasrt_(char *, integer *, doublereal *,
            integer *, ftnlen);
    integer lendsv;
    doublereal ssfmin;
    integer nmaxit, icompz;
    doublereal ssfmax;


/*  -- LAPACK computational routine -- */
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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --d__;
    --e;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;

    /* Function Body */
    *info = 0;

    if (lsame_(compz, (char *)"N", (ftnlen)1, (ftnlen)1)) {
        icompz = 0;
    } else if (lsame_(compz, (char *)"V", (ftnlen)1, (ftnlen)1)) {
        icompz = 1;
    } else if (lsame_(compz, (char *)"I", (ftnlen)1, (ftnlen)1)) {
        icompz = 2;
    } else {
        icompz = -1;
    }
    if (icompz < 0) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*ldz < 1 || icompz > 0 && *ldz < max(1,*n)) {
        *info = -6;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DSTEQR", &i__1, (ftnlen)6);
        return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
        return 0;
    }

    if (*n == 1) {
        if (icompz == 2) {
            z__[z_dim1 + 1] = 1.;
        }
        return 0;
    }

/*     Determine the unit roundoff and over/underflow thresholds. */

    eps = dlamch_((char *)"E", (ftnlen)1);
/* Computing 2nd power */
    d__1 = eps;
    eps2 = d__1 * d__1;
    safmin = dlamch_((char *)"S", (ftnlen)1);
    safmax = 1. / safmin;
    ssfmax = sqrt(safmax) / 3.;
    ssfmin = sqrt(safmin) / eps2;

/*     Compute the eigenvalues and eigenvectors of the tridiagonal */
/*     matrix. */

    if (icompz == 2) {
        dlaset_((char *)"Full", n, n, &c_b9, &c_b10, &z__[z_offset], ldz, (ftnlen)4);
    }

    nmaxit = *n * 30;
    jtot = 0;

/*     Determine where the matrix splits and choose QL or QR iteration */
/*     for each block, according to whether top or bottom diagonal */
/*     element is smaller. */

    l1 = 1;
    nm1 = *n - 1;

L10:
    if (l1 > *n) {
        goto L160;
    }
    if (l1 > 1) {
        e[l1 - 1] = 0.;
    }
    if (l1 <= nm1) {
        i__1 = nm1;
        for (m = l1; m <= i__1; ++m) {
            tst = (d__1 = e[m], abs(d__1));
            if (tst == 0.) {
                goto L30;
            }
            if (tst <= sqrt((d__1 = d__[m], abs(d__1))) * sqrt((d__2 = d__[m
                    + 1], abs(d__2))) * eps) {
                e[m] = 0.;
                goto L30;
            }
/* L20: */
        }
    }
    m = *n;

L30:
    l = l1;
    lsv = l;
    lend = m;
    lendsv = lend;
    l1 = m + 1;
    if (lend == l) {
        goto L10;
    }

/*     Scale submatrix in rows and columns L to LEND */

    i__1 = lend - l + 1;
    anorm = dlanst_((char *)"M", &i__1, &d__[l], &e[l], (ftnlen)1);
    iscale = 0;
    if (anorm == 0.) {
        goto L10;
    }
    if (anorm > ssfmax) {
        iscale = 1;
        i__1 = lend - l + 1;
        dlascl_((char *)"G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l], n,
                info, (ftnlen)1);
        i__1 = lend - l;
        dlascl_((char *)"G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n,
                info, (ftnlen)1);
    } else if (anorm < ssfmin) {
        iscale = 2;
        i__1 = lend - l + 1;
        dlascl_((char *)"G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l], n,
                info, (ftnlen)1);
        i__1 = lend - l;
        dlascl_((char *)"G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n,
                info, (ftnlen)1);
    }

/*     Choose between QL and QR iteration */

    if ((d__1 = d__[lend], abs(d__1)) < (d__2 = d__[l], abs(d__2))) {
        lend = lsv;
        l = lendsv;
    }

    if (lend > l) {

/*        QL Iteration */

/*        Look for small subdiagonal element. */

L40:
        if (l != lend) {
            lendm1 = lend - 1;
            i__1 = lendm1;
            for (m = l; m <= i__1; ++m) {
/* Computing 2nd power */
                d__2 = (d__1 = e[m], abs(d__1));
                tst = d__2 * d__2;
                if (tst <= eps2 * (d__1 = d__[m], abs(d__1)) * (d__2 = d__[m
                        + 1], abs(d__2)) + safmin) {
                    goto L60;
                }
/* L50: */
            }
        }

        m = lend;

L60:
        if (m < lend) {
            e[m] = 0.;
        }
        p = d__[l];
        if (m == l) {
            goto L80;
        }

/*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2 */
/*        to compute its eigensystem. */

        if (m == l + 1) {
            if (icompz > 0) {
                dlaev2_(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2, &c__, &s);
                work[l] = c__;
                work[*n - 1 + l] = s;
                dlasr_((char *)"R", (char *)"V", (char *)"B", n, &c__2, &work[l], &work[*n - 1 + l], &
                        z__[l * z_dim1 + 1], ldz, (ftnlen)1, (ftnlen)1, (
                        ftnlen)1);
            } else {
                dlae2_(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2);
            }
            d__[l] = rt1;
            d__[l + 1] = rt2;
            e[l] = 0.;
            l += 2;
            if (l <= lend) {
                goto L40;
            }
            goto L140;
        }

        if (jtot == nmaxit) {
            goto L140;
        }
        ++jtot;

/*        Form shift. */

        g = (d__[l + 1] - p) / (e[l] * 2.);
        r__ = dlapy2_(&g, &c_b10);
        g = d__[m] - p + e[l] / (g + d_sign(&r__, &g));

        s = 1.;
        c__ = 1.;
        p = 0.;

/*        Inner loop */

        mm1 = m - 1;
        i__1 = l;
        for (i__ = mm1; i__ >= i__1; --i__) {
            f = s * e[i__];
            b = c__ * e[i__];
            dlartg_(&g, &f, &c__, &s, &r__);
            if (i__ != m - 1) {
                e[i__ + 1] = r__;
            }
            g = d__[i__ + 1] - p;
            r__ = (d__[i__] - g) * s + c__ * 2. * b;
            p = s * r__;
            d__[i__ + 1] = g + p;
            g = c__ * r__ - b;

/*           If eigenvectors are desired, then save rotations. */

            if (icompz > 0) {
                work[i__] = c__;
                work[*n - 1 + i__] = -s;
            }

/* L70: */
        }

/*        If eigenvectors are desired, then apply saved rotations. */

        if (icompz > 0) {
            mm = m - l + 1;
            dlasr_((char *)"R", (char *)"V", (char *)"B", n, &mm, &work[l], &work[*n - 1 + l], &z__[l
                    * z_dim1 + 1], ldz, (ftnlen)1, (ftnlen)1, (ftnlen)1);
        }

        d__[l] -= p;
        e[l] = g;
        goto L40;

/*        Eigenvalue found. */

L80:
        d__[l] = p;

        ++l;
        if (l <= lend) {
            goto L40;
        }
        goto L140;

    } else {

/*        QR Iteration */

/*        Look for small superdiagonal element. */

L90:
        if (l != lend) {
            lendp1 = lend + 1;
            i__1 = lendp1;
            for (m = l; m >= i__1; --m) {
/* Computing 2nd power */
                d__2 = (d__1 = e[m - 1], abs(d__1));
                tst = d__2 * d__2;
                if (tst <= eps2 * (d__1 = d__[m], abs(d__1)) * (d__2 = d__[m
                        - 1], abs(d__2)) + safmin) {
                    goto L110;
                }
/* L100: */
            }
        }

        m = lend;

L110:
        if (m > lend) {
            e[m - 1] = 0.;
        }
        p = d__[l];
        if (m == l) {
            goto L130;
        }

/*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2 */
/*        to compute its eigensystem. */

        if (m == l - 1) {
            if (icompz > 0) {
                dlaev2_(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2, &c__, &s)
                        ;
                work[m] = c__;
                work[*n - 1 + m] = s;
                dlasr_((char *)"R", (char *)"V", (char *)"F", n, &c__2, &work[m], &work[*n - 1 + m], &
                        z__[(l - 1) * z_dim1 + 1], ldz, (ftnlen)1, (ftnlen)1,
                        (ftnlen)1);
            } else {
                dlae2_(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2);
            }
            d__[l - 1] = rt1;
            d__[l] = rt2;
            e[l - 1] = 0.;
            l += -2;
            if (l >= lend) {
                goto L90;
            }
            goto L140;
        }

        if (jtot == nmaxit) {
            goto L140;
        }
        ++jtot;

/*        Form shift. */

        g = (d__[l - 1] - p) / (e[l - 1] * 2.);
        r__ = dlapy2_(&g, &c_b10);
        g = d__[m] - p + e[l - 1] / (g + d_sign(&r__, &g));

        s = 1.;
        c__ = 1.;
        p = 0.;

/*        Inner loop */

        lm1 = l - 1;
        i__1 = lm1;
        for (i__ = m; i__ <= i__1; ++i__) {
            f = s * e[i__];
            b = c__ * e[i__];
            dlartg_(&g, &f, &c__, &s, &r__);
            if (i__ != m) {
                e[i__ - 1] = r__;
            }
            g = d__[i__] - p;
            r__ = (d__[i__ + 1] - g) * s + c__ * 2. * b;
            p = s * r__;
            d__[i__] = g + p;
            g = c__ * r__ - b;

/*           If eigenvectors are desired, then save rotations. */

            if (icompz > 0) {
                work[i__] = c__;
                work[*n - 1 + i__] = s;
            }

/* L120: */
        }

/*        If eigenvectors are desired, then apply saved rotations. */

        if (icompz > 0) {
            mm = l - m + 1;
            dlasr_((char *)"R", (char *)"V", (char *)"F", n, &mm, &work[m], &work[*n - 1 + m], &z__[m
                    * z_dim1 + 1], ldz, (ftnlen)1, (ftnlen)1, (ftnlen)1);
        }

        d__[l] -= p;
        e[lm1] = g;
        goto L90;

/*        Eigenvalue found. */

L130:
        d__[l] = p;

        --l;
        if (l >= lend) {
            goto L90;
        }
        goto L140;

    }

/*     Undo scaling if necessary */

L140:
    if (iscale == 1) {
        i__1 = lendsv - lsv + 1;
        dlascl_((char *)"G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d__[lsv],
                n, info, (ftnlen)1);
        i__1 = lendsv - lsv;
        dlascl_((char *)"G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &e[lsv], n,
                info, (ftnlen)1);
    } else if (iscale == 2) {
        i__1 = lendsv - lsv + 1;
        dlascl_((char *)"G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d__[lsv],
                n, info, (ftnlen)1);
        i__1 = lendsv - lsv;
        dlascl_((char *)"G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &e[lsv], n,
                info, (ftnlen)1);
    }

/*     Check for no convergence to an eigenvalue after a total */
/*     of N*MAXIT iterations. */

    if (jtot < nmaxit) {
        goto L10;
    }
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
        if (e[i__] != 0.) {
            ++(*info);
        }
/* L150: */
    }
    goto L190;

/*     Order eigenvalues and eigenvectors. */

L160:
    if (icompz == 0) {

/*        Use Quick Sort */

        dlasrt_((char *)"I", n, &d__[1], info, (ftnlen)1);

    } else {

/*        Use Selection Sort to minimize swaps of eigenvectors */

        i__1 = *n;
        for (ii = 2; ii <= i__1; ++ii) {
            i__ = ii - 1;
            k = i__;
            p = d__[i__];
            i__2 = *n;
            for (j = ii; j <= i__2; ++j) {
                if (d__[j] < p) {
                    k = j;
                    p = d__[j];
                }
/* L170: */
            }
            if (k != i__) {
                d__[k] = d__[i__];
                d__[i__] = p;
                dswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[k * z_dim1 + 1],
                         &c__1);
            }
/* L180: */
        }
    }

L190:
    return 0;

/*     End of DSTEQR */

} /* dsteqr_ */

#ifdef __cplusplus
        }
#endif
