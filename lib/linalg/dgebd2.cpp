/* fortran/dgebd2.f -- translated by f2c (version 20200916).
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

static integer c__1 = 1;

/* > \brief \b DGEBD2 reduces a general matrix to bidiagonal form using an unblocked algorithm. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEBD2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgebd2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgebd2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgebd2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEBD2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAUP( * ), */
/*      $                   TAUQ( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGEBD2 reduces a real general m by n matrix A to upper or lower */
/* > bidiagonal form B by an orthogonal transformation: Q**T * A * P = B. */
/* > */
/* > If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows in the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns in the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the m by n general matrix to be reduced. */
/* >          On exit, */
/* >          if m >= n, the diagonal and the first superdiagonal are */
/* >            overwritten with the upper bidiagonal matrix B; the */
/* >            elements below the diagonal, with the array TAUQ, represent */
/* >            the orthogonal matrix Q as a product of elementary */
/* >            reflectors, and the elements above the first superdiagonal, */
/* >            with the array TAUP, represent the orthogonal matrix P as */
/* >            a product of elementary reflectors; */
/* >          if m < n, the diagonal and the first subdiagonal are */
/* >            overwritten with the lower bidiagonal matrix B; the */
/* >            elements below the first subdiagonal, with the array TAUQ, */
/* >            represent the orthogonal matrix Q as a product of */
/* >            elementary reflectors, and the elements above the diagonal, */
/* >            with the array TAUP, represent the orthogonal matrix P as */
/* >            a product of elementary reflectors. */
/* >          See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (min(M,N)) */
/* >          The diagonal elements of the bidiagonal matrix B: */
/* >          D(i) = A(i,i). */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (min(M,N)-1) */
/* >          The off-diagonal elements of the bidiagonal matrix B: */
/* >          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1; */
/* >          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUQ */
/* > \verbatim */
/* >          TAUQ is DOUBLE PRECISION array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors which */
/* >          represent the orthogonal matrix Q. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP */
/* > \verbatim */
/* >          TAUP is DOUBLE PRECISION array, dimension (min(M,N)) */
/* >          The scalar factors of the elementary reflectors which */
/* >          represent the orthogonal matrix P. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (max(M,N)) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit. */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doubleGEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrices Q and P are represented as products of elementary */
/* >  reflectors: */
/* > */
/* >  If m >= n, */
/* > */
/* >     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1) */
/* > */
/* >  Each H(i) and G(i) has the form: */
/* > */
/* >     H(i) = I - tauq * v * v**T  and G(i) = I - taup * u * u**T */
/* > */
/* >  where tauq and taup are real scalars, and v and u are real vectors; */
/* >  v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i); */
/* >  u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n); */
/* >  tauq is stored in TAUQ(i) and taup in TAUP(i). */
/* > */
/* >  If m < n, */
/* > */
/* >     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m) */
/* > */
/* >  Each H(i) and G(i) has the form: */
/* > */
/* >     H(i) = I - tauq * v * v**T  and G(i) = I - taup * u * u**T */
/* > */
/* >  where tauq and taup are real scalars, and v and u are real vectors; */
/* >  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i); */
/* >  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n); */
/* >  tauq is stored in TAUQ(i) and taup in TAUP(i). */
/* > */
/* >  The contents of A on exit are illustrated by the following examples: */
/* > */
/* >  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n): */
/* > */
/* >    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 ) */
/* >    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 ) */
/* >    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 ) */
/* >    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 ) */
/* >    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 ) */
/* >    (  v1  v2  v3  v4  v5 ) */
/* > */
/* >  where d and e denote diagonal and off-diagonal elements of B, vi */
/* >  denotes an element of the vector defining H(i), and ui an element of */
/* >  the vector defining G(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dgebd2_(integer *m, integer *n, doublereal *a, integer *
        lda, doublereal *d__, doublereal *e, doublereal *tauq, doublereal *
        taup, doublereal *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    integer i__;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *,
            doublereal *, integer *, doublereal *, doublereal *, integer *,
            doublereal *, ftnlen), dlarfg_(integer *, doublereal *,
            doublereal *, integer *, doublereal *), xerbla_(char *, integer *,
             ftnlen);


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tauq;
    --taup;
    --work;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*lda < max(1,*m)) {
        *info = -4;
    }
    if (*info < 0) {
        i__1 = -(*info);
        xerbla_((char *)"DGEBD2", &i__1, (ftnlen)6);
        return 0;
    }

    if (*m >= *n) {

/*        Reduce to upper bidiagonal form */

        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {

/*           Generate elementary reflector H(i) to annihilate A(i+1:m,i) */

            i__2 = *m - i__ + 1;
/* Computing MIN */
            i__3 = i__ + 1;
            dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[min(i__3,*m) + i__ *
                    a_dim1], &c__1, &tauq[i__]);
            d__[i__] = a[i__ + i__ * a_dim1];
            a[i__ + i__ * a_dim1] = 1.;

/*           Apply H(i) to A(i:m,i+1:n) from the left */

            if (i__ < *n) {
                i__2 = *m - i__ + 1;
                i__3 = *n - i__;
                dlarf_((char *)"Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &
                        tauq[i__], &a[i__ + (i__ + 1) * a_dim1], lda, &work[1]
                        , (ftnlen)4);
            }
            a[i__ + i__ * a_dim1] = d__[i__];

            if (i__ < *n) {

/*              Generate elementary reflector G(i) to annihilate */
/*              A(i,i+2:n) */

                i__2 = *n - i__;
/* Computing MIN */
                i__3 = i__ + 2;
                dlarfg_(&i__2, &a[i__ + (i__ + 1) * a_dim1], &a[i__ + min(
                        i__3,*n) * a_dim1], lda, &taup[i__]);
                e[i__] = a[i__ + (i__ + 1) * a_dim1];
                a[i__ + (i__ + 1) * a_dim1] = 1.;

/*              Apply G(i) to A(i+1:m,i+1:n) from the right */

                i__2 = *m - i__;
                i__3 = *n - i__;
                dlarf_((char *)"Right", &i__2, &i__3, &a[i__ + (i__ + 1) * a_dim1],
                        lda, &taup[i__], &a[i__ + 1 + (i__ + 1) * a_dim1],
                        lda, &work[1], (ftnlen)5);
                a[i__ + (i__ + 1) * a_dim1] = e[i__];
            } else {
                taup[i__] = 0.;
            }
/* L10: */
        }
    } else {

/*        Reduce to lower bidiagonal form */

        i__1 = *m;
        for (i__ = 1; i__ <= i__1; ++i__) {

/*           Generate elementary reflector G(i) to annihilate A(i,i+1:n) */

            i__2 = *n - i__ + 1;
/* Computing MIN */
            i__3 = i__ + 1;
            dlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[i__ + min(i__3,*n) *
                    a_dim1], lda, &taup[i__]);
            d__[i__] = a[i__ + i__ * a_dim1];
            a[i__ + i__ * a_dim1] = 1.;

/*           Apply G(i) to A(i+1:m,i:n) from the right */

            if (i__ < *m) {
                i__2 = *m - i__;
                i__3 = *n - i__ + 1;
                dlarf_((char *)"Right", &i__2, &i__3, &a[i__ + i__ * a_dim1], lda, &
                        taup[i__], &a[i__ + 1 + i__ * a_dim1], lda, &work[1],
                        (ftnlen)5);
            }
            a[i__ + i__ * a_dim1] = d__[i__];

            if (i__ < *m) {

/*              Generate elementary reflector H(i) to annihilate */
/*              A(i+2:m,i) */

                i__2 = *m - i__;
/* Computing MIN */
                i__3 = i__ + 2;
                dlarfg_(&i__2, &a[i__ + 1 + i__ * a_dim1], &a[min(i__3,*m) +
                        i__ * a_dim1], &c__1, &tauq[i__]);
                e[i__] = a[i__ + 1 + i__ * a_dim1];
                a[i__ + 1 + i__ * a_dim1] = 1.;

/*              Apply H(i) to A(i+1:m,i+1:n) from the left */

                i__2 = *m - i__;
                i__3 = *n - i__;
                dlarf_((char *)"Left", &i__2, &i__3, &a[i__ + 1 + i__ * a_dim1], &
                        c__1, &tauq[i__], &a[i__ + 1 + (i__ + 1) * a_dim1],
                        lda, &work[1], (ftnlen)4);
                a[i__ + 1 + i__ * a_dim1] = e[i__];
            } else {
                tauq[i__] = 0.;
            }
/* L20: */
        }
    }
    return 0;

/*     End of DGEBD2 */

} /* dgebd2_ */

#ifdef __cplusplus
        }
#endif
