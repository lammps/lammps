/* fortran/dgetrf.f -- translated by f2c (version 20200916).
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
static integer c_n1 = -1;
static doublereal c_b16 = 1.;
static doublereal c_b19 = -1.;

/* > \brief \b DGETRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGETRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGETRF computes an LU factorization of a general M-by-N matrix A */
/* > using partial pivoting with row interchanges. */
/* > */
/* > The factorization has the form */
/* >    A = P * L * U */
/* > where P is a permutation matrix, L is lower triangular with unit */
/* > diagonal elements (lower trapezoidal if m > n), and U is upper */
/* > triangular (upper trapezoidal if m < n). */
/* > */
/* > This is the right-looking Level 3 BLAS version of the algorithm. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix to be factored. */
/* >          On exit, the factors L and U from the factorization */
/* >          A = P*L*U; the unit diagonal elements of L are not stored. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (min(M,N)) */
/* >          The pivot indices; for 1 <= i <= min(M,N), row i of the */
/* >          matrix was interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization */
/* >                has been completed, but the factor U is exactly */
/* >                singular, and division by zero will occur if it is used */
/* >                to solve a system of equations. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doubleGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int dgetrf_(integer *m, integer *n, doublereal *a, integer *
        lda, integer *ipiv, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    integer i__, j, jb, nb;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *,
            integer *, doublereal *, doublereal *, integer *, doublereal *,
            integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    integer iinfo;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *,
            integer *, integer *, doublereal *, doublereal *, integer *,
            doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), xerbla_(
            char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *,
            integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlaswp_(integer *, doublereal *, integer *,
            integer *, integer *, integer *, integer *), dgetrf2_(integer *,
            integer *, doublereal *, integer *, integer *, integer *);


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
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*lda < max(1,*m)) {
        *info = -4;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DGETRF", &i__1, (ftnlen)6);
        return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
        return 0;
    }

/*     Determine the block size for this environment. */

    nb = ilaenv_(&c__1, (char *)"DGETRF", (char *)" ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)
            1);
    if (nb <= 1 || nb >= min(*m,*n)) {

/*        Use unblocked code. */

        dgetrf2_(m, n, &a[a_offset], lda, &ipiv[1], info);
    } else {

/*        Use blocked code. */

        i__1 = min(*m,*n);
        i__2 = nb;
        for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
            i__3 = min(*m,*n) - j + 1;
            jb = min(i__3,nb);

/*           Factor diagonal and subdiagonal blocks and test for exact */
/*           singularity. */

            i__3 = *m - j + 1;
            dgetrf2_(&i__3, &jb, &a[j + j * a_dim1], lda, &ipiv[j], &iinfo);

/*           Adjust INFO and the pivot indices. */

            if (*info == 0 && iinfo > 0) {
                *info = iinfo + j - 1;
            }
/* Computing MIN */
            i__4 = *m, i__5 = j + jb - 1;
            i__3 = min(i__4,i__5);
            for (i__ = j; i__ <= i__3; ++i__) {
                ipiv[i__] = j - 1 + ipiv[i__];
/* L10: */
            }

/*           Apply interchanges to columns 1:J-1. */

            i__3 = j - 1;
            i__4 = j + jb - 1;
            dlaswp_(&i__3, &a[a_offset], lda, &j, &i__4, &ipiv[1], &c__1);

            if (j + jb <= *n) {

/*              Apply interchanges to columns J+JB:N. */

                i__3 = *n - j - jb + 1;
                i__4 = j + jb - 1;
                dlaswp_(&i__3, &a[(j + jb) * a_dim1 + 1], lda, &j, &i__4, &
                        ipiv[1], &c__1);

/*              Compute block row of U. */

                i__3 = *n - j - jb + 1;
                dtrsm_((char *)"Left", (char *)"Lower", (char *)"No transpose", (char *)"Unit", &jb, &i__3, &
                        c_b16, &a[j + j * a_dim1], lda, &a[j + (j + jb) *
                        a_dim1], lda, (ftnlen)4, (ftnlen)5, (ftnlen)12, (
                        ftnlen)4);
                if (j + jb <= *m) {

/*                 Update trailing submatrix. */

                    i__3 = *m - j - jb + 1;
                    i__4 = *n - j - jb + 1;
                    dgemm_((char *)"No transpose", (char *)"No transpose", &i__3, &i__4, &jb,
                            &c_b19, &a[j + jb + j * a_dim1], lda, &a[j + (j +
                            jb) * a_dim1], lda, &c_b16, &a[j + jb + (j + jb) *
                             a_dim1], lda, (ftnlen)12, (ftnlen)12);
                }
            }
/* L20: */
        }
    }
    return 0;

/*     End of DGETRF */

} /* dgetrf_ */

#ifdef __cplusplus
        }
#endif
