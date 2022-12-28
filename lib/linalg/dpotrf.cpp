/* fortran/dpotrf.f -- translated by f2c (version 20200916).
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
static doublereal c_b13 = -1.;
static doublereal c_b14 = 1.;

/* > \brief \b DPOTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPOTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpotrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpotrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpotrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPOTRF computes the Cholesky factorization of a real symmetric */
/* > positive definite matrix A. */
/* > */
/* > The factorization has the form */
/* >    A = U**T * U,  if UPLO = 'U', or */
/* >    A = L  * L**T,  if UPLO = 'L', */
/* > where U is an upper triangular matrix and L is lower triangular. */
/* > */
/* > This is the block version of the algorithm, calling Level 3 BLAS. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored; */
/* >          = 'L':  Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
/* >          N-by-N upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading N-by-N lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* > */
/* >          On exit, if INFO = 0, the factor U or L from the Cholesky */
/* >          factorization A = U**T*U or A = L*L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, the leading minor of order i is not */
/* >                positive definite, and the factorization could not be */
/* >                completed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doublePOcomputational */

/*  ===================================================================== */
/* Subroutine */ int dpotrf_(char *uplo, integer *n, doublereal *a, integer *
        lda, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    integer j, jb, nb;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *,
            integer *, doublereal *, doublereal *, integer *, doublereal *,
            integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *,
            integer *, integer *, doublereal *, doublereal *, integer *,
            doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    logical upper;
    extern /* Subroutine */ int dsyrk_(char *, char *, integer *, integer *,
            doublereal *, doublereal *, integer *, doublereal *, doublereal *,
             integer *, ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *,
            integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dpotrf2_(char *, integer *, doublereal *,
            integer *, integer *, ftnlen);


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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);
    if (! upper && ! lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*lda < max(1,*n)) {
        *info = -4;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DPOTRF", &i__1, (ftnlen)6);
        return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
        return 0;
    }

/*     Determine the block size for this environment. */

    nb = ilaenv_(&c__1, (char *)"DPOTRF", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
            ftnlen)1);
    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code. */

        dpotrf2_(uplo, n, &a[a_offset], lda, info, (ftnlen)1);
    } else {

/*        Use blocked code. */

        if (upper) {

/*           Compute the Cholesky factorization A = U**T*U. */

            i__1 = *n;
            i__2 = nb;
            for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {

/*              Update and factorize the current diagonal block and test */
/*              for non-positive-definiteness. */

/* Computing MIN */
                i__3 = nb, i__4 = *n - j + 1;
                jb = min(i__3,i__4);
                i__3 = j - 1;
                dsyrk_((char *)"Upper", (char *)"Transpose", &jb, &i__3, &c_b13, &a[j *
                        a_dim1 + 1], lda, &c_b14, &a[j + j * a_dim1], lda, (
                        ftnlen)5, (ftnlen)9);
                dpotrf2_((char *)"Upper", &jb, &a[j + j * a_dim1], lda, info, (ftnlen)
                        5);
                if (*info != 0) {
                    goto L30;
                }
                if (j + jb <= *n) {

/*                 Compute the current block row. */

                    i__3 = *n - j - jb + 1;
                    i__4 = j - 1;
                    dgemm_((char *)"Transpose", (char *)"No transpose", &jb, &i__3, &i__4, &
                            c_b13, &a[j * a_dim1 + 1], lda, &a[(j + jb) *
                            a_dim1 + 1], lda, &c_b14, &a[j + (j + jb) *
                            a_dim1], lda, (ftnlen)9, (ftnlen)12);
                    i__3 = *n - j - jb + 1;
                    dtrsm_((char *)"Left", (char *)"Upper", (char *)"Transpose", (char *)"Non-unit", &jb, &
                            i__3, &c_b14, &a[j + j * a_dim1], lda, &a[j + (j
                            + jb) * a_dim1], lda, (ftnlen)4, (ftnlen)5, (
                            ftnlen)9, (ftnlen)8);
                }
/* L10: */
            }

        } else {

/*           Compute the Cholesky factorization A = L*L**T. */

            i__2 = *n;
            i__1 = nb;
            for (j = 1; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {

/*              Update and factorize the current diagonal block and test */
/*              for non-positive-definiteness. */

/* Computing MIN */
                i__3 = nb, i__4 = *n - j + 1;
                jb = min(i__3,i__4);
                i__3 = j - 1;
                dsyrk_((char *)"Lower", (char *)"No transpose", &jb, &i__3, &c_b13, &a[j +
                        a_dim1], lda, &c_b14, &a[j + j * a_dim1], lda, (
                        ftnlen)5, (ftnlen)12);
                dpotrf2_((char *)"Lower", &jb, &a[j + j * a_dim1], lda, info, (ftnlen)
                        5);
                if (*info != 0) {
                    goto L30;
                }
                if (j + jb <= *n) {

/*                 Compute the current block column. */

                    i__3 = *n - j - jb + 1;
                    i__4 = j - 1;
                    dgemm_((char *)"No transpose", (char *)"Transpose", &i__3, &jb, &i__4, &
                            c_b13, &a[j + jb + a_dim1], lda, &a[j + a_dim1],
                            lda, &c_b14, &a[j + jb + j * a_dim1], lda, (
                            ftnlen)12, (ftnlen)9);
                    i__3 = *n - j - jb + 1;
                    dtrsm_((char *)"Right", (char *)"Lower", (char *)"Transpose", (char *)"Non-unit", &i__3, &
                            jb, &c_b14, &a[j + j * a_dim1], lda, &a[j + jb +
                            j * a_dim1], lda, (ftnlen)5, (ftnlen)5, (ftnlen)9,
                             (ftnlen)8);
                }
/* L20: */
            }
        }
    }
    goto L40;

L30:
    *info = *info + j - 1;

L40:
    return 0;

/*     End of DPOTRF */

} /* dpotrf_ */

#ifdef __cplusplus
        }
#endif
