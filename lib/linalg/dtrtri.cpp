/* fortran/dtrtri.f -- translated by f2c (version 20200916).
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
static integer c__2 = 2;
static doublereal c_b18 = 1.;
static doublereal c_b22 = -1.;

/* > \brief \b DTRTRI */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTRTRI + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrtri.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrtri.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrtri.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, UPLO */
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
/* > DTRTRI computes the inverse of a real upper or lower triangular */
/* > matrix A. */
/* > */
/* > This is the Level 3 BLAS version of the algorithm. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  A is upper triangular; */
/* >          = 'L':  A is lower triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >          = 'N':  A is non-unit triangular; */
/* >          = 'U':  A is unit triangular. */
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
/* >          On entry, the triangular matrix A.  If UPLO = 'U', the */
/* >          leading N-by-N upper triangular part of the array A contains */
/* >          the upper triangular matrix, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading N-by-N lower triangular part of the array A contains */
/* >          the lower triangular matrix, and the strictly upper */
/* >          triangular part of A is not referenced.  If DIAG = 'U', the */
/* >          diagonal elements of A are also not referenced and are */
/* >          assumed to be 1. */
/* >          On exit, the (triangular) inverse of the original matrix, in */
/* >          the same storage format. */
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
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* >          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular */
/* >               matrix is singular and its inverse can not be computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dtrtri_(char *uplo, char *diag, integer *n, doublereal *
        a, integer *lda, integer *info, ftnlen uplo_len, ftnlen diag_len)
{
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, i__1, i__2[2], i__3, i__4, i__5;
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_lmp_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    integer j, jb, nb, nn;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dtrmm_(char *, char *, char *, char *,
            integer *, integer *, doublereal *, doublereal *, integer *,
            doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dtrsm_(
            char *, char *, char *, char *, integer *, integer *, doublereal *
            , doublereal *, integer *, doublereal *, integer *, ftnlen,
            ftnlen, ftnlen, ftnlen);
    logical upper;
    extern /* Subroutine */ int dtrti2_(char *, char *, integer *, doublereal
            *, integer *, integer *, ftnlen, ftnlen), xerbla_(char *, integer
            *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *,
            integer *, integer *, ftnlen, ftnlen);
    logical nounit;


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
    nounit = lsame_(diag, (char *)"N", (ftnlen)1, (ftnlen)1);
    if (! upper && ! lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (! nounit && ! lsame_(diag, (char *)"U", (ftnlen)1, (ftnlen)1)) {
        *info = -2;
    } else if (*n < 0) {
        *info = -3;
    } else if (*lda < max(1,*n)) {
        *info = -5;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DTRTRI", &i__1, (ftnlen)6);
        return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
        return 0;
    }

/*     Check for singularity if non-unit. */

    if (nounit) {
        i__1 = *n;
        for (*info = 1; *info <= i__1; ++(*info)) {
            if (a[*info + *info * a_dim1] == 0.) {
                return 0;
            }
/* L10: */
        }
        *info = 0;
    }

/*     Determine the block size for this environment. */

/* Writing concatenation */
    i__2[0] = 1, a__1[0] = uplo;
    i__2[1] = 1, a__1[1] = diag;
    s_lmp_cat(ch__1, a__1, i__2, &c__2, (ftnlen)2);
    nb = ilaenv_(&c__1, (char *)"DTRTRI", ch__1, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
            ftnlen)2);
    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code */

        dtrti2_(uplo, diag, n, &a[a_offset], lda, info, (ftnlen)1, (ftnlen)1);
    } else {

/*        Use blocked code */

        if (upper) {

/*           Compute inverse of upper triangular matrix */

            i__1 = *n;
            i__3 = nb;
            for (j = 1; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {
/* Computing MIN */
                i__4 = nb, i__5 = *n - j + 1;
                jb = min(i__4,i__5);

/*              Compute rows 1:j-1 of current block column */

                i__4 = j - 1;
                dtrmm_((char *)"Left", (char *)"Upper", (char *)"No transpose", diag, &i__4, &jb, &
                        c_b18, &a[a_offset], lda, &a[j * a_dim1 + 1], lda, (
                        ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)1);
                i__4 = j - 1;
                dtrsm_((char *)"Right", (char *)"Upper", (char *)"No transpose", diag, &i__4, &jb, &
                        c_b22, &a[j + j * a_dim1], lda, &a[j * a_dim1 + 1],
                        lda, (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)1);

/*              Compute inverse of current diagonal block */

                dtrti2_((char *)"Upper", diag, &jb, &a[j + j * a_dim1], lda, info, (
                        ftnlen)5, (ftnlen)1);
/* L20: */
            }
        } else {

/*           Compute inverse of lower triangular matrix */

            nn = (*n - 1) / nb * nb + 1;
            i__3 = -nb;
            for (j = nn; i__3 < 0 ? j >= 1 : j <= 1; j += i__3) {
/* Computing MIN */
                i__1 = nb, i__4 = *n - j + 1;
                jb = min(i__1,i__4);
                if (j + jb <= *n) {

/*                 Compute rows j+jb:n of current block column */

                    i__1 = *n - j - jb + 1;
                    dtrmm_((char *)"Left", (char *)"Lower", (char *)"No transpose", diag, &i__1, &jb,
                            &c_b18, &a[j + jb + (j + jb) * a_dim1], lda, &a[j
                            + jb + j * a_dim1], lda, (ftnlen)4, (ftnlen)5, (
                            ftnlen)12, (ftnlen)1);
                    i__1 = *n - j - jb + 1;
                    dtrsm_((char *)"Right", (char *)"Lower", (char *)"No transpose", diag, &i__1, &jb,
                             &c_b22, &a[j + j * a_dim1], lda, &a[j + jb + j *
                            a_dim1], lda, (ftnlen)5, (ftnlen)5, (ftnlen)12, (
                            ftnlen)1);
                }

/*              Compute inverse of current diagonal block */

                dtrti2_((char *)"Lower", diag, &jb, &a[j + j * a_dim1], lda, info, (
                        ftnlen)5, (ftnlen)1);
/* L30: */
            }
        }
    }

    return 0;

/*     End of DTRTRI */

} /* dtrtri_ */

#ifdef __cplusplus
        }
#endif
