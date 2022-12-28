/* fortran/zhetd2.f -- translated by f2c (version 20200916).
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

static doublecomplex c_b2 = {0.,0.};
static integer c__1 = 1;

/* > \brief \b ZHETD2 reduces a Hermitian matrix to real symmetric tridiagonal form by an unitary similarity t
ransformation (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHETD2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetd2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetd2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetd2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHETD2( UPLO, N, A, LDA, D, E, TAU, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   D( * ), E( * ) */
/*       COMPLEX*16         A( LDA, * ), TAU( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHETD2 reduces a complex Hermitian matrix A to real symmetric */
/* > tridiagonal form T by a unitary similarity transformation: */
/* > Q**H * A * Q = T. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the upper or lower triangular part of the */
/* >          Hermitian matrix A is stored: */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading */
/* >          n-by-n upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading n-by-n lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* >          On exit, if UPLO = 'U', the diagonal and first superdiagonal */
/* >          of A are overwritten by the corresponding elements of the */
/* >          tridiagonal matrix T, and the elements above the first */
/* >          superdiagonal, with the array TAU, represent the unitary */
/* >          matrix Q as a product of elementary reflectors; if UPLO */
/* >          = 'L', the diagonal and first subdiagonal of A are over- */
/* >          written by the corresponding elements of the tridiagonal */
/* >          matrix T, and the elements below the first subdiagonal, with */
/* >          the array TAU, represent the unitary matrix Q as a product */
/* >          of elementary reflectors. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The diagonal elements of the tridiagonal matrix T: */
/* >          D(i) = A(i,i). */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >          The off-diagonal elements of the tridiagonal matrix T: */
/* >          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX*16 array, dimension (N-1) */
/* >          The scalar factors of the elementary reflectors (see Further */
/* >          Details). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complex16HEcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  If UPLO = 'U', the matrix Q is represented as a product of elementary */
/* >  reflectors */
/* > */
/* >     Q = H(n-1) . . . H(2) H(1). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a complex vector with */
/* >  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in */
/* >  A(1:i-1,i+1), and tau in TAU(i). */
/* > */
/* >  If UPLO = 'L', the matrix Q is represented as a product of elementary */
/* >  reflectors */
/* > */
/* >     Q = H(1) H(2) . . . H(n-1). */
/* > */
/* >  Each H(i) has the form */
/* > */
/* >     H(i) = I - tau * v * v**H */
/* > */
/* >  where tau is a complex scalar, and v is a complex vector with */
/* >  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i), */
/* >  and tau in TAU(i). */
/* > */
/* >  The contents of A on exit are illustrated by the following examples */
/* >  with n = 5: */
/* > */
/* >  if UPLO = 'U':                       if UPLO = 'L': */
/* > */
/* >    (  d   e   v2  v3  v4 )              (  d                  ) */
/* >    (      d   e   v3  v4 )              (  e   d              ) */
/* >    (          d   e   v4 )              (  v1  e   d          ) */
/* >    (              d   e  )              (  v1  v2  e   d      ) */
/* >    (                  d  )              (  v1  v2  v3  e   d  ) */
/* > */
/* >  where d and e denote diagonal and off-diagonal elements of T, and vi */
/* >  denotes an element of the vector defining H(i). */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zhetd2_(char *uplo, integer *n, doublecomplex *a,
        integer *lda, doublereal *d__, doublereal *e, doublecomplex *tau,
        integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Local variables */
    integer i__;
    doublecomplex taui;
    extern /* Subroutine */ int zher2_(char *, integer *, doublecomplex *,
            doublecomplex *, integer *, doublecomplex *, integer *,
            doublecomplex *, integer *, ftnlen);
    doublecomplex alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *,
            doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int zhemv_(char *, integer *, doublecomplex *,
            doublecomplex *, integer *, doublecomplex *, integer *,
            doublecomplex *, doublecomplex *, integer *, ftnlen);
    logical upper;
    extern /* Subroutine */ int zaxpy_(integer *, doublecomplex *,
            doublecomplex *, integer *, doublecomplex *, integer *), xerbla_(
            char *, integer *, ftnlen), zlarfg_(integer *, doublecomplex *,
            doublecomplex *, integer *, doublecomplex *);


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

/*     Test the input parameters */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tau;

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
        xerbla_((char *)"ZHETD2", &i__1, (ftnlen)6);
        return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
        return 0;
    }

    if (upper) {

/*        Reduce the upper triangle of A */

        i__1 = *n + *n * a_dim1;
        i__2 = *n + *n * a_dim1;
        d__1 = a[i__2].r;
        a[i__1].r = d__1, a[i__1].i = 0.;
        for (i__ = *n - 1; i__ >= 1; --i__) {

/*           Generate elementary reflector H(i) = I - tau * v * v**H */
/*           to annihilate A(1:i-1,i+1) */

            i__1 = i__ + (i__ + 1) * a_dim1;
            alpha.r = a[i__1].r, alpha.i = a[i__1].i;
            zlarfg_(&i__, &alpha, &a[(i__ + 1) * a_dim1 + 1], &c__1, &taui);
            e[i__] = alpha.r;

            if (taui.r != 0. || taui.i != 0.) {

/*              Apply H(i) from both sides to A(1:i,1:i) */

                i__1 = i__ + (i__ + 1) * a_dim1;
                a[i__1].r = 1., a[i__1].i = 0.;

/*              Compute  x := tau * A * v  storing x in TAU(1:i) */

                zhemv_(uplo, &i__, &taui, &a[a_offset], lda, &a[(i__ + 1) *
                        a_dim1 + 1], &c__1, &c_b2, &tau[1], &c__1, (ftnlen)1);

/*              Compute  w := x - 1/2 * tau * (x**H * v) * v */

                z__3.r = -.5, z__3.i = -0.;
                z__2.r = z__3.r * taui.r - z__3.i * taui.i, z__2.i = z__3.r *
                        taui.i + z__3.i * taui.r;
                zdotc_(&z__4, &i__, &tau[1], &c__1, &a[(i__ + 1) * a_dim1 + 1]
                        , &c__1);
                z__1.r = z__2.r * z__4.r - z__2.i * z__4.i, z__1.i = z__2.r *
                        z__4.i + z__2.i * z__4.r;
                alpha.r = z__1.r, alpha.i = z__1.i;
                zaxpy_(&i__, &alpha, &a[(i__ + 1) * a_dim1 + 1], &c__1, &tau[
                        1], &c__1);

/*              Apply the transformation as a rank-2 update: */
/*                 A := A - v * w**H - w * v**H */

                z__1.r = -1., z__1.i = -0.;
                zher2_(uplo, &i__, &z__1, &a[(i__ + 1) * a_dim1 + 1], &c__1, &
                        tau[1], &c__1, &a[a_offset], lda, (ftnlen)1);

            } else {
                i__1 = i__ + i__ * a_dim1;
                i__2 = i__ + i__ * a_dim1;
                d__1 = a[i__2].r;
                a[i__1].r = d__1, a[i__1].i = 0.;
            }
            i__1 = i__ + (i__ + 1) * a_dim1;
            i__2 = i__;
            a[i__1].r = e[i__2], a[i__1].i = 0.;
            i__1 = i__ + 1 + (i__ + 1) * a_dim1;
            d__[i__ + 1] = a[i__1].r;
            i__1 = i__;
            tau[i__1].r = taui.r, tau[i__1].i = taui.i;
/* L10: */
        }
        i__1 = a_dim1 + 1;
        d__[1] = a[i__1].r;
    } else {

/*        Reduce the lower triangle of A */

        i__1 = a_dim1 + 1;
        i__2 = a_dim1 + 1;
        d__1 = a[i__2].r;
        a[i__1].r = d__1, a[i__1].i = 0.;
        i__1 = *n - 1;
        for (i__ = 1; i__ <= i__1; ++i__) {

/*           Generate elementary reflector H(i) = I - tau * v * v**H */
/*           to annihilate A(i+2:n,i) */

            i__2 = i__ + 1 + i__ * a_dim1;
            alpha.r = a[i__2].r, alpha.i = a[i__2].i;
            i__2 = *n - i__;
/* Computing MIN */
            i__3 = i__ + 2;
            zlarfg_(&i__2, &alpha, &a[min(i__3,*n) + i__ * a_dim1], &c__1, &
                    taui);
            e[i__] = alpha.r;

            if (taui.r != 0. || taui.i != 0.) {

/*              Apply H(i) from both sides to A(i+1:n,i+1:n) */

                i__2 = i__ + 1 + i__ * a_dim1;
                a[i__2].r = 1., a[i__2].i = 0.;

/*              Compute  x := tau * A * v  storing y in TAU(i:n-1) */

                i__2 = *n - i__;
                zhemv_(uplo, &i__2, &taui, &a[i__ + 1 + (i__ + 1) * a_dim1],
                        lda, &a[i__ + 1 + i__ * a_dim1], &c__1, &c_b2, &tau[
                        i__], &c__1, (ftnlen)1);

/*              Compute  w := x - 1/2 * tau * (x**H * v) * v */

                z__3.r = -.5, z__3.i = -0.;
                z__2.r = z__3.r * taui.r - z__3.i * taui.i, z__2.i = z__3.r *
                        taui.i + z__3.i * taui.r;
                i__2 = *n - i__;
                zdotc_(&z__4, &i__2, &tau[i__], &c__1, &a[i__ + 1 + i__ *
                        a_dim1], &c__1);
                z__1.r = z__2.r * z__4.r - z__2.i * z__4.i, z__1.i = z__2.r *
                        z__4.i + z__2.i * z__4.r;
                alpha.r = z__1.r, alpha.i = z__1.i;
                i__2 = *n - i__;
                zaxpy_(&i__2, &alpha, &a[i__ + 1 + i__ * a_dim1], &c__1, &tau[
                        i__], &c__1);

/*              Apply the transformation as a rank-2 update: */
/*                 A := A - v * w**H - w * v**H */

                i__2 = *n - i__;
                z__1.r = -1., z__1.i = -0.;
                zher2_(uplo, &i__2, &z__1, &a[i__ + 1 + i__ * a_dim1], &c__1,
                        &tau[i__], &c__1, &a[i__ + 1 + (i__ + 1) * a_dim1],
                        lda, (ftnlen)1);

            } else {
                i__2 = i__ + 1 + (i__ + 1) * a_dim1;
                i__3 = i__ + 1 + (i__ + 1) * a_dim1;
                d__1 = a[i__3].r;
                a[i__2].r = d__1, a[i__2].i = 0.;
            }
            i__2 = i__ + 1 + i__ * a_dim1;
            i__3 = i__;
            a[i__2].r = e[i__3], a[i__2].i = 0.;
            i__2 = i__ + i__ * a_dim1;
            d__[i__] = a[i__2].r;
            i__2 = i__;
            tau[i__2].r = taui.r, tau[i__2].i = taui.i;
/* L20: */
        }
        i__1 = *n + *n * a_dim1;
        d__[*n] = a[i__1].r;
    }

    return 0;

/*     End of ZHETD2 */

} /* zhetd2_ */

#ifdef __cplusplus
        }
#endif
