/* fortran/zungl2.f -- translated by f2c (version 20200916).
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

/* > \brief \b ZUNGL2 generates all or part of the unitary matrix Q from an LQ factorization determined by cge
lqf (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZUNGL2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zungl2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zungl2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zungl2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZUNGL2( M, N, K, A, LDA, TAU, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, K, LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZUNGL2 generates an m-by-n complex matrix Q with orthonormal rows, */
/* > which is defined as the first m rows of a product of k elementary */
/* > reflectors of order n */
/* > */
/* >       Q  =  H(k)**H . . . H(2)**H H(1)**H */
/* > */
/* > as returned by ZGELQF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix Q. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix Q. N >= M. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of elementary reflectors whose product defines the */
/* >          matrix Q. M >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the i-th row must contain the vector which defines */
/* >          the elementary reflector H(i), for i = 1,2,...,k, as returned */
/* >          by ZGELQF in the first k rows of its array argument A. */
/* >          On exit, the m by n matrix Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The first dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX*16 array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by ZGELQF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (M) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument has an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zungl2_(integer *m, integer *n, integer *k,
        doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
        work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    integer i__, j, l;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *,
            doublecomplex *, integer *), zlarf_(char *, integer *, integer *,
            doublecomplex *, integer *, doublecomplex *, doublecomplex *,
            integer *, doublecomplex *, ftnlen), xerbla_(char *, integer *,
            ftnlen), zlacgv_(integer *, doublecomplex *, integer *);


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

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
        *info = -1;
    } else if (*n < *m) {
        *info = -2;
    } else if (*k < 0 || *k > *m) {
        *info = -3;
    } else if (*lda < max(1,*m)) {
        *info = -5;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"ZUNGL2", &i__1, (ftnlen)6);
        return 0;
    }

/*     Quick return if possible */

    if (*m <= 0) {
        return 0;
    }

    if (*k < *m) {

/*        Initialise rows k+1:m to rows of the unit matrix */

        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__2 = *m;
            for (l = *k + 1; l <= i__2; ++l) {
                i__3 = l + j * a_dim1;
                a[i__3].r = 0., a[i__3].i = 0.;
/* L10: */
            }
            if (j > *k && j <= *m) {
                i__2 = j + j * a_dim1;
                a[i__2].r = 1., a[i__2].i = 0.;
            }
/* L20: */
        }
    }

    for (i__ = *k; i__ >= 1; --i__) {

/*        Apply H(i)**H to A(i:m,i:n) from the right */

        if (i__ < *n) {
            i__1 = *n - i__;
            zlacgv_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda);
            if (i__ < *m) {
                i__1 = i__ + i__ * a_dim1;
                a[i__1].r = 1., a[i__1].i = 0.;
                i__1 = *m - i__;
                i__2 = *n - i__ + 1;
                d_cnjg(&z__1, &tau[i__]);
                zlarf_((char *)"Right", &i__1, &i__2, &a[i__ + i__ * a_dim1], lda, &
                        z__1, &a[i__ + 1 + i__ * a_dim1], lda, &work[1], (
                        ftnlen)5);
            }
            i__1 = *n - i__;
            i__2 = i__;
            z__1.r = -tau[i__2].r, z__1.i = -tau[i__2].i;
            zscal_(&i__1, &z__1, &a[i__ + (i__ + 1) * a_dim1], lda);
            i__1 = *n - i__;
            zlacgv_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda);
        }
        i__1 = i__ + i__ * a_dim1;
        d_cnjg(&z__2, &tau[i__]);
        z__1.r = 1. - z__2.r, z__1.i = 0. - z__2.i;
        a[i__1].r = z__1.r, a[i__1].i = z__1.i;

/*        Set A(i,1:i-1) to zero */

        i__1 = i__ - 1;
        for (l = 1; l <= i__1; ++l) {
            i__2 = i__ + l * a_dim1;
            a[i__2].r = 0., a[i__2].i = 0.;
/* L30: */
        }
/* L40: */
    }
    return 0;

/*     End of ZUNGL2 */

} /* zungl2_ */

#ifdef __cplusplus
        }
#endif
