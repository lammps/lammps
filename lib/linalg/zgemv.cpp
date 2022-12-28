/* fortran/zgemv.f -- translated by f2c (version 20200916).
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

/* > \brief \b ZGEMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       COMPLEX*16 ALPHA,BETA */
/*       INTEGER INCX,INCY,LDA,M,N */
/*       CHARACTER TRANS */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 A(LDA,*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEMV  performs one of the matrix-vector operations */
/* > */
/* >    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or */
/* > */
/* >    y := alpha*A**H*x + beta*y, */
/* > */
/* > where alpha and beta are scalars, x and y are vectors and A is an */
/* > m by n matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >           On entry, TRANS specifies the operation to be performed as */
/* >           follows: */
/* > */
/* >              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y. */
/* > */
/* >              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y. */
/* > */
/* >              TRANS = 'C' or 'c'   y := alpha*A**H*x + beta*y. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >           On entry, M specifies the number of rows of the matrix A. */
/* >           M must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the number of columns of the matrix A. */
/* >           N must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is COMPLEX*16 */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension ( LDA, N ) */
/* >           Before entry, the leading m by n part of the array A must */
/* >           contain the matrix of coefficients. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           max( 1, m ). */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array, dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n' */
/* >           and at least */
/* >           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise. */
/* >           Before entry, the incremented array X must contain the */
/* >           vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >           On entry, INCX specifies the increment for the elements of */
/* >           X. INCX must not be zero. */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is COMPLEX*16 */
/* >           On entry, BETA specifies the scalar beta. When BETA is */
/* >           supplied as zero then Y need not be set on input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is COMPLEX*16 array, dimension at least */
/* >           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n' */
/* >           and at least */
/* >           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise. */
/* >           Before entry with BETA non-zero, the incremented array Y */
/* >           must contain the vector y. On exit, Y is overwritten by the */
/* >           updated vector y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >           On entry, INCY specifies the increment for the elements of */
/* >           Y. INCY must not be zero. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complex16_blas_level2 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Level 2 Blas routine. */
/* >  The vector and matrix arguments are not referenced when N = 0, or M = 0 */
/* > */
/* >  -- Written on 22-October-1986. */
/* >     Jack Dongarra, Argonne National Lab. */
/* >     Jeremy Du Croz, Nag Central Office. */
/* >     Sven Hammarling, Nag Central Office. */
/* >     Richard Hanson, Sandia National Labs. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zgemv_(char *trans, integer *m, integer *n,
        doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *
        x, integer *incx, doublecomplex *beta, doublecomplex *y, integer *
        incy, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void d_lmp_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    integer i__, j, ix, iy, jx, jy, kx, ky, info;
    doublecomplex temp;
    integer lenx, leny;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    logical noconj;


/*  -- Reference BLAS level2 routine -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
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

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
    info = 0;
    if (! lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, (char *)"T", (
            ftnlen)1, (ftnlen)1) && ! lsame_(trans, (char *)"C", (ftnlen)1, (ftnlen)1)
            ) {
        info = 1;
    } else if (*m < 0) {
        info = 2;
    } else if (*n < 0) {
        info = 3;
    } else if (*lda < max(1,*m)) {
        info = 6;
    } else if (*incx == 0) {
        info = 8;
    } else if (*incy == 0) {
        info = 11;
    }
    if (info != 0) {
        xerbla_((char *)"ZGEMV ", &info, (ftnlen)6);
        return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || alpha->r == 0. && alpha->i == 0. && (beta->r ==
            1. && beta->i == 0.)) {
        return 0;
    }

    noconj = lsame_(trans, (char *)"T", (ftnlen)1, (ftnlen)1);

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
/*     up the start points in  X  and  Y. */

    if (lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1)) {
        lenx = *n;
        leny = *m;
    } else {
        lenx = *m;
        leny = *n;
    }
    if (*incx > 0) {
        kx = 1;
    } else {
        kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
        ky = 1;
    } else {
        ky = 1 - (leny - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

/*     First form  y := beta*y. */

    if (beta->r != 1. || beta->i != 0.) {
        if (*incy == 1) {
            if (beta->r == 0. && beta->i == 0.) {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    i__2 = i__;
                    y[i__2].r = 0., y[i__2].i = 0.;
/* L10: */
                }
            } else {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    i__2 = i__;
                    i__3 = i__;
                    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i,
                            z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
                            .r;
                    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
/* L20: */
                }
            }
        } else {
            iy = ky;
            if (beta->r == 0. && beta->i == 0.) {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    i__2 = iy;
                    y[i__2].r = 0., y[i__2].i = 0.;
                    iy += *incy;
/* L30: */
                }
            } else {
                i__1 = leny;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    i__2 = iy;
                    i__3 = iy;
                    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i,
                            z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
                            .r;
                    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
                    iy += *incy;
/* L40: */
                }
            }
        }
    }
    if (alpha->r == 0. && alpha->i == 0.) {
        return 0;
    }
    if (lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1)) {

/*        Form  y := alpha*A*x + y. */

        jx = kx;
        if (*incy == 1) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = jx;
                z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
                         alpha->r * x[i__2].i + alpha->i * x[i__2].r;
                temp.r = z__1.r, temp.i = z__1.i;
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    i__3 = i__;
                    i__4 = i__;
                    i__5 = i__ + j * a_dim1;
                    z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, z__2.i =
                             temp.r * a[i__5].i + temp.i * a[i__5].r;
                    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
                    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
/* L50: */
                }
                jx += *incx;
/* L60: */
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = jx;
                z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
                         alpha->r * x[i__2].i + alpha->i * x[i__2].r;
                temp.r = z__1.r, temp.i = z__1.i;
                iy = ky;
                i__2 = *m;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    i__3 = iy;
                    i__4 = iy;
                    i__5 = i__ + j * a_dim1;
                    z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, z__2.i =
                             temp.r * a[i__5].i + temp.i * a[i__5].r;
                    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
                    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
                    iy += *incy;
/* L70: */
                }
                jx += *incx;
/* L80: */
            }
        }
    } else {

/*        Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y. */

        jy = ky;
        if (*incx == 1) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                temp.r = 0., temp.i = 0.;
                if (noconj) {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * a_dim1;
                        i__4 = i__;
                        z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4]
                                .i, z__2.i = a[i__3].r * x[i__4].i + a[i__3]
                                .i * x[i__4].r;
                        z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                        temp.r = z__1.r, temp.i = z__1.i;
/* L90: */
                    }
                } else {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        d_lmp_cnjg(&z__3, &a[i__ + j * a_dim1]);
                        i__3 = i__;
                        z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i,
                                z__2.i = z__3.r * x[i__3].i + z__3.i * x[i__3]
                                .r;
                        z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                        temp.r = z__1.r, temp.i = z__1.i;
/* L100: */
                    }
                }
                i__2 = jy;
                i__3 = jy;
                z__2.r = alpha->r * temp.r - alpha->i * temp.i, z__2.i =
                        alpha->r * temp.i + alpha->i * temp.r;
                z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
                y[i__2].r = z__1.r, y[i__2].i = z__1.i;
                jy += *incy;
/* L110: */
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                temp.r = 0., temp.i = 0.;
                ix = kx;
                if (noconj) {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        i__3 = i__ + j * a_dim1;
                        i__4 = ix;
                        z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4]
                                .i, z__2.i = a[i__3].r * x[i__4].i + a[i__3]
                                .i * x[i__4].r;
                        z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                        temp.r = z__1.r, temp.i = z__1.i;
                        ix += *incx;
/* L120: */
                    }
                } else {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        d_lmp_cnjg(&z__3, &a[i__ + j * a_dim1]);
                        i__3 = ix;
                        z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i,
                                z__2.i = z__3.r * x[i__3].i + z__3.i * x[i__3]
                                .r;
                        z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
                        temp.r = z__1.r, temp.i = z__1.i;
                        ix += *incx;
/* L130: */
                    }
                }
                i__2 = jy;
                i__3 = jy;
                z__2.r = alpha->r * temp.r - alpha->i * temp.i, z__2.i =
                        alpha->r * temp.i + alpha->i * temp.r;
                z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
                y[i__2].r = z__1.r, y[i__2].i = z__1.i;
                jy += *incy;
/* L140: */
            }
        }
    }

    return 0;

/*     End of ZGEMV */

} /* zgemv_ */

#ifdef __cplusplus
        }
#endif
