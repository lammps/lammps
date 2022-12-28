/* fortran/dtrmv.f -- translated by f2c (version 20200916).
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

/* > \brief \b DTRMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,LDA,N */
/*       CHARACTER DIAG,TRANS,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION A(LDA,*),X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTRMV  performs one of the matrix-vector operations */
/* > */
/* >    x := A*x,   or   x := A**T*x, */
/* > */
/* > where x is an n element vector and  A is an n by n unit, or non-unit, */
/* > upper or lower triangular matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On entry, UPLO specifies whether the matrix is an upper or */
/* >           lower triangular matrix as follows: */
/* > */
/* >              UPLO = 'U' or 'u'   A is an upper triangular matrix. */
/* > */
/* >              UPLO = 'L' or 'l'   A is a lower triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >           On entry, TRANS specifies the operation to be performed as */
/* >           follows: */
/* > */
/* >              TRANS = 'N' or 'n'   x := A*x. */
/* > */
/* >              TRANS = 'T' or 't'   x := A**T*x. */
/* > */
/* >              TRANS = 'C' or 'c'   x := A**T*x. */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >           On entry, DIAG specifies whether or not A is unit */
/* >           triangular as follows: */
/* > */
/* >              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */
/* > */
/* >              DIAG = 'N' or 'n'   A is not assumed to be unit */
/* >                                  triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the order of the matrix A. */
/* >           N must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension ( LDA, N ) */
/* >           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/* >           upper triangular part of the array A must contain the upper */
/* >           triangular matrix and the strictly lower triangular part of */
/* >           A is not referenced. */
/* >           Before entry with UPLO = 'L' or 'l', the leading n by n */
/* >           lower triangular part of the array A must contain the lower */
/* >           triangular matrix and the strictly upper triangular part of */
/* >           A is not referenced. */
/* >           Note that when  DIAG = 'U' or 'u', the diagonal elements of */
/* >           A are not referenced either, but are assumed to be unity. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           max( 1, n ). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCX ) ). */
/* >           Before entry, the incremented array X must contain the n */
/* >           element vector x. On exit, X is overwritten with the */
/* >           transformed vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >           On entry, INCX specifies the increment for the elements of */
/* >           X. INCX must not be zero. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup double_blas_level2 */

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
/* Subroutine */ int dtrmv_(char *uplo, char *trans, char *diag, integer *n,
        doublereal *a, integer *lda, doublereal *x, integer *incx, ftnlen
        uplo_len, ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    integer i__, j, ix, jx, kx, info;
    doublereal temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    logical nounit;


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

    /* Function Body */
    info = 0;
    if (! lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, (char *)"L", (
            ftnlen)1, (ftnlen)1)) {
        info = 1;
    } else if (! lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans,
            (char *)"T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, (char *)"C", (ftnlen)1, (
            ftnlen)1)) {
        info = 2;
    } else if (! lsame_(diag, (char *)"U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag,
            (char *)"N", (ftnlen)1, (ftnlen)1)) {
        info = 3;
    } else if (*n < 0) {
        info = 4;
    } else if (*lda < max(1,*n)) {
        info = 6;
    } else if (*incx == 0) {
        info = 8;
    }
    if (info != 0) {
        xerbla_((char *)"DTRMV ", &info, (ftnlen)6);
        return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
        return 0;
    }

    nounit = lsame_(diag, (char *)"N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
        kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
        kx = 1;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

    if (lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x := A*x. */

        if (lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1)) {
            if (*incx == 1) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    if (x[j] != 0.) {
                        temp = x[j];
                        i__2 = j - 1;
                        for (i__ = 1; i__ <= i__2; ++i__) {
                            x[i__] += temp * a[i__ + j * a_dim1];
/* L10: */
                        }
                        if (nounit) {
                            x[j] *= a[j + j * a_dim1];
                        }
                    }
/* L20: */
                }
            } else {
                jx = kx;
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    if (x[jx] != 0.) {
                        temp = x[jx];
                        ix = kx;
                        i__2 = j - 1;
                        for (i__ = 1; i__ <= i__2; ++i__) {
                            x[ix] += temp * a[i__ + j * a_dim1];
                            ix += *incx;
/* L30: */
                        }
                        if (nounit) {
                            x[jx] *= a[j + j * a_dim1];
                        }
                    }
                    jx += *incx;
/* L40: */
                }
            }
        } else {
            if (*incx == 1) {
                for (j = *n; j >= 1; --j) {
                    if (x[j] != 0.) {
                        temp = x[j];
                        i__1 = j + 1;
                        for (i__ = *n; i__ >= i__1; --i__) {
                            x[i__] += temp * a[i__ + j * a_dim1];
/* L50: */
                        }
                        if (nounit) {
                            x[j] *= a[j + j * a_dim1];
                        }
                    }
/* L60: */
                }
            } else {
                kx += (*n - 1) * *incx;
                jx = kx;
                for (j = *n; j >= 1; --j) {
                    if (x[jx] != 0.) {
                        temp = x[jx];
                        ix = kx;
                        i__1 = j + 1;
                        for (i__ = *n; i__ >= i__1; --i__) {
                            x[ix] += temp * a[i__ + j * a_dim1];
                            ix -= *incx;
/* L70: */
                        }
                        if (nounit) {
                            x[jx] *= a[j + j * a_dim1];
                        }
                    }
                    jx -= *incx;
/* L80: */
                }
            }
        }
    } else {

/*        Form  x := A**T*x. */

        if (lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1)) {
            if (*incx == 1) {
                for (j = *n; j >= 1; --j) {
                    temp = x[j];
                    if (nounit) {
                        temp *= a[j + j * a_dim1];
                    }
                    for (i__ = j - 1; i__ >= 1; --i__) {
                        temp += a[i__ + j * a_dim1] * x[i__];
/* L90: */
                    }
                    x[j] = temp;
/* L100: */
                }
            } else {
                jx = kx + (*n - 1) * *incx;
                for (j = *n; j >= 1; --j) {
                    temp = x[jx];
                    ix = jx;
                    if (nounit) {
                        temp *= a[j + j * a_dim1];
                    }
                    for (i__ = j - 1; i__ >= 1; --i__) {
                        ix -= *incx;
                        temp += a[i__ + j * a_dim1] * x[ix];
/* L110: */
                    }
                    x[jx] = temp;
                    jx -= *incx;
/* L120: */
                }
            }
        } else {
            if (*incx == 1) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    temp = x[j];
                    if (nounit) {
                        temp *= a[j + j * a_dim1];
                    }
                    i__2 = *n;
                    for (i__ = j + 1; i__ <= i__2; ++i__) {
                        temp += a[i__ + j * a_dim1] * x[i__];
/* L130: */
                    }
                    x[j] = temp;
/* L140: */
                }
            } else {
                jx = kx;
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    temp = x[jx];
                    ix = jx;
                    if (nounit) {
                        temp *= a[j + j * a_dim1];
                    }
                    i__2 = *n;
                    for (i__ = j + 1; i__ <= i__2; ++i__) {
                        ix += *incx;
                        temp += a[i__ + j * a_dim1] * x[ix];
/* L150: */
                    }
                    x[jx] = temp;
                    jx += *incx;
/* L160: */
                }
            }
        }
    }

    return 0;

/*     End of DTRMV */

} /* dtrmv_ */

#ifdef __cplusplus
        }
#endif
