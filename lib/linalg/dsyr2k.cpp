/* fortran/dsyr2k.f -- translated by f2c (version 20200916).
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

/* > \brief \b DSYR2K */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION ALPHA,BETA */
/*       INTEGER K,LDA,LDB,LDC,N */
/*       CHARACTER TRANS,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYR2K  performs one of the symmetric rank 2k operations */
/* > */
/* >    C := alpha*A*B**T + alpha*B*A**T + beta*C, */
/* > */
/* > or */
/* > */
/* >    C := alpha*A**T*B + alpha*B**T*A + beta*C, */
/* > */
/* > where  alpha and beta  are scalars, C is an  n by n  symmetric matrix */
/* > and  A and B  are  n by k  matrices  in the  first  case  and  k by n */
/* > matrices in the second case. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On  entry,   UPLO  specifies  whether  the  upper  or  lower */
/* >           triangular  part  of the  array  C  is to be  referenced  as */
/* >           follows: */
/* > */
/* >              UPLO = 'U' or 'u'   Only the  upper triangular part of  C */
/* >                                  is to be referenced. */
/* > */
/* >              UPLO = 'L' or 'l'   Only the  lower triangular part of  C */
/* >                                  is to be referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >           On entry,  TRANS  specifies the operation to be performed as */
/* >           follows: */
/* > */
/* >              TRANS = 'N' or 'n'   C := alpha*A*B**T + alpha*B*A**T + */
/* >                                        beta*C. */
/* > */
/* >              TRANS = 'T' or 't'   C := alpha*A**T*B + alpha*B**T*A + */
/* >                                        beta*C. */
/* > */
/* >              TRANS = 'C' or 'c'   C := alpha*A**T*B + alpha*B**T*A + */
/* >                                        beta*C. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry,  N specifies the order of the matrix C.  N must be */
/* >           at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >           On entry with  TRANS = 'N' or 'n',  K  specifies  the number */
/* >           of  columns  of the  matrices  A and B,  and on  entry  with */
/* >           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number */
/* >           of rows of the matrices  A and B.  K must be at least  zero. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is DOUBLE PRECISION. */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension ( LDA, ka ), where ka is */
/* >           k  when  TRANS = 'N' or 'n',  and is  n  otherwise. */
/* >           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k */
/* >           part of the array  A  must contain the matrix  A,  otherwise */
/* >           the leading  k by n  part of the array  A  must contain  the */
/* >           matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n' */
/* >           then  LDA must be at least  max( 1, n ), otherwise  LDA must */
/* >           be at least  max( 1, k ). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension ( LDB, kb ), where kb is */
/* >           k  when  TRANS = 'N' or 'n',  and is  n  otherwise. */
/* >           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k */
/* >           part of the array  B  must contain the matrix  B,  otherwise */
/* >           the leading  k by n  part of the array  B  must contain  the */
/* >           matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >           On entry, LDB specifies the first dimension of B as declared */
/* >           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n' */
/* >           then  LDB must be at least  max( 1, n ), otherwise  LDB must */
/* >           be at least  max( 1, k ). */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is DOUBLE PRECISION. */
/* >           On entry, BETA specifies the scalar beta. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension ( LDC, N ) */
/* >           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n */
/* >           upper triangular part of the array C must contain the upper */
/* >           triangular part  of the  symmetric matrix  and the strictly */
/* >           lower triangular part of C is not referenced.  On exit, the */
/* >           upper triangular part of the array  C is overwritten by the */
/* >           upper triangular part of the updated matrix. */
/* >           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n */
/* >           lower triangular part of the array C must contain the lower */
/* >           triangular part  of the  symmetric matrix  and the strictly */
/* >           upper triangular part of C is not referenced.  On exit, the */
/* >           lower triangular part of the array  C is overwritten by the */
/* >           lower triangular part of the updated matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >           On entry, LDC specifies the first dimension of C as declared */
/* >           in  the  calling  (sub)  program.   LDC  must  be  at  least */
/* >           max( 1, n ). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup double_blas_level3 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Level 3 Blas routine. */
/* > */
/* > */
/* >  -- Written on 8-February-1989. */
/* >     Jack Dongarra, Argonne National Laboratory. */
/* >     Iain Duff, AERE Harwell. */
/* >     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/* >     Sven Hammarling, Numerical Algorithms Group Ltd. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dsyr2k_(char *uplo, char *trans, integer *n, integer *k,
        doublereal *alpha, doublereal *a, integer *lda, doublereal *b,
        integer *ldb, doublereal *beta, doublereal *c__, integer *ldc, ftnlen
        uplo_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2,
            i__3;

    /* Local variables */
    integer i__, j, l, info;
    doublereal temp1, temp2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer nrowa;
    logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- Reference BLAS level3 routine -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Parameters .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */
    if (lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1)) {
        nrowa = *n;
    } else {
        nrowa = *k;
    }
    upper = lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1);

    info = 0;
    if (! upper && ! lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1)) {
        info = 1;
    } else if (! lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans,
            (char *)"T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, (char *)"C", (ftnlen)1, (
            ftnlen)1)) {
        info = 2;
    } else if (*n < 0) {
        info = 3;
    } else if (*k < 0) {
        info = 4;
    } else if (*lda < max(1,nrowa)) {
        info = 7;
    } else if (*ldb < max(1,nrowa)) {
        info = 9;
    } else if (*ldc < max(1,*n)) {
        info = 12;
    }
    if (info != 0) {
        xerbla_((char *)"DSYR2K", &info, (ftnlen)6);
        return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
        return 0;
    }

/*     And when  alpha.eq.zero. */

    if (*alpha == 0.) {
        if (upper) {
            if (*beta == 0.) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = 0.;
/* L10: */
                    }
/* L20: */
                }
            } else {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L30: */
                    }
/* L40: */
                }
            }
        } else {
            if (*beta == 0.) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = 0.;
/* L50: */
                    }
/* L60: */
                }
            } else {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L70: */
                    }
/* L80: */
                }
            }
        }
        return 0;
    }

/*     Start the operations. */

    if (lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1)) {

/*        Form  C := alpha*A*B**T + alpha*B*A**T + C. */

        if (upper) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (*beta == 0.) {
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = 0.;
/* L90: */
                    }
                } else if (*beta != 1.) {
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L100: */
                    }
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l) {
                    if (a[j + l * a_dim1] != 0. || b[j + l * b_dim1] != 0.) {
                        temp1 = *alpha * b[j + l * b_dim1];
                        temp2 = *alpha * a[j + l * a_dim1];
                        i__3 = j;
                        for (i__ = 1; i__ <= i__3; ++i__) {
                            c__[i__ + j * c_dim1] = c__[i__ + j * c_dim1] + a[
                                    i__ + l * a_dim1] * temp1 + b[i__ + l *
                                    b_dim1] * temp2;
/* L110: */
                        }
                    }
/* L120: */
                }
/* L130: */
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                if (*beta == 0.) {
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = 0.;
/* L140: */
                    }
                } else if (*beta != 1.) {
                    i__2 = *n;
                    for (i__ = j; i__ <= i__2; ++i__) {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
/* L150: */
                    }
                }
                i__2 = *k;
                for (l = 1; l <= i__2; ++l) {
                    if (a[j + l * a_dim1] != 0. || b[j + l * b_dim1] != 0.) {
                        temp1 = *alpha * b[j + l * b_dim1];
                        temp2 = *alpha * a[j + l * a_dim1];
                        i__3 = *n;
                        for (i__ = j; i__ <= i__3; ++i__) {
                            c__[i__ + j * c_dim1] = c__[i__ + j * c_dim1] + a[
                                    i__ + l * a_dim1] * temp1 + b[i__ + l *
                                    b_dim1] * temp2;
/* L160: */
                        }
                    }
/* L170: */
                }
/* L180: */
            }
        }
    } else {

/*        Form  C := alpha*A**T*B + alpha*B**T*A + C. */

        if (upper) {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = j;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp1 = 0.;
                    temp2 = 0.;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l) {
                        temp1 += a[l + i__ * a_dim1] * b[l + j * b_dim1];
                        temp2 += b[l + i__ * b_dim1] * a[l + j * a_dim1];
/* L190: */
                    }
                    if (*beta == 0.) {
                        c__[i__ + j * c_dim1] = *alpha * temp1 + *alpha *
                                temp2;
                    } else {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1]
                                + *alpha * temp1 + *alpha * temp2;
                    }
/* L200: */
                }
/* L210: */
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= i__1; ++j) {
                i__2 = *n;
                for (i__ = j; i__ <= i__2; ++i__) {
                    temp1 = 0.;
                    temp2 = 0.;
                    i__3 = *k;
                    for (l = 1; l <= i__3; ++l) {
                        temp1 += a[l + i__ * a_dim1] * b[l + j * b_dim1];
                        temp2 += b[l + i__ * b_dim1] * a[l + j * a_dim1];
/* L220: */
                    }
                    if (*beta == 0.) {
                        c__[i__ + j * c_dim1] = *alpha * temp1 + *alpha *
                                temp2;
                    } else {
                        c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1]
                                + *alpha * temp1 + *alpha * temp2;
                    }
/* L230: */
                }
/* L240: */
            }
        }
    }

    return 0;

/*     End of DSYR2K */

} /* dsyr2k_ */

#ifdef __cplusplus
        }
#endif
