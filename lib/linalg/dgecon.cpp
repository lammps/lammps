/* fortran/dgecon.f -- translated by f2c (version 20200916).
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

/* > \brief \b DGECON */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGECON + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgecon.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgecon.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgecon.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          NORM */
/*       INTEGER            INFO, LDA, N */
/*       DOUBLE PRECISION   ANORM, RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGECON estimates the reciprocal of the condition number of a general */
/* > real matrix A, in either the 1-norm or the infinity-norm, using */
/* > the LU factorization computed by DGETRF. */
/* > */
/* > An estimate is obtained for norm(inv(A)), and the reciprocal of the */
/* > condition number is computed as */
/* >    RCOND = 1 / ( norm(A) * norm(inv(A)) ). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] NORM */
/* > \verbatim */
/* >          NORM is CHARACTER*1 */
/* >          Specifies whether the 1-norm condition number or the */
/* >          infinity-norm condition number is required: */
/* >          = '1' or 'O':  1-norm; */
/* >          = 'I':         Infinity-norm. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          The factors L and U from the factorization A = P*L*U */
/* >          as computed by DGETRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] ANORM */
/* > \verbatim */
/* >          ANORM is DOUBLE PRECISION */
/* >          If NORM = '1' or 'O', the 1-norm of the original matrix A. */
/* >          If NORM = 'I', the infinity-norm of the original matrix A. */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* >          RCOND is DOUBLE PRECISION */
/* >          The reciprocal of the condition number of the matrix A, */
/* >          computed as RCOND = 1/(norm(A) * norm(inv(A))). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (4*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doubleGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int dgecon_(char *norm, integer *n, doublereal *a, integer *
	lda, doublereal *anorm, doublereal *rcond, doublereal *work, integer *
	iwork, integer *info, ftnlen norm_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1;

    /* Local variables */
    doublereal sl;
    integer ix;
    doublereal su;
    integer kase, kase1;
    doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer isave[3];
    extern /* Subroutine */ int drscl_(integer *, doublereal *, doublereal *, 
	    integer *), dlacn2_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    doublereal ainvnm;
    extern /* Subroutine */ int dlatrs_(char *, char *, char *, char *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    logical onenrm;
    char normin[1];
    doublereal smlnum;


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
/*     .. Local Arrays .. */
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
    --work;
    --iwork;

    /* Function Body */
    *info = 0;
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, (char *)"O", (ftnlen)1, (
	    ftnlen)1);
    if (! onenrm && ! lsame_(norm, (char *)"I", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    } else if (*anorm < 0.) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_((char *)"DGECON", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible */

    *rcond = 0.;
    if (*n == 0) {
	*rcond = 1.;
	return 0;
    } else if (*anorm == 0.) {
	return 0;
    }

    smlnum = dlamch_((char *)"Safe minimum", (ftnlen)12);

/*     Estimate the norm of inv(A). */

    ainvnm = 0.;
    *(unsigned char *)normin = 'N';
    if (onenrm) {
	kase1 = 1;
    } else {
	kase1 = 2;
    }
    kase = 0;
L10:
    dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
    if (kase != 0) {
	if (kase == kase1) {

/*           Multiply by inv(L). */

	    dlatrs_((char *)"Lower", (char *)"No transpose", (char *)"Unit", normin, n, &a[a_offset], 
		    lda, &work[1], &sl, &work[(*n << 1) + 1], info, (ftnlen)5,
		     (ftnlen)12, (ftnlen)4, (ftnlen)1);

/*           Multiply by inv(U). */

	    dlatrs_((char *)"Upper", (char *)"No transpose", (char *)"Non-unit", normin, n, &a[
		    a_offset], lda, &work[1], &su, &work[*n * 3 + 1], info, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
	} else {

/*           Multiply by inv(U**T). */

	    dlatrs_((char *)"Upper", (char *)"Transpose", (char *)"Non-unit", normin, n, &a[a_offset],
		     lda, &work[1], &su, &work[*n * 3 + 1], info, (ftnlen)5, (
		    ftnlen)9, (ftnlen)8, (ftnlen)1);

/*           Multiply by inv(L**T). */

	    dlatrs_((char *)"Lower", (char *)"Transpose", (char *)"Unit", normin, n, &a[a_offset], 
		    lda, &work[1], &sl, &work[(*n << 1) + 1], info, (ftnlen)5,
		     (ftnlen)9, (ftnlen)4, (ftnlen)1);
	}

/*        Divide X by 1/(SL*SU) if doing so will not cause overflow. */

	scale = sl * su;
	*(unsigned char *)normin = 'Y';
	if (scale != 1.) {
	    ix = idamax_(n, &work[1], &c__1);
	    if (scale < (d__1 = work[ix], abs(d__1)) * smlnum || scale == 0.) 
		    {
		goto L20;
	    }
	    drscl_(n, &scale, &work[1], &c__1);
	}
	goto L10;
    }

/*     Compute the estimate of the reciprocal condition number. */

    if (ainvnm != 0.) {
	*rcond = 1. / ainvnm / *anorm;
    }

L20:
    return 0;

/*     End of DGECON */

} /* dgecon_ */

#ifdef __cplusplus
	}
#endif
