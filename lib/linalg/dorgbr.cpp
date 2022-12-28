/* fortran/dorgbr.f -- translated by f2c (version 20200916).
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

static integer c_n1 = -1;

/* > \brief \b DORGBR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DORGBR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgbr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgbr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgbr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          VECT */
/*       INTEGER            INFO, K, LDA, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DORGBR generates one of the real orthogonal matrices Q or P**T */
/* > determined by DGEBRD when reducing a real matrix A to bidiagonal */
/* > form: A = Q * B * P**T.  Q and P**T are defined as products of */
/* > elementary reflectors H(i) or G(i) respectively. */
/* > */
/* > If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q */
/* > is of order M: */
/* > if m >= k, Q = H(1) H(2) . . . H(k) and DORGBR returns the first n */
/* > columns of Q, where m >= n >= k; */
/* > if m < k, Q = H(1) H(2) . . . H(m-1) and DORGBR returns Q as an */
/* > M-by-M matrix. */
/* > */
/* > If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**T */
/* > is of order N: */
/* > if k < n, P**T = G(k) . . . G(2) G(1) and DORGBR returns the first m */
/* > rows of P**T, where n >= m >= k; */
/* > if k >= n, P**T = G(n-1) . . . G(2) G(1) and DORGBR returns P**T as */
/* > an N-by-N matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] VECT */
/* > \verbatim */
/* >          VECT is CHARACTER*1 */
/* >          Specifies whether the matrix Q or the matrix P**T is */
/* >          required, as defined in the transformation applied by DGEBRD: */
/* >          = 'Q':  generate Q; */
/* >          = 'P':  generate P**T. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix Q or P**T to be returned. */
/* >          M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix Q or P**T to be returned. */
/* >          N >= 0. */
/* >          If VECT = 'Q', M >= N >= min(M,K); */
/* >          if VECT = 'P', N >= M >= min(N,K). */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          If VECT = 'Q', the number of columns in the original M-by-K */
/* >          matrix reduced by DGEBRD. */
/* >          If VECT = 'P', the number of rows in the original K-by-N */
/* >          matrix reduced by DGEBRD. */
/* >          K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the vectors which define the elementary reflectors, */
/* >          as returned by DGEBRD. */
/* >          On exit, the M-by-N matrix Q or P**T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is DOUBLE PRECISION array, dimension */
/* >                                (min(M,K)) if VECT = 'Q' */
/* >                                (min(N,K)) if VECT = 'P' */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i) or G(i), which determines Q or P**T, as */
/* >          returned by DGEBRD in its array argument TAUQ or TAUP. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. LWORK >= max(1,min(M,N)). */
/* >          For optimum performance LWORK >= min(M,N)*NB, where NB */
/* >          is the optimal blocksize. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
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

/* > \ingroup doubleGBcomputational */

/*  ===================================================================== */
/* Subroutine */ int dorgbr_(char *vect, integer *m, integer *n, integer *k, 
	doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	integer *lwork, integer *info, ftnlen vect_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    integer i__, j, mn;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer iinfo;
    logical wantq;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dorglq_(
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *), dorgqr_(
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *);
    integer lwkopt;
    logical lquery;


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

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    wantq = lsame_(vect, (char *)"Q", (ftnlen)1, (ftnlen)1);
    mn = min(*m,*n);
    lquery = *lwork == -1;
    if (! wantq && ! lsame_(vect, (char *)"P", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*n < 0 || wantq && (*n > *m || *n < min(*m,*k)) || ! wantq && (
	    *m > *n || *m < min(*n,*k))) {
	*info = -3;
    } else if (*k < 0) {
	*info = -4;
    } else if (*lda < max(1,*m)) {
	*info = -6;
    } else if (*lwork < max(1,mn) && ! lquery) {
	*info = -9;
    }

    if (*info == 0) {
	work[1] = 1.;
	if (wantq) {
	    if (*m >= *k) {
		dorgqr_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, 
			&iinfo);
	    } else {
		if (*m > 1) {
		    i__1 = *m - 1;
		    i__2 = *m - 1;
		    i__3 = *m - 1;
		    dorgqr_(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], &
			    work[1], &c_n1, &iinfo);
		}
	    }
	} else {
	    if (*k < *n) {
		dorglq_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, 
			&iinfo);
	    } else {
		if (*n > 1) {
		    i__1 = *n - 1;
		    i__2 = *n - 1;
		    i__3 = *n - 1;
		    dorglq_(&i__1, &i__2, &i__3, &a[a_offset], lda, &tau[1], &
			    work[1], &c_n1, &iinfo);
		}
	    }
	}
	lwkopt = (integer) work[1];
	lwkopt = max(lwkopt,mn);
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_((char *)"DORGBR", &i__1, (ftnlen)6);
	return 0;
    } else if (lquery) {
	work[1] = (doublereal) lwkopt;
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	work[1] = 1.;
	return 0;
    }

    if (wantq) {

/*        Form Q, determined by a call to DGEBRD to reduce an m-by-k */
/*        matrix */

	if (*m >= *k) {

/*           If m >= k, assume m >= n >= k */

	    dorgqr_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], lwork, &
		    iinfo);

	} else {

/*           If m < k, assume m = n */

/*           Shift the vectors which define the elementary reflectors one */
/*           column to the right, and set the first row and column of Q */
/*           to those of the unit matrix */

	    for (j = *m; j >= 2; --j) {
		a[j * a_dim1 + 1] = 0.;
		i__1 = *m;
		for (i__ = j + 1; i__ <= i__1; ++i__) {
		    a[i__ + j * a_dim1] = a[i__ + (j - 1) * a_dim1];
/* L10: */
		}
/* L20: */
	    }
	    a[a_dim1 + 1] = 1.;
	    i__1 = *m;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		a[i__ + a_dim1] = 0.;
/* L30: */
	    }
	    if (*m > 1) {

/*              Form Q(2:m,2:m) */

		i__1 = *m - 1;
		i__2 = *m - 1;
		i__3 = *m - 1;
		dorgqr_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[
			1], &work[1], lwork, &iinfo);
	    }
	}
    } else {

/*        Form P**T, determined by a call to DGEBRD to reduce a k-by-n */
/*        matrix */

	if (*k < *n) {

/*           If k < n, assume k <= m <= n */

	    dorglq_(m, n, k, &a[a_offset], lda, &tau[1], &work[1], lwork, &
		    iinfo);

	} else {

/*           If k >= n, assume m = n */

/*           Shift the vectors which define the elementary reflectors one */
/*           row downward, and set the first row and column of P**T to */
/*           those of the unit matrix */

	    a[a_dim1 + 1] = 1.;
	    i__1 = *n;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		a[i__ + a_dim1] = 0.;
/* L40: */
	    }
	    i__1 = *n;
	    for (j = 2; j <= i__1; ++j) {
		for (i__ = j - 1; i__ >= 2; --i__) {
		    a[i__ + j * a_dim1] = a[i__ - 1 + j * a_dim1];
/* L50: */
		}
		a[j * a_dim1 + 1] = 0.;
/* L60: */
	    }
	    if (*n > 1) {

/*              Form P**T(2:n,2:n) */

		i__1 = *n - 1;
		i__2 = *n - 1;
		i__3 = *n - 1;
		dorglq_(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[
			1], &work[1], lwork, &iinfo);
	    }
	}
    }
    work[1] = (doublereal) lwkopt;
    return 0;

/*     End of DORGBR */

} /* dorgbr_ */

#ifdef __cplusplus
	}
#endif
