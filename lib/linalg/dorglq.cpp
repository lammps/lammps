/* fortran/dorglq.f -- translated by f2c (version 20200916).
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
static integer c__3 = 3;
static integer c__2 = 2;

/* > \brief \b DORGLQ */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DORGLQ + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorglq.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorglq.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorglq.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
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
/* > DORGLQ generates an M-by-N real matrix Q with orthonormal rows, */
/* > which is defined as the first M rows of a product of K elementary */
/* > reflectors of order N */
/* > */
/* >       Q  =  H(k) . . . H(2) H(1) */
/* > */
/* > as returned by DGELQF. */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the i-th row must contain the vector which defines */
/* >          the elementary reflector H(i), for i = 1,2,...,k, as returned */
/* >          by DGELQF in the first k rows of its array argument A. */
/* >          On exit, the M-by-N matrix Q. */
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
/* >          TAU is DOUBLE PRECISION array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by DGELQF. */
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
/* >          The dimension of the array WORK. LWORK >= max(1,M). */
/* >          For optimum performance LWORK >= M*NB, where NB is */
/* >          the optimal blocksize. */
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
/* >          < 0:  if INFO = -i, the i-th argument has an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dorglq_(integer *m, integer *n, integer *k, doublereal *
        a, integer *lda, doublereal *tau, doublereal *work, integer *lwork,
        integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    integer i__, j, l, ib, nb, ki, kk, nx, iws, nbmin, iinfo;
    extern /* Subroutine */ int dorgl2_(integer *, integer *, integer *,
            doublereal *, integer *, doublereal *, doublereal *, integer *),
            dlarfb_(char *, char *, char *, char *, integer *, integer *,
            integer *, doublereal *, integer *, doublereal *, integer *,
            doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen,
            ftnlen, ftnlen), dlarft_(char *, char *, integer *, integer *,
            doublereal *, integer *, doublereal *, doublereal *, integer *,
            ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *,
            integer *, integer *, ftnlen, ftnlen);
    integer ldwork, lwkopt;
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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
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
    nb = ilaenv_(&c__1, (char *)"DORGLQ", (char *)" ", m, n, k, &c_n1, (ftnlen)6, (ftnlen)1);
    lwkopt = max(1,*m) * nb;
    work[1] = (doublereal) lwkopt;
    lquery = *lwork == -1;
    if (*m < 0) {
        *info = -1;
    } else if (*n < *m) {
        *info = -2;
    } else if (*k < 0 || *k > *m) {
        *info = -3;
    } else if (*lda < max(1,*m)) {
        *info = -5;
    } else if (*lwork < max(1,*m) && ! lquery) {
        *info = -8;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DORGLQ", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }

/*     Quick return if possible */

    if (*m <= 0) {
        work[1] = 1.;
        return 0;
    }

    nbmin = 2;
    nx = 0;
    iws = *m;
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
        i__1 = 0, i__2 = ilaenv_(&c__3, (char *)"DORGLQ", (char *)" ", m, n, k, &c_n1, (
                ftnlen)6, (ftnlen)1);
        nx = max(i__1,i__2);
        if (nx < *k) {

/*           Determine if workspace is large enough for blocked code. */

            ldwork = *m;
            iws = ldwork * nb;
            if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

                nb = *lwork / ldwork;
/* Computing MAX */
                i__1 = 2, i__2 = ilaenv_(&c__2, (char *)"DORGLQ", (char *)" ", m, n, k, &c_n1,
                         (ftnlen)6, (ftnlen)1);
                nbmin = max(i__1,i__2);
            }
        }
    }

    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the last block. */
/*        The first kk rows are handled by the block method. */

        ki = (*k - nx - 1) / nb * nb;
/* Computing MIN */
        i__1 = *k, i__2 = ki + nb;
        kk = min(i__1,i__2);

/*        Set A(kk+1:m,1:kk) to zero. */

        i__1 = kk;
        for (j = 1; j <= i__1; ++j) {
            i__2 = *m;
            for (i__ = kk + 1; i__ <= i__2; ++i__) {
                a[i__ + j * a_dim1] = 0.;
/* L10: */
            }
/* L20: */
        }
    } else {
        kk = 0;
    }

/*     Use unblocked code for the last or only block. */

    if (kk < *m) {
        i__1 = *m - kk;
        i__2 = *n - kk;
        i__3 = *k - kk;
        dorgl2_(&i__1, &i__2, &i__3, &a[kk + 1 + (kk + 1) * a_dim1], lda, &
                tau[kk + 1], &work[1], &iinfo);
    }

    if (kk > 0) {

/*        Use blocked code */

        i__1 = -nb;
        for (i__ = ki + 1; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
            i__2 = nb, i__3 = *k - i__ + 1;
            ib = min(i__2,i__3);
            if (i__ + ib <= *m) {

/*              Form the triangular factor of the block reflector */
/*              H = H(i) H(i+1) . . . H(i+ib-1) */

                i__2 = *n - i__ + 1;
                dlarft_((char *)"Forward", (char *)"Rowwise", &i__2, &ib, &a[i__ + i__ *
                        a_dim1], lda, &tau[i__], &work[1], &ldwork, (ftnlen)7,
                         (ftnlen)7);

/*              Apply H**T to A(i+ib:m,i:n) from the right */

                i__2 = *m - i__ - ib + 1;
                i__3 = *n - i__ + 1;
                dlarfb_((char *)"Right", (char *)"Transpose", (char *)"Forward", (char *)"Rowwise", &i__2, &
                        i__3, &ib, &a[i__ + i__ * a_dim1], lda, &work[1], &
                        ldwork, &a[i__ + ib + i__ * a_dim1], lda, &work[ib +
                        1], &ldwork, (ftnlen)5, (ftnlen)9, (ftnlen)7, (ftnlen)
                        7);
            }

/*           Apply H**T to columns i:n of current block */

            i__2 = *n - i__ + 1;
            dorgl2_(&ib, &i__2, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &
                    work[1], &iinfo);

/*           Set columns 1:i-1 of current block to zero */

            i__2 = i__ - 1;
            for (j = 1; j <= i__2; ++j) {
                i__3 = i__ + ib - 1;
                for (l = i__; l <= i__3; ++l) {
                    a[l + j * a_dim1] = 0.;
/* L30: */
                }
/* L40: */
            }
/* L50: */
        }
    }

    work[1] = (doublereal) iws;
    return 0;

/*     End of DORGLQ */

} /* dorglq_ */

#ifdef __cplusplus
        }
#endif
