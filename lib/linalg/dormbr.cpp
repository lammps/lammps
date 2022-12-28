/* fortran/dormbr.f -- translated by f2c (version 20200916).
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

/* > \brief \b DORMBR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DORMBR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormbr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormbr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormbr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DORMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C, */
/*                          LDC, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS, VECT */
/*       INTEGER            INFO, K, LDA, LDC, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > If VECT = 'Q', DORMBR overwrites the general real M-by-N matrix C */
/* > with */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'T':      Q**T * C       C * Q**T */
/* > */
/* > If VECT = 'P', DORMBR overwrites the general real M-by-N matrix C */
/* > with */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      P * C          C * P */
/* > TRANS = 'T':      P**T * C       C * P**T */
/* > */
/* > Here Q and P**T are the orthogonal matrices determined by DGEBRD when */
/* > reducing a real matrix A to bidiagonal form: A = Q * B * P**T. Q and */
/* > P**T are defined as products of elementary reflectors H(i) and G(i) */
/* > respectively. */
/* > */
/* > Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the */
/* > order of the orthogonal matrix Q or P**T that is applied. */
/* > */
/* > If VECT = 'Q', A is assumed to have been an NQ-by-K matrix: */
/* > if nq >= k, Q = H(1) H(2) . . . H(k); */
/* > if nq < k, Q = H(1) H(2) . . . H(nq-1). */
/* > */
/* > If VECT = 'P', A is assumed to have been a K-by-NQ matrix: */
/* > if k < nq, P = G(1) G(2) . . . G(k); */
/* > if k >= nq, P = G(1) G(2) . . . G(nq-1). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] VECT */
/* > \verbatim */
/* >          VECT is CHARACTER*1 */
/* >          = 'Q': apply Q or Q**T; */
/* >          = 'P': apply P or P**T. */
/* > \endverbatim */
/* > */
/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q, Q**T, P or P**T from the Left; */
/* >          = 'R': apply Q, Q**T, P or P**T from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  No transpose, apply Q  or P; */
/* >          = 'T':  Transpose, apply Q**T or P**T. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix C. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix C. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          If VECT = 'Q', the number of columns in the original */
/* >          matrix reduced by DGEBRD. */
/* >          If VECT = 'P', the number of rows in the original */
/* >          matrix reduced by DGEBRD. */
/* >          K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension */
/* >                                (LDA,min(nq,K)) if VECT = 'Q' */
/* >                                (LDA,nq)        if VECT = 'P' */
/* >          The vectors which define the elementary reflectors H(i) and */
/* >          G(i), whose products determine the matrices Q and P, as */
/* >          returned by DGEBRD. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. */
/* >          If VECT = 'Q', LDA >= max(1,nq); */
/* >          if VECT = 'P', LDA >= max(1,min(nq,K)). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is DOUBLE PRECISION array, dimension (min(nq,K)) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i) or G(i) which determines Q or P, as returned */
/* >          by DGEBRD in the array argument TAUQ or TAUP. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
/* >          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q */
/* >          or P*C or P**T*C or C*P or C*P**T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >          The leading dimension of the array C. LDC >= max(1,M). */
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
/* >          The dimension of the array WORK. */
/* >          If SIDE = 'L', LWORK >= max(1,N); */
/* >          if SIDE = 'R', LWORK >= max(1,M). */
/* >          For optimum performance LWORK >= N*NB if SIDE = 'L', and */
/* >          LWORK >= M*NB if SIDE = 'R', where NB is the optimal */
/* >          blocksize. */
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

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dormbr_(char *vect, char *side, char *trans, integer *m,
        integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau,
        doublereal *c__, integer *ldc, doublereal *work, integer *lwork,
        integer *info, ftnlen vect_len, ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3[2];
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_lmp_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    integer i1, i2, nb, mi, ni, nq, nw;
    logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer iinfo;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *,
            integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dormlq_(char *, char *, integer *, integer *,
            integer *, doublereal *, integer *, doublereal *, doublereal *,
            integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    logical notran;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *,
            integer *, doublereal *, integer *, doublereal *, doublereal *,
            integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    logical applyq;
    char transt[1];
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
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    /* Function Body */
    *info = 0;
    applyq = lsame_(vect, (char *)"Q", (ftnlen)1, (ftnlen)1);
    left = lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1);
    notran = lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1);
    lquery = *lwork == -1;

/*     NQ is the order of Q or P and NW is the minimum dimension of WORK */

    if (left) {
        nq = *m;
        nw = max(1,*n);
    } else {
        nq = *n;
        nw = max(1,*m);
    }
    if (! applyq && ! lsame_(vect, (char *)"P", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (! left && ! lsame_(side, (char *)"R", (ftnlen)1, (ftnlen)1)) {
        *info = -2;
    } else if (! notran && ! lsame_(trans, (char *)"T", (ftnlen)1, (ftnlen)1)) {
        *info = -3;
    } else if (*m < 0) {
        *info = -4;
    } else if (*n < 0) {
        *info = -5;
    } else if (*k < 0) {
        *info = -6;
    } else /* if(complicated condition) */ {
/* Computing MAX */
        i__1 = 1, i__2 = min(nq,*k);
        if (applyq && *lda < max(1,nq) || ! applyq && *lda < max(i__1,i__2)) {
            *info = -8;
        } else if (*ldc < max(1,*m)) {
            *info = -11;
        } else if (*lwork < nw && ! lquery) {
            *info = -13;
        }
    }

    if (*info == 0) {
        if (applyq) {
            if (left) {
/* Writing concatenation */
                i__3[0] = 1, a__1[0] = side;
                i__3[1] = 1, a__1[1] = trans;
                s_lmp_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
                i__1 = *m - 1;
                i__2 = *m - 1;
                nb = ilaenv_(&c__1, (char *)"DORMQR", ch__1, &i__1, n, &i__2, &c_n1, (
                        ftnlen)6, (ftnlen)2);
            } else {
/* Writing concatenation */
                i__3[0] = 1, a__1[0] = side;
                i__3[1] = 1, a__1[1] = trans;
                s_lmp_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
                i__1 = *n - 1;
                i__2 = *n - 1;
                nb = ilaenv_(&c__1, (char *)"DORMQR", ch__1, m, &i__1, &i__2, &c_n1, (
                        ftnlen)6, (ftnlen)2);
            }
        } else {
            if (left) {
/* Writing concatenation */
                i__3[0] = 1, a__1[0] = side;
                i__3[1] = 1, a__1[1] = trans;
                s_lmp_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
                i__1 = *m - 1;
                i__2 = *m - 1;
                nb = ilaenv_(&c__1, (char *)"DORMLQ", ch__1, &i__1, n, &i__2, &c_n1, (
                        ftnlen)6, (ftnlen)2);
            } else {
/* Writing concatenation */
                i__3[0] = 1, a__1[0] = side;
                i__3[1] = 1, a__1[1] = trans;
                s_lmp_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
                i__1 = *n - 1;
                i__2 = *n - 1;
                nb = ilaenv_(&c__1, (char *)"DORMLQ", ch__1, m, &i__1, &i__2, &c_n1, (
                        ftnlen)6, (ftnlen)2);
            }
        }
        lwkopt = nw * nb;
        work[1] = (doublereal) lwkopt;
    }

    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DORMBR", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }

/*     Quick return if possible */

    work[1] = 1.;
    if (*m == 0 || *n == 0) {
        return 0;
    }

    if (applyq) {

/*        Apply Q */

        if (nq >= *k) {

/*           Q was determined by a call to DGEBRD with nq >= k */

            dormqr_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
                    c_offset], ldc, &work[1], lwork, &iinfo, (ftnlen)1, (
                    ftnlen)1);
        } else if (nq > 1) {

/*           Q was determined by a call to DGEBRD with nq < k */

            if (left) {
                mi = *m - 1;
                ni = *n;
                i1 = 2;
                i2 = 1;
            } else {
                mi = *m;
                ni = *n - 1;
                i1 = 1;
                i2 = 2;
            }
            i__1 = nq - 1;
            dormqr_(side, trans, &mi, &ni, &i__1, &a[a_dim1 + 2], lda, &tau[1]
                    , &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &iinfo, (
                    ftnlen)1, (ftnlen)1);
        }
    } else {

/*        Apply P */

        if (notran) {
            *(unsigned char *)transt = 'T';
        } else {
            *(unsigned char *)transt = 'N';
        }
        if (nq > *k) {

/*           P was determined by a call to DGEBRD with nq > k */

            dormlq_(side, transt, m, n, k, &a[a_offset], lda, &tau[1], &c__[
                    c_offset], ldc, &work[1], lwork, &iinfo, (ftnlen)1, (
                    ftnlen)1);
        } else if (nq > 1) {

/*           P was determined by a call to DGEBRD with nq <= k */

            if (left) {
                mi = *m - 1;
                ni = *n;
                i1 = 2;
                i2 = 1;
            } else {
                mi = *m;
                ni = *n - 1;
                i1 = 1;
                i2 = 2;
            }
            i__1 = nq - 1;
            dormlq_(side, transt, &mi, &ni, &i__1, &a[(a_dim1 << 1) + 1], lda,
                     &tau[1], &c__[i1 + i2 * c_dim1], ldc, &work[1], lwork, &
                    iinfo, (ftnlen)1, (ftnlen)1);
        }
    }
    work[1] = (doublereal) lwkopt;
    return 0;

/*     End of DORMBR */

} /* dormbr_ */

#ifdef __cplusplus
        }
#endif
