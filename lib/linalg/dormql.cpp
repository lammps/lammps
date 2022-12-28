/* fortran/dormql.f -- translated by f2c (version 20200916).
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
static integer c__65 = 65;

/* > \brief \b DORMQL */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DORMQL + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormql.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormql.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormql.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DORMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, */
/*                          WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE, TRANS */
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
/* > DORMQL overwrites the general real M-by-N matrix C with */
/* > */
/* >                 SIDE = 'L'     SIDE = 'R' */
/* > TRANS = 'N':      Q * C          C * Q */
/* > TRANS = 'T':      Q**T * C       C * Q**T */
/* > */
/* > where Q is a real orthogonal matrix defined as the product of k */
/* > elementary reflectors */
/* > */
/* >       Q = H(k) . . . H(2) H(1) */
/* > */
/* > as returned by DGEQLF. Q is of order M if SIDE = 'L' and of order N */
/* > if SIDE = 'R'. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply Q or Q**T from the Left; */
/* >          = 'R': apply Q or Q**T from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N':  No transpose, apply Q; */
/* >          = 'T':  Transpose, apply Q**T. */
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
/* >          The number of elementary reflectors whose product defines */
/* >          the matrix Q. */
/* >          If SIDE = 'L', M >= K >= 0; */
/* >          if SIDE = 'R', N >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,K) */
/* >          The i-th column must contain the vector which defines the */
/* >          elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* >          DGEQLF in the last k columns of its array argument A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. */
/* >          If SIDE = 'L', LDA >= max(1,M); */
/* >          if SIDE = 'R', LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is DOUBLE PRECISION array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by DGEQLF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
/* >          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. */
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
/* >          For good performance, LWORK should generally be larger. */
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
/* Subroutine */ int dormql_(char *side, char *trans, integer *m, integer *n,
        integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
        c__, integer *ldc, doublereal *work, integer *lwork, integer *info,
        ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3[2], i__4,
            i__5;
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    integer i__, i1, i2, i3, ib, nb, mi, ni, nq, nw, iwt;
    logical left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer nbmin, iinfo;
    extern /* Subroutine */ int dorm2l_(char *, char *, integer *, integer *,
            integer *, doublereal *, integer *, doublereal *, doublereal *,
            integer *, doublereal *, integer *, ftnlen, ftnlen), dlarfb_(char
            *, char *, char *, char *, integer *, integer *, integer *,
            doublereal *, integer *, doublereal *, integer *, doublereal *,
            integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen,
            ftnlen), dlarft_(char *, char *, integer *, integer *, doublereal
            *, integer *, doublereal *, doublereal *, integer *, ftnlen,
            ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *,
            integer *, integer *, ftnlen, ftnlen);
    logical notran;
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
    left = lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1);
    notran = lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1);
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

    if (left) {
        nq = *m;
        nw = max(1,*n);
    } else {
        nq = *n;
        nw = max(1,*m);
    }
    if (! left && ! lsame_(side, (char *)"R", (ftnlen)1, (ftnlen)1)) {
        *info = -1;
    } else if (! notran && ! lsame_(trans, (char *)"T", (ftnlen)1, (ftnlen)1)) {
        *info = -2;
    } else if (*m < 0) {
        *info = -3;
    } else if (*n < 0) {
        *info = -4;
    } else if (*k < 0 || *k > nq) {
        *info = -5;
    } else if (*lda < max(1,nq)) {
        *info = -7;
    } else if (*ldc < max(1,*m)) {
        *info = -10;
    } else if (*lwork < nw && ! lquery) {
        *info = -12;
    }

    if (*info == 0) {

/*        Compute the workspace requirements */

        if (*m == 0 || *n == 0) {
            lwkopt = 1;
        } else {
/* Computing MIN */
/* Writing concatenation */
            i__3[0] = 1, a__1[0] = side;
            i__3[1] = 1, a__1[1] = trans;
            s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
            i__1 = 64, i__2 = ilaenv_(&c__1, (char *)"DORMQL", ch__1, m, n, k, &c_n1,
                    (ftnlen)6, (ftnlen)2);
            nb = min(i__1,i__2);
            lwkopt = nw * nb + 4160;
        }
        work[1] = (doublereal) lwkopt;
    }

    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DORMQL", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
        return 0;
    }

    nbmin = 2;
    ldwork = nw;
    if (nb > 1 && nb < *k) {
        if (*lwork < lwkopt) {
            nb = (*lwork - 4160) / ldwork;
/* Computing MAX */
/* Writing concatenation */
            i__3[0] = 1, a__1[0] = side;
            i__3[1] = 1, a__1[1] = trans;
            s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
            i__1 = 2, i__2 = ilaenv_(&c__2, (char *)"DORMQL", ch__1, m, n, k, &c_n1, (
                    ftnlen)6, (ftnlen)2);
            nbmin = max(i__1,i__2);
        }
    }

    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

        dorm2l_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
                c_offset], ldc, &work[1], &iinfo, (ftnlen)1, (ftnlen)1);
    } else {

/*        Use blocked code */

        iwt = nw * nb + 1;
        if (left && notran || ! left && ! notran) {
            i1 = 1;
            i2 = *k;
            i3 = nb;
        } else {
            i1 = (*k - 1) / nb * nb + 1;
            i2 = 1;
            i3 = -nb;
        }

        if (left) {
            ni = *n;
        } else {
            mi = *m;
        }

        i__1 = i2;
        i__2 = i3;
        for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
            i__4 = nb, i__5 = *k - i__ + 1;
            ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector */
/*           H = H(i+ib-1) . . . H(i+1) H(i) */

            i__4 = nq - *k + i__ + ib - 1;
            dlarft_((char *)"Backward", (char *)"Columnwise", &i__4, &ib, &a[i__ * a_dim1 + 1]
                    , lda, &tau[i__], &work[iwt], &c__65, (ftnlen)8, (ftnlen)
                    10);
            if (left) {

/*              H or H**T is applied to C(1:m-k+i+ib-1,1:n) */

                mi = *m - *k + i__ + ib - 1;
            } else {

/*              H or H**T is applied to C(1:m,1:n-k+i+ib-1) */

                ni = *n - *k + i__ + ib - 1;
            }

/*           Apply H or H**T */

            dlarfb_(side, trans, (char *)"Backward", (char *)"Columnwise", &mi, &ni, &ib, &a[
                    i__ * a_dim1 + 1], lda, &work[iwt], &c__65, &c__[c_offset]
                    , ldc, &work[1], &ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)8,
                     (ftnlen)10);
/* L10: */
        }
    }
    work[1] = (doublereal) lwkopt;
    return 0;

/*     End of DORMQL */

} /* dormql_ */

#ifdef __cplusplus
        }
#endif
