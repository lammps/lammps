/* fortran/dgesvd.f -- translated by f2c (version 20200916).
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

static integer c__6 = 6;
static integer c__0 = 0;
static integer c__2 = 2;
static integer c_n1 = -1;
static doublereal c_b57 = 0.;
static integer c__1 = 1;
static doublereal c_b79 = 1.;

/* > \brief <b> DGESVD computes the singular value decomposition (SVD) for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGESVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesvd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesvd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesvd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, */
/*                          WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBU, JOBVT */
/*       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ), */
/*      $                   VT( LDVT, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGESVD computes the singular value decomposition (SVD) of a real */
/* > M-by-N matrix A, optionally computing the left and/or right singular */
/* > vectors. The SVD is written */
/* > */
/* >      A = U * SIGMA * transpose(V) */
/* > */
/* > where SIGMA is an M-by-N matrix which is zero except for its */
/* > min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and */
/* > V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA */
/* > are the singular values of A; they are real and non-negative, and */
/* > are returned in descending order.  The first min(m,n) columns of */
/* > U and V are the left and right singular vectors of A. */
/* > */
/* > Note that the routine returns V**T, not V. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBU */
/* > \verbatim */
/* >          JOBU is CHARACTER*1 */
/* >          Specifies options for computing all or part of the matrix U: */
/* >          = 'A':  all M columns of U are returned in array U: */
/* >          = 'S':  the first min(m,n) columns of U (the left singular */
/* >                  vectors) are returned in the array U; */
/* >          = 'O':  the first min(m,n) columns of U (the left singular */
/* >                  vectors) are overwritten on the array A; */
/* >          = 'N':  no columns of U (no left singular vectors) are */
/* >                  computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBVT */
/* > \verbatim */
/* >          JOBVT is CHARACTER*1 */
/* >          Specifies options for computing all or part of the matrix */
/* >          V**T: */
/* >          = 'A':  all N rows of V**T are returned in the array VT; */
/* >          = 'S':  the first min(m,n) rows of V**T (the right singular */
/* >                  vectors) are returned in the array VT; */
/* >          = 'O':  the first min(m,n) rows of V**T (the right singular */
/* >                  vectors) are overwritten on the array A; */
/* >          = 'N':  no rows of V**T (no right singular vectors) are */
/* >                  computed. */
/* > */
/* >          JOBVT and JOBU cannot both be 'O'. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the input matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the input matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, */
/* >          if JOBU = 'O',  A is overwritten with the first min(m,n) */
/* >                          columns of U (the left singular vectors, */
/* >                          stored columnwise); */
/* >          if JOBVT = 'O', A is overwritten with the first min(m,n) */
/* >                          rows of V**T (the right singular vectors, */
/* >                          stored rowwise); */
/* >          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A */
/* >                          are destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is DOUBLE PRECISION array, dimension (min(M,N)) */
/* >          The singular values of A, sorted so that S(i) >= S(i+1). */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is DOUBLE PRECISION array, dimension (LDU,UCOL) */
/* >          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'. */
/* >          If JOBU = 'A', U contains the M-by-M orthogonal matrix U; */
/* >          if JOBU = 'S', U contains the first min(m,n) columns of U */
/* >          (the left singular vectors, stored columnwise); */
/* >          if JOBU = 'N' or 'O', U is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* >          LDU is INTEGER */
/* >          The leading dimension of the array U.  LDU >= 1; if */
/* >          JOBU = 'S' or 'A', LDU >= M. */
/* > \endverbatim */
/* > */
/* > \param[out] VT */
/* > \verbatim */
/* >          VT is DOUBLE PRECISION array, dimension (LDVT,N) */
/* >          If JOBVT = 'A', VT contains the N-by-N orthogonal matrix */
/* >          V**T; */
/* >          if JOBVT = 'S', VT contains the first min(m,n) rows of */
/* >          V**T (the right singular vectors, stored rowwise); */
/* >          if JOBVT = 'N' or 'O', VT is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVT */
/* > \verbatim */
/* >          LDVT is INTEGER */
/* >          The leading dimension of the array VT.  LDVT >= 1; if */
/* >          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK; */
/* >          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged */
/* >          superdiagonal elements of an upper bidiagonal matrix B */
/* >          whose diagonal is in S (not necessarily sorted). B */
/* >          satisfies A = U * B * VT, so it has the same singular values */
/* >          as A, and singular vectors related by U and VT. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          LWORK >= MAX(1,5*MIN(M,N)) for the paths (see comments inside code): */
/* >             - PATH 1  (M much larger than N, JOBU='N') */
/* >             - PATH 1t (N much larger than M, JOBVT='N') */
/* >          LWORK >= MAX(1,3*MIN(M,N) + MAX(M,N),5*MIN(M,N)) for the other paths */
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
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  if DBDSQR did not converge, INFO specifies how many */
/* >                superdiagonals of an intermediate bidiagonal form B */
/* >                did not converge to zero. See the description of WORK */
/* >                above for details. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doubleGEsing */

/*  ===================================================================== */
/* Subroutine */ int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n,
        doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
        ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork,
        integer *info, ftnlen jobu_len, ftnlen jobvt_len)
{
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, u_dim1, u_offset, vt_dim1, vt_offset, i__1[2],
            i__2, i__3, i__4;
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    integer i__, ie, ir, iu, blk, ncu;
    doublereal dum[1], eps;
    integer nru, iscl;
    doublereal anrm;
    integer ierr, itau, ncvt, nrvt, lwork_dgebrd__, lwork_dgelqf__,
            lwork_dgeqrf__;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *,
            integer *, doublereal *, doublereal *, integer *, doublereal *,
            integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer chunk, minmn, wrkbl, itaup, itauq, mnthr, iwork;
    logical wntua, wntva, wntun, wntuo, wntvn, wntvo, wntus, wntvs;
    extern /* Subroutine */ int dgebrd_(integer *, integer *, doublereal *,
            integer *, doublereal *, doublereal *, doublereal *, doublereal *,
             doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *,
            integer *, doublereal *, integer *, doublereal *, ftnlen);
    integer bdspac;
    extern /* Subroutine */ int dgelqf_(integer *, integer *, doublereal *,
            integer *, doublereal *, doublereal *, integer *, integer *),
            dlascl_(char *, integer *, integer *, doublereal *, doublereal *,
            integer *, integer *, doublereal *, integer *, integer *, ftnlen),
             dgeqrf_(integer *, integer *, doublereal *, integer *,
            doublereal *, doublereal *, integer *, integer *), dlacpy_(char *,
             integer *, integer *, doublereal *, integer *, doublereal *,
            integer *, ftnlen), dlaset_(char *, integer *, integer *,
            doublereal *, doublereal *, doublereal *, integer *, ftnlen),
            dbdsqr_(char *, integer *, integer *, integer *, integer *,
            doublereal *, doublereal *, doublereal *, integer *, doublereal *,
             integer *, doublereal *, integer *, doublereal *, integer *,
            ftnlen), dorgbr_(char *, integer *, integer *, integer *,
            doublereal *, integer *, doublereal *, doublereal *, integer *,
            integer *, ftnlen);
    doublereal bignum;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *,
            integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dormbr_(char *, char *, char *, integer *,
            integer *, integer *, doublereal *, integer *, doublereal *,
            doublereal *, integer *, doublereal *, integer *, integer *,
            ftnlen, ftnlen, ftnlen), dorglq_(integer *, integer *, integer *,
            doublereal *, integer *, doublereal *, doublereal *, integer *,
            integer *), dorgqr_(integer *, integer *, integer *, doublereal *,
             integer *, doublereal *, doublereal *, integer *, integer *);
    integer ldwrkr, minwrk, ldwrku, maxwrk;
    doublereal smlnum;
    logical lquery, wntuas, wntvas;
    integer lwork_dorgbr_p__, lwork_dorgbr_q__, lwork_dorglq_m__,
            lwork_dorglq_n__, lwork_dorgqr_m__, lwork_dorgqr_n__;


/*  -- LAPACK driver routine -- */
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
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --s;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --work;

    /* Function Body */
    *info = 0;
    minmn = min(*m,*n);
    wntua = lsame_(jobu, (char *)"A", (ftnlen)1, (ftnlen)1);
    wntus = lsame_(jobu, (char *)"S", (ftnlen)1, (ftnlen)1);
    wntuas = wntua || wntus;
    wntuo = lsame_(jobu, (char *)"O", (ftnlen)1, (ftnlen)1);
    wntun = lsame_(jobu, (char *)"N", (ftnlen)1, (ftnlen)1);
    wntva = lsame_(jobvt, (char *)"A", (ftnlen)1, (ftnlen)1);
    wntvs = lsame_(jobvt, (char *)"S", (ftnlen)1, (ftnlen)1);
    wntvas = wntva || wntvs;
    wntvo = lsame_(jobvt, (char *)"O", (ftnlen)1, (ftnlen)1);
    wntvn = lsame_(jobvt, (char *)"N", (ftnlen)1, (ftnlen)1);
    lquery = *lwork == -1;

    if (! (wntua || wntus || wntuo || wntun)) {
        *info = -1;
    } else if (! (wntva || wntvs || wntvo || wntvn) || wntvo && wntuo) {
        *info = -2;
    } else if (*m < 0) {
        *info = -3;
    } else if (*n < 0) {
        *info = -4;
    } else if (*lda < max(1,*m)) {
        *info = -6;
    } else if (*ldu < 1 || wntuas && *ldu < *m) {
        *info = -9;
    } else if (*ldvt < 1 || wntva && *ldvt < *n || wntvs && *ldvt < minmn) {
        *info = -11;
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning (char *)"Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV.) */

    if (*info == 0) {
        minwrk = 1;
        maxwrk = 1;
        if (*m >= *n && minmn > 0) {

/*           Compute space needed for DBDSQR */

/* Writing concatenation */
            i__1[0] = 1, a__1[0] = jobu;
            i__1[1] = 1, a__1[1] = jobvt;
            s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
            mnthr = ilaenv_(&c__6, (char *)"DGESVD", ch__1, m, n, &c__0, &c__0, (
                    ftnlen)6, (ftnlen)2);
            bdspac = *n * 5;
/*           Compute space needed for DGEQRF */
            dgeqrf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
            lwork_dgeqrf__ = (integer) dum[0];
/*           Compute space needed for DORGQR */
            dorgqr_(m, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
            lwork_dorgqr_n__ = (integer) dum[0];
            dorgqr_(m, m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
            lwork_dorgqr_m__ = (integer) dum[0];
/*           Compute space needed for DGEBRD */
            dgebrd_(n, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &c_n1,
                     &ierr);
            lwork_dgebrd__ = (integer) dum[0];
/*           Compute space needed for DORGBR P */
            dorgbr_((char *)"P", n, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr, (
                    ftnlen)1);
            lwork_dorgbr_p__ = (integer) dum[0];
/*           Compute space needed for DORGBR Q */
            dorgbr_((char *)"Q", n, n, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr, (
                    ftnlen)1);
            lwork_dorgbr_q__ = (integer) dum[0];

            if (*m >= mnthr) {
                if (wntun) {

/*                 Path 1 (M much larger than N, JOBU='N') */

                    maxwrk = *n + lwork_dgeqrf__;
/* Computing MAX */
                    i__2 = maxwrk, i__3 = *n * 3 + lwork_dgebrd__;
                    maxwrk = max(i__2,i__3);
                    if (wntvo || wntvas) {
/* Computing MAX */
                        i__2 = maxwrk, i__3 = *n * 3 + lwork_dorgbr_p__;
                        maxwrk = max(i__2,i__3);
                    }
                    maxwrk = max(maxwrk,bdspac);
/* Computing MAX */
                    i__2 = *n << 2;
                    minwrk = max(i__2,bdspac);
                } else if (wntuo && wntvn) {

/*                 Path 2 (M much larger than N, JOBU='O', JOBVT='N') */

                    wrkbl = *n + lwork_dgeqrf__;
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_n__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2,i__3);
                    wrkbl = max(wrkbl,bdspac);
/* Computing MAX */
                    i__2 = *n * *n + wrkbl, i__3 = *n * *n + *m * *n + *n;
                    maxwrk = max(i__2,i__3);
/* Computing MAX */
                    i__2 = *n * 3 + *m;
                    minwrk = max(i__2,bdspac);
                } else if (wntuo && wntvas) {

/*                 Path 3 (M much larger than N, JOBU='O', JOBVT='S' or */
/*                 'A') */

                    wrkbl = *n + lwork_dgeqrf__;
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_n__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2,i__3);
                    wrkbl = max(wrkbl,bdspac);
/* Computing MAX */
                    i__2 = *n * *n + wrkbl, i__3 = *n * *n + *m * *n + *n;
                    maxwrk = max(i__2,i__3);
/* Computing MAX */
                    i__2 = *n * 3 + *m;
                    minwrk = max(i__2,bdspac);
                } else if (wntus && wntvn) {

/*                 Path 4 (M much larger than N, JOBU='S', JOBVT='N') */

                    wrkbl = *n + lwork_dgeqrf__;
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_n__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2,i__3);
                    wrkbl = max(wrkbl,bdspac);
                    maxwrk = *n * *n + wrkbl;
/* Computing MAX */
                    i__2 = *n * 3 + *m;
                    minwrk = max(i__2,bdspac);
                } else if (wntus && wntvo) {

/*                 Path 5 (M much larger than N, JOBU='S', JOBVT='O') */

                    wrkbl = *n + lwork_dgeqrf__;
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_n__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2,i__3);
                    wrkbl = max(wrkbl,bdspac);
                    maxwrk = (*n << 1) * *n + wrkbl;
/* Computing MAX */
                    i__2 = *n * 3 + *m;
                    minwrk = max(i__2,bdspac);
                } else if (wntus && wntvas) {

/*                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' or */
/*                 'A') */

                    wrkbl = *n + lwork_dgeqrf__;
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_n__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2,i__3);
                    wrkbl = max(wrkbl,bdspac);
                    maxwrk = *n * *n + wrkbl;
/* Computing MAX */
                    i__2 = *n * 3 + *m;
                    minwrk = max(i__2,bdspac);
                } else if (wntua && wntvn) {

/*                 Path 7 (M much larger than N, JOBU='A', JOBVT='N') */

                    wrkbl = *n + lwork_dgeqrf__;
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_m__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2,i__3);
                    wrkbl = max(wrkbl,bdspac);
                    maxwrk = *n * *n + wrkbl;
/* Computing MAX */
                    i__2 = *n * 3 + *m;
                    minwrk = max(i__2,bdspac);
                } else if (wntua && wntvo) {

/*                 Path 8 (M much larger than N, JOBU='A', JOBVT='O') */

                    wrkbl = *n + lwork_dgeqrf__;
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_m__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2,i__3);
                    wrkbl = max(wrkbl,bdspac);
                    maxwrk = (*n << 1) * *n + wrkbl;
/* Computing MAX */
                    i__2 = *n * 3 + *m;
                    minwrk = max(i__2,bdspac);
                } else if (wntua && wntvas) {

/*                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' or */
/*                 'A') */

                    wrkbl = *n + lwork_dgeqrf__;
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n + lwork_dorgqr_m__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *n * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2,i__3);
                    wrkbl = max(wrkbl,bdspac);
                    maxwrk = *n * *n + wrkbl;
/* Computing MAX */
                    i__2 = *n * 3 + *m;
                    minwrk = max(i__2,bdspac);
                }
            } else {

/*              Path 10 (M at least N, but not much larger) */

                dgebrd_(m, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &
                        c_n1, &ierr);
                lwork_dgebrd__ = (integer) dum[0];
                maxwrk = *n * 3 + lwork_dgebrd__;
                if (wntus || wntuo) {
                    dorgbr_((char *)"Q", m, n, n, &a[a_offset], lda, dum, dum, &c_n1,
                            &ierr, (ftnlen)1);
                    lwork_dorgbr_q__ = (integer) dum[0];
/* Computing MAX */
                    i__2 = maxwrk, i__3 = *n * 3 + lwork_dorgbr_q__;
                    maxwrk = max(i__2,i__3);
                }
                if (wntua) {
                    dorgbr_((char *)"Q", m, m, n, &a[a_offset], lda, dum, dum, &c_n1,
                            &ierr, (ftnlen)1);
                    lwork_dorgbr_q__ = (integer) dum[0];
/* Computing MAX */
                    i__2 = maxwrk, i__3 = *n * 3 + lwork_dorgbr_q__;
                    maxwrk = max(i__2,i__3);
                }
                if (! wntvn) {
/* Computing MAX */
                    i__2 = maxwrk, i__3 = *n * 3 + lwork_dorgbr_p__;
                    maxwrk = max(i__2,i__3);
                }
                maxwrk = max(maxwrk,bdspac);
/* Computing MAX */
                i__2 = *n * 3 + *m;
                minwrk = max(i__2,bdspac);
            }
        } else if (minmn > 0) {

/*           Compute space needed for DBDSQR */

/* Writing concatenation */
            i__1[0] = 1, a__1[0] = jobu;
            i__1[1] = 1, a__1[1] = jobvt;
            s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
            mnthr = ilaenv_(&c__6, (char *)"DGESVD", ch__1, m, n, &c__0, &c__0, (
                    ftnlen)6, (ftnlen)2);
            bdspac = *m * 5;
/*           Compute space needed for DGELQF */
            dgelqf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
            lwork_dgelqf__ = (integer) dum[0];
/*           Compute space needed for DORGLQ */
            dorglq_(n, n, m, dum, n, dum, dum, &c_n1, &ierr);
            lwork_dorglq_n__ = (integer) dum[0];
            dorglq_(m, n, m, &a[a_offset], lda, dum, dum, &c_n1, &ierr);
            lwork_dorglq_m__ = (integer) dum[0];
/*           Compute space needed for DGEBRD */
            dgebrd_(m, m, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &c_n1,
                     &ierr);
            lwork_dgebrd__ = (integer) dum[0];
/*            Compute space needed for DORGBR P */
            dorgbr_((char *)"P", m, m, m, &a[a_offset], n, dum, dum, &c_n1, &ierr, (
                    ftnlen)1);
            lwork_dorgbr_p__ = (integer) dum[0];
/*           Compute space needed for DORGBR Q */
            dorgbr_((char *)"Q", m, m, m, &a[a_offset], n, dum, dum, &c_n1, &ierr, (
                    ftnlen)1);
            lwork_dorgbr_q__ = (integer) dum[0];
            if (*n >= mnthr) {
                if (wntvn) {

/*                 Path 1t(N much larger than M, JOBVT='N') */

                    maxwrk = *m + lwork_dgelqf__;
/* Computing MAX */
                    i__2 = maxwrk, i__3 = *m * 3 + lwork_dgebrd__;
                    maxwrk = max(i__2,i__3);
                    if (wntuo || wntuas) {
/* Computing MAX */
                        i__2 = maxwrk, i__3 = *m * 3 + lwork_dorgbr_q__;
                        maxwrk = max(i__2,i__3);
                    }
                    maxwrk = max(maxwrk,bdspac);
/* Computing MAX */
                    i__2 = *m << 2;
                    minwrk = max(i__2,bdspac);
                } else if (wntvo && wntun) {

/*                 Path 2t(N much larger than M, JOBU='N', JOBVT='O') */

                    wrkbl = *m + lwork_dgelqf__;
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m + lwork_dorglq_m__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2,i__3);
                    wrkbl = max(wrkbl,bdspac);
/* Computing MAX */
                    i__2 = *m * *m + wrkbl, i__3 = *m * *m + *m * *n + *m;
                    maxwrk = max(i__2,i__3);
/* Computing MAX */
                    i__2 = *m * 3 + *n;
                    minwrk = max(i__2,bdspac);
                } else if (wntvo && wntuas) {

/*                 Path 3t(N much larger than M, JOBU='S' or 'A', */
/*                 JOBVT='O') */

                    wrkbl = *m + lwork_dgelqf__;
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m + lwork_dorglq_m__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2,i__3);
                    wrkbl = max(wrkbl,bdspac);
/* Computing MAX */
                    i__2 = *m * *m + wrkbl, i__3 = *m * *m + *m * *n + *m;
                    maxwrk = max(i__2,i__3);
/* Computing MAX */
                    i__2 = *m * 3 + *n;
                    minwrk = max(i__2,bdspac);
                } else if (wntvs && wntun) {

/*                 Path 4t(N much larger than M, JOBU='N', JOBVT='S') */

                    wrkbl = *m + lwork_dgelqf__;
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m + lwork_dorglq_m__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2,i__3);
                    wrkbl = max(wrkbl,bdspac);
                    maxwrk = *m * *m + wrkbl;
/* Computing MAX */
                    i__2 = *m * 3 + *n;
                    minwrk = max(i__2,bdspac);
                } else if (wntvs && wntuo) {

/*                 Path 5t(N much larger than M, JOBU='O', JOBVT='S') */

                    wrkbl = *m + lwork_dgelqf__;
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m + lwork_dorglq_m__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2,i__3);
                    wrkbl = max(wrkbl,bdspac);
                    maxwrk = (*m << 1) * *m + wrkbl;
/* Computing MAX */
                    i__2 = *m * 3 + *n;
                    minwrk = max(i__2,bdspac);
                } else if (wntvs && wntuas) {

/*                 Path 6t(N much larger than M, JOBU='S' or 'A', */
/*                 JOBVT='S') */

                    wrkbl = *m + lwork_dgelqf__;
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m + lwork_dorglq_m__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2,i__3);
                    wrkbl = max(wrkbl,bdspac);
                    maxwrk = *m * *m + wrkbl;
/* Computing MAX */
                    i__2 = *m * 3 + *n;
                    minwrk = max(i__2,bdspac);
                } else if (wntva && wntun) {

/*                 Path 7t(N much larger than M, JOBU='N', JOBVT='A') */

                    wrkbl = *m + lwork_dgelqf__;
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m + lwork_dorglq_n__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2,i__3);
                    wrkbl = max(wrkbl,bdspac);
                    maxwrk = *m * *m + wrkbl;
/* Computing MAX */
                    i__2 = *m * 3 + *n;
                    minwrk = max(i__2,bdspac);
                } else if (wntva && wntuo) {

/*                 Path 8t(N much larger than M, JOBU='O', JOBVT='A') */

                    wrkbl = *m + lwork_dgelqf__;
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m + lwork_dorglq_n__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2,i__3);
                    wrkbl = max(wrkbl,bdspac);
                    maxwrk = (*m << 1) * *m + wrkbl;
/* Computing MAX */
                    i__2 = *m * 3 + *n;
                    minwrk = max(i__2,bdspac);
                } else if (wntva && wntuas) {

/*                 Path 9t(N much larger than M, JOBU='S' or 'A', */
/*                 JOBVT='A') */

                    wrkbl = *m + lwork_dgelqf__;
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m + lwork_dorglq_n__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dgebrd__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_p__;
                    wrkbl = max(i__2,i__3);
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *m * 3 + lwork_dorgbr_q__;
                    wrkbl = max(i__2,i__3);
                    wrkbl = max(wrkbl,bdspac);
                    maxwrk = *m * *m + wrkbl;
/* Computing MAX */
                    i__2 = *m * 3 + *n;
                    minwrk = max(i__2,bdspac);
                }
            } else {

/*              Path 10t(N greater than M, but not much larger) */

                dgebrd_(m, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &
                        c_n1, &ierr);
                lwork_dgebrd__ = (integer) dum[0];
                maxwrk = *m * 3 + lwork_dgebrd__;
                if (wntvs || wntvo) {
/*                Compute space needed for DORGBR P */
                    dorgbr_((char *)"P", m, n, m, &a[a_offset], n, dum, dum, &c_n1, &
                            ierr, (ftnlen)1);
                    lwork_dorgbr_p__ = (integer) dum[0];
/* Computing MAX */
                    i__2 = maxwrk, i__3 = *m * 3 + lwork_dorgbr_p__;
                    maxwrk = max(i__2,i__3);
                }
                if (wntva) {
                    dorgbr_((char *)"P", n, n, m, &a[a_offset], n, dum, dum, &c_n1, &
                            ierr, (ftnlen)1);
                    lwork_dorgbr_p__ = (integer) dum[0];
/* Computing MAX */
                    i__2 = maxwrk, i__3 = *m * 3 + lwork_dorgbr_p__;
                    maxwrk = max(i__2,i__3);
                }
                if (! wntun) {
/* Computing MAX */
                    i__2 = maxwrk, i__3 = *m * 3 + lwork_dorgbr_q__;
                    maxwrk = max(i__2,i__3);
                }
                maxwrk = max(maxwrk,bdspac);
/* Computing MAX */
                i__2 = *m * 3 + *n;
                minwrk = max(i__2,bdspac);
            }
        }
        maxwrk = max(maxwrk,minwrk);
        work[1] = (doublereal) maxwrk;

        if (*lwork < minwrk && ! lquery) {
            *info = -13;
        }
    }

    if (*info != 0) {
        i__2 = -(*info);
        xerbla_((char *)"DGESVD", &i__2, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
        return 0;
    }

/*     Get machine constants */

    eps = dlamch_((char *)"P", (ftnlen)1);
    smlnum = sqrt(dlamch_((char *)"S", (ftnlen)1)) / eps;
    bignum = 1. / smlnum;

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

    anrm = dlange_((char *)"M", m, n, &a[a_offset], lda, dum, (ftnlen)1);
    iscl = 0;
    if (anrm > 0. && anrm < smlnum) {
        iscl = 1;
        dlascl_((char *)"G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, &
                ierr, (ftnlen)1);
    } else if (anrm > bignum) {
        iscl = 1;
        dlascl_((char *)"G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, &
                ierr, (ftnlen)1);
    }

    if (*m >= *n) {

/*        A has at least as many rows as columns. If A has sufficiently */
/*        more rows than columns, first reduce using the QR */
/*        decomposition (if sufficient workspace available) */

        if (*m >= mnthr) {

            if (wntun) {

/*              Path 1 (M much larger than N, JOBU='N') */
/*              No left singular vectors to be computed */

                itau = 1;
                iwork = itau + *n;

/*              Compute A=Q*R */
/*              (Workspace: need 2*N, prefer N + N*NB) */

                i__2 = *lwork - iwork + 1;
                dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &
                        i__2, &ierr);

/*              Zero out below R */

                if (*n > 1) {
                    i__2 = *n - 1;
                    i__3 = *n - 1;
                    dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &a[a_dim1 + 2],
                             lda, (ftnlen)1);
                }
                ie = 1;
                itauq = ie + *n;
                itaup = itauq + *n;
                iwork = itaup + *n;

/*              Bidiagonalize R in A */
/*              (Workspace: need 4*N, prefer 3*N + 2*N*NB) */

                i__2 = *lwork - iwork + 1;
                dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &work[
                        itauq], &work[itaup], &work[iwork], &i__2, &ierr);
                ncvt = 0;
                if (wntvo || wntvas) {

/*                 If right singular vectors desired, generate P'. */
/*                 (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB) */

                    i__2 = *lwork - iwork + 1;
                    dorgbr_((char *)"P", n, n, n, &a[a_offset], lda, &work[itaup], &
                            work[iwork], &i__2, &ierr, (ftnlen)1);
                    ncvt = *n;
                }
                iwork = ie + *n;

/*              Perform bidiagonal QR iteration, computing right */
/*              singular vectors of A in A if desired */
/*              (Workspace: need BDSPAC) */

                dbdsqr_((char *)"U", n, &ncvt, &c__0, &c__0, &s[1], &work[ie], &a[
                        a_offset], lda, dum, &c__1, dum, &c__1, &work[iwork],
                        info, (ftnlen)1);

/*              If right singular vectors desired in VT, copy them there */

                if (wntvas) {
                    dlacpy_((char *)"F", n, n, &a[a_offset], lda, &vt[vt_offset],
                            ldvt, (ftnlen)1);
                }

            } else if (wntuo && wntvn) {

/*              Path 2 (M much larger than N, JOBU='O', JOBVT='N') */
/*              N left singular vectors to be overwritten on A and */
/*              no right singular vectors to be computed */

/* Computing MAX */
                i__2 = *n << 2;
                if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                 Sufficient workspace for a fast algorithm */

                    ir = 1;
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *lda * *n + *n;
                    if (*lwork >= max(i__2,i__3) + *lda * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) is LDA by N */

                        ldwrku = *lda;
                        ldwrkr = *lda;
                    } else /* if(complicated condition) */ {
/* Computing MAX */
                        i__2 = wrkbl, i__3 = *lda * *n + *n;
                        if (*lwork >= max(i__2,i__3) + *n * *n) {

/*                    WORK(IU) is LDA by N, WORK(IR) is N by N */

                            ldwrku = *lda;
                            ldwrkr = *n;
                        } else {

/*                    WORK(IU) is LDWRKU by N, WORK(IR) is N by N */

                            ldwrku = (*lwork - *n * *n - *n) / *n;
                            ldwrkr = *n;
                        }
                    }
                    itau = ir + ldwrkr * *n;
                    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (Workspace: need N*N + 2*N, prefer N*N + N + N*NB) */

                    i__2 = *lwork - iwork + 1;
                    dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
                            , &i__2, &ierr);

/*                 Copy R to WORK(IR) and zero out below it */

                    dlacpy_((char *)"U", n, n, &a[a_offset], lda, &work[ir], &ldwrkr,
                            (ftnlen)1);
                    i__2 = *n - 1;
                    i__3 = *n - 1;
                    dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &work[ir + 1],
                            &ldwrkr, (ftnlen)1);

/*                 Generate Q in A */
/*                 (Workspace: need N*N + 2*N, prefer N*N + N + N*NB) */

                    i__2 = *lwork - iwork + 1;
                    dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
                            iwork], &i__2, &ierr);
                    ie = itau;
                    itauq = ie + *n;
                    itaup = itauq + *n;
                    iwork = itaup + *n;

/*                 Bidiagonalize R in WORK(IR) */
/*                 (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB) */

                    i__2 = *lwork - iwork + 1;
                    dgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &work[
                            itauq], &work[itaup], &work[iwork], &i__2, &ierr);

/*                 Generate left vectors bidiagonalizing R */
/*                 (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB) */

                    i__2 = *lwork - iwork + 1;
                    dorgbr_((char *)"Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &
                            work[iwork], &i__2, &ierr, (ftnlen)1);
                    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of R in WORK(IR) */
/*                 (Workspace: need N*N + BDSPAC) */

                    dbdsqr_((char *)"U", n, &c__0, n, &c__0, &s[1], &work[ie], dum, &
                            c__1, &work[ir], &ldwrkr, dum, &c__1, &work[iwork]
                            , info, (ftnlen)1);
                    iu = ie + *n;

/*                 Multiply Q in A by left singular vectors of R in */
/*                 WORK(IR), storing result in WORK(IU) and copying to A */
/*                 (Workspace: need N*N + 2*N, prefer N*N + M*N + N) */

                    i__2 = *m;
                    i__3 = ldwrku;
                    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
                             i__3) {
/* Computing MIN */
                        i__4 = *m - i__ + 1;
                        chunk = min(i__4,ldwrku);
                        dgemm_((char *)"N", (char *)"N", &chunk, n, n, &c_b79, &a[i__ +
                                a_dim1], lda, &work[ir], &ldwrkr, &c_b57, &
                                work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
                        dlacpy_((char *)"F", &chunk, n, &work[iu], &ldwrku, &a[i__ +
                                a_dim1], lda, (ftnlen)1);
/* L10: */
                    }

                } else {

/*                 Insufficient workspace for a fast algorithm */

                    ie = 1;
                    itauq = ie + *n;
                    itaup = itauq + *n;
                    iwork = itaup + *n;

/*                 Bidiagonalize A */
/*                 (Workspace: need 3*N + M, prefer 3*N + (M + N)*NB) */

                    i__3 = *lwork - iwork + 1;
                    dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[
                            itauq], &work[itaup], &work[iwork], &i__3, &ierr);

/*                 Generate left vectors bidiagonalizing A */
/*                 (Workspace: need 4*N, prefer 3*N + N*NB) */

                    i__3 = *lwork - iwork + 1;
                    dorgbr_((char *)"Q", m, n, n, &a[a_offset], lda, &work[itauq], &
                            work[iwork], &i__3, &ierr, (ftnlen)1);
                    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in A */
/*                 (Workspace: need BDSPAC) */

                    dbdsqr_((char *)"U", n, &c__0, m, &c__0, &s[1], &work[ie], dum, &
                            c__1, &a[a_offset], lda, dum, &c__1, &work[iwork],
                             info, (ftnlen)1);

                }

            } else if (wntuo && wntvas) {

/*              Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A') */
/*              N left singular vectors to be overwritten on A and */
/*              N right singular vectors to be computed in VT */

/* Computing MAX */
                i__3 = *n << 2;
                if (*lwork >= *n * *n + max(i__3,bdspac)) {

/*                 Sufficient workspace for a fast algorithm */

                    ir = 1;
/* Computing MAX */
                    i__3 = wrkbl, i__2 = *lda * *n + *n;
                    if (*lwork >= max(i__3,i__2) + *lda * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by N */

                        ldwrku = *lda;
                        ldwrkr = *lda;
                    } else /* if(complicated condition) */ {
/* Computing MAX */
                        i__3 = wrkbl, i__2 = *lda * *n + *n;
                        if (*lwork >= max(i__3,i__2) + *n * *n) {

/*                    WORK(IU) is LDA by N and WORK(IR) is N by N */

                            ldwrku = *lda;
                            ldwrkr = *n;
                        } else {

/*                    WORK(IU) is LDWRKU by N and WORK(IR) is N by N */

                            ldwrku = (*lwork - *n * *n - *n) / *n;
                            ldwrkr = *n;
                        }
                    }
                    itau = ir + ldwrkr * *n;
                    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (Workspace: need N*N + 2*N, prefer N*N + N + N*NB) */

                    i__3 = *lwork - iwork + 1;
                    dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
                            , &i__3, &ierr);

/*                 Copy R to VT, zeroing out below it */

                    dlacpy_((char *)"U", n, n, &a[a_offset], lda, &vt[vt_offset],
                            ldvt, (ftnlen)1);
                    if (*n > 1) {
                        i__3 = *n - 1;
                        i__2 = *n - 1;
                        dlaset_((char *)"L", &i__3, &i__2, &c_b57, &c_b57, &vt[
                                vt_dim1 + 2], ldvt, (ftnlen)1);
                    }

/*                 Generate Q in A */
/*                 (Workspace: need N*N + 2*N, prefer N*N + N + N*NB) */

                    i__3 = *lwork - iwork + 1;
                    dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
                            iwork], &i__3, &ierr);
                    ie = itau;
                    itauq = ie + *n;
                    itaup = itauq + *n;
                    iwork = itaup + *n;

/*                 Bidiagonalize R in VT, copying result to WORK(IR) */
/*                 (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB) */

                    i__3 = *lwork - iwork + 1;
                    dgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie], &
                            work[itauq], &work[itaup], &work[iwork], &i__3, &
                            ierr);
                    dlacpy_((char *)"L", n, n, &vt[vt_offset], ldvt, &work[ir], &
                            ldwrkr, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing R in WORK(IR) */
/*                 (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB) */

                    i__3 = *lwork - iwork + 1;
                    dorgbr_((char *)"Q", n, n, n, &work[ir], &ldwrkr, &work[itauq], &
                            work[iwork], &i__3, &ierr, (ftnlen)1);

/*                 Generate right vectors bidiagonalizing R in VT */
/*                 (Workspace: need N*N + 4*N-1, prefer N*N + 3*N + (N-1)*NB) */

                    i__3 = *lwork - iwork + 1;
                    dorgbr_((char *)"P", n, n, n, &vt[vt_offset], ldvt, &work[itaup],
                            &work[iwork], &i__3, &ierr, (ftnlen)1);
                    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of R in WORK(IR) and computing right */
/*                 singular vectors of R in VT */
/*                 (Workspace: need N*N + BDSPAC) */

                    dbdsqr_((char *)"U", n, n, n, &c__0, &s[1], &work[ie], &vt[
                            vt_offset], ldvt, &work[ir], &ldwrkr, dum, &c__1,
                            &work[iwork], info, (ftnlen)1);
                    iu = ie + *n;

/*                 Multiply Q in A by left singular vectors of R in */
/*                 WORK(IR), storing result in WORK(IU) and copying to A */
/*                 (Workspace: need N*N + 2*N, prefer N*N + M*N + N) */

                    i__3 = *m;
                    i__2 = ldwrku;
                    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ +=
                             i__2) {
/* Computing MIN */
                        i__4 = *m - i__ + 1;
                        chunk = min(i__4,ldwrku);
                        dgemm_((char *)"N", (char *)"N", &chunk, n, n, &c_b79, &a[i__ +
                                a_dim1], lda, &work[ir], &ldwrkr, &c_b57, &
                                work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
                        dlacpy_((char *)"F", &chunk, n, &work[iu], &ldwrku, &a[i__ +
                                a_dim1], lda, (ftnlen)1);
/* L20: */
                    }

                } else {

/*                 Insufficient workspace for a fast algorithm */

                    itau = 1;
                    iwork = itau + *n;

/*                 Compute A=Q*R */
/*                 (Workspace: need 2*N, prefer N + N*NB) */

                    i__2 = *lwork - iwork + 1;
                    dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
                            , &i__2, &ierr);

/*                 Copy R to VT, zeroing out below it */

                    dlacpy_((char *)"U", n, n, &a[a_offset], lda, &vt[vt_offset],
                            ldvt, (ftnlen)1);
                    if (*n > 1) {
                        i__2 = *n - 1;
                        i__3 = *n - 1;
                        dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &vt[
                                vt_dim1 + 2], ldvt, (ftnlen)1);
                    }

/*                 Generate Q in A */
/*                 (Workspace: need 2*N, prefer N + N*NB) */

                    i__2 = *lwork - iwork + 1;
                    dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &work[
                            iwork], &i__2, &ierr);
                    ie = itau;
                    itauq = ie + *n;
                    itaup = itauq + *n;
                    iwork = itaup + *n;

/*                 Bidiagonalize R in VT */
/*                 (Workspace: need 4*N, prefer 3*N + 2*N*NB) */

                    i__2 = *lwork - iwork + 1;
                    dgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie], &
                            work[itauq], &work[itaup], &work[iwork], &i__2, &
                            ierr);

/*                 Multiply Q in A by left vectors bidiagonalizing R */
/*                 (Workspace: need 3*N + M, prefer 3*N + M*NB) */

                    i__2 = *lwork - iwork + 1;
                    dormbr_((char *)"Q", (char *)"R", (char *)"N", m, n, n, &vt[vt_offset], ldvt, &
                            work[itauq], &a[a_offset], lda, &work[iwork], &
                            i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Generate right vectors bidiagonalizing R in VT */
/*                 (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB) */

                    i__2 = *lwork - iwork + 1;
                    dorgbr_((char *)"P", n, n, n, &vt[vt_offset], ldvt, &work[itaup],
                            &work[iwork], &i__2, &ierr, (ftnlen)1);
                    iwork = ie + *n;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in A and computing right */
/*                 singular vectors of A in VT */
/*                 (Workspace: need BDSPAC) */

                    dbdsqr_((char *)"U", n, n, m, &c__0, &s[1], &work[ie], &vt[
                            vt_offset], ldvt, &a[a_offset], lda, dum, &c__1, &
                            work[iwork], info, (ftnlen)1);

                }

            } else if (wntus) {

                if (wntvn) {

/*                 Path 4 (M much larger than N, JOBU='S', JOBVT='N') */
/*                 N left singular vectors to be computed in U and */
/*                 no right singular vectors to be computed */

/* Computing MAX */
                    i__2 = *n << 2;
                    if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

                        ir = 1;
                        if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

                            ldwrkr = *lda;
                        } else {

/*                       WORK(IR) is N by N */

                            ldwrkr = *n;
                        }
                        itau = ir + ldwrkr * *n;
                        iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);

/*                    Copy R to WORK(IR), zeroing out below it */

                        dlacpy_((char *)"U", n, n, &a[a_offset], lda, &work[ir], &
                                ldwrkr, (ftnlen)1);
                        i__2 = *n - 1;
                        i__3 = *n - 1;
                        dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &work[ir +
                                1], &ldwrkr, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &
                                work[iwork], &i__2, &ierr);
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR) */
/*                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &
                                work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);

/*                    Generate left vectors bidiagonalizing R in WORK(IR) */
/*                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", n, n, n, &work[ir], &ldwrkr, &work[itauq]
                                , &work[iwork], &i__2, &ierr, (ftnlen)1);
                        iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IR) */
/*                    (Workspace: need N*N + BDSPAC) */

                        dbdsqr_((char *)"U", n, &c__0, n, &c__0, &s[1], &work[ie],
                                dum, &c__1, &work[ir], &ldwrkr, dum, &c__1, &
                                work[iwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IR), storing result in U */
/*                    (Workspace: need N*N) */

                        dgemm_((char *)"N", (char *)"N", m, n, n, &c_b79, &a[a_offset], lda, &
                                work[ir], &ldwrkr, &c_b57, &u[u_offset], ldu,
                                (ftnlen)1, (ftnlen)1);

                    } else {

/*                    Insufficient workspace for a fast algorithm */

                        itau = 1;
                        iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, n, &a[a_offset], lda, &u[u_offset],
                                ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need 2*N, prefer N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
                                work[iwork], &i__2, &ierr);
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;

/*                    Zero out below R in A */

                        if (*n > 1) {
                            i__2 = *n - 1;
                            i__3 = *n - 1;
                            dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &a[
                                    a_dim1 + 2], lda, (ftnlen)1);
                        }

/*                    Bidiagonalize R in A */
/*                    (Workspace: need 4*N, prefer 3*N + 2*N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &
                                work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);

/*                    Multiply Q in U by left vectors bidiagonalizing R */
/*                    (Workspace: need 3*N + M, prefer 3*N + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"Q", (char *)"R", (char *)"N", m, n, n, &a[a_offset], lda, &
                                work[itauq], &u[u_offset], ldu, &work[iwork],
                                &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
                                ;
                        iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U */
/*                    (Workspace: need BDSPAC) */

                        dbdsqr_((char *)"U", n, &c__0, m, &c__0, &s[1], &work[ie],
                                dum, &c__1, &u[u_offset], ldu, dum, &c__1, &
                                work[iwork], info, (ftnlen)1);

                    }

                } else if (wntvo) {

/*                 Path 5 (M much larger than N, JOBU='S', JOBVT='O') */
/*                 N left singular vectors to be computed in U and */
/*                 N right singular vectors to be overwritten on A */

/* Computing MAX */
                    i__2 = *n << 2;
                    if (*lwork >= (*n << 1) * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

                        iu = 1;
                        if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N */

                            ldwrku = *lda;
                            ir = iu + ldwrku * *n;
                            ldwrkr = *lda;
                        } else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is N by N */

                            ldwrku = *lda;
                            ir = iu + ldwrku * *n;
                            ldwrkr = *n;
                        } else {

/*                       WORK(IU) is N by N and WORK(IR) is N by N */

                            ldwrku = *n;
                            ir = iu + ldwrku * *n;
                            ldwrkr = *n;
                        }
                        itau = ir + ldwrkr * *n;
                        iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (Workspace: need 2*N*N + 2*N, prefer 2*N*N + N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

                        dlacpy_((char *)"U", n, n, &a[a_offset], lda, &work[iu], &
                                ldwrku, (ftnlen)1);
                        i__2 = *n - 1;
                        i__3 = *n - 1;
                        dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &work[iu +
                                1], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need 2*N*N + 2*N, prefer 2*N*N + N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &
                                work[iwork], &i__2, &ierr);
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (Workspace: need 2*N*N + 4*N, */
/*                                prefer 2*N*N+3*N+2*N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &
                                work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);
                        dlacpy_((char *)"U", n, n, &work[iu], &ldwrku, &work[ir], &
                                ldwrkr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need 2*N*N + 4*N, prefer 2*N*N + 3*N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
                                , &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need 2*N*N + 4*N-1, */
/*                                prefer 2*N*N+3*N+(N-1)*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", n, n, n, &work[ir], &ldwrkr, &work[itaup]
                                , &work[iwork], &i__2, &ierr, (ftnlen)1);
                        iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in WORK(IR) */
/*                    (Workspace: need 2*N*N + BDSPAC) */

                        dbdsqr_((char *)"U", n, n, n, &c__0, &s[1], &work[ie], &work[
                                ir], &ldwrkr, &work[iu], &ldwrku, dum, &c__1,
                                &work[iwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IU), storing result in U */
/*                    (Workspace: need N*N) */

                        dgemm_((char *)"N", (char *)"N", m, n, n, &c_b79, &a[a_offset], lda, &
                                work[iu], &ldwrku, &c_b57, &u[u_offset], ldu,
                                (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of R to A */
/*                    (Workspace: need N*N) */

                        dlacpy_((char *)"F", n, n, &work[ir], &ldwrkr, &a[a_offset],
                                lda, (ftnlen)1);

                    } else {

/*                    Insufficient workspace for a fast algorithm */

                        itau = 1;
                        iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, n, &a[a_offset], lda, &u[u_offset],
                                ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need 2*N, prefer N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
                                work[iwork], &i__2, &ierr);
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;

/*                    Zero out below R in A */

                        if (*n > 1) {
                            i__2 = *n - 1;
                            i__3 = *n - 1;
                            dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &a[
                                    a_dim1 + 2], lda, (ftnlen)1);
                        }

/*                    Bidiagonalize R in A */
/*                    (Workspace: need 4*N, prefer 3*N + 2*N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &
                                work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);

/*                    Multiply Q in U by left vectors bidiagonalizing R */
/*                    (Workspace: need 3*N + M, prefer 3*N + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"Q", (char *)"R", (char *)"N", m, n, n, &a[a_offset], lda, &
                                work[itauq], &u[u_offset], ldu, &work[iwork],
                                &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
                                ;

/*                    Generate right vectors bidiagonalizing R in A */
/*                    (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", n, n, n, &a[a_offset], lda, &work[itaup],
                                 &work[iwork], &i__2, &ierr, (ftnlen)1);
                        iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in A */
/*                    (Workspace: need BDSPAC) */

                        dbdsqr_((char *)"U", n, n, m, &c__0, &s[1], &work[ie], &a[
                                a_offset], lda, &u[u_offset], ldu, dum, &c__1,
                                 &work[iwork], info, (ftnlen)1);

                    }

                } else if (wntvas) {

/*                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' */
/*                         or 'A') */
/*                 N left singular vectors to be computed in U and */
/*                 N right singular vectors to be computed in VT */

/* Computing MAX */
                    i__2 = *n << 2;
                    if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

                        iu = 1;
                        if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

                            ldwrku = *lda;
                        } else {

/*                       WORK(IU) is N by N */

                            ldwrku = *n;
                        }
                        itau = iu + ldwrku * *n;
                        iwork = itau + *n;

/*                    Compute A=Q*R */
/*                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

                        dlacpy_((char *)"U", n, n, &a[a_offset], lda, &work[iu], &
                                ldwrku, (ftnlen)1);
                        i__2 = *n - 1;
                        i__3 = *n - 1;
                        dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &work[iu +
                                1], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, n, n, &a[a_offset], lda, &work[itau], &
                                work[iwork], &i__2, &ierr);
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to VT */
/*                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &
                                work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);
                        dlacpy_((char *)"U", n, n, &work[iu], &ldwrku, &vt[vt_offset],
                                 ldvt, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
                                , &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (Workspace: need N*N + 4*N-1, */
/*                                prefer N*N+3*N+(N-1)*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", n, n, n, &vt[vt_offset], ldvt, &work[
                                itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
                                ;
                        iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in VT */
/*                    (Workspace: need N*N + BDSPAC) */

                        dbdsqr_((char *)"U", n, n, n, &c__0, &s[1], &work[ie], &vt[
                                vt_offset], ldvt, &work[iu], &ldwrku, dum, &
                                c__1, &work[iwork], info, (ftnlen)1);

/*                    Multiply Q in A by left singular vectors of R in */
/*                    WORK(IU), storing result in U */
/*                    (Workspace: need N*N) */

                        dgemm_((char *)"N", (char *)"N", m, n, n, &c_b79, &a[a_offset], lda, &
                                work[iu], &ldwrku, &c_b57, &u[u_offset], ldu,
                                (ftnlen)1, (ftnlen)1);

                    } else {

/*                    Insufficient workspace for a fast algorithm */

                        itau = 1;
                        iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, n, &a[a_offset], lda, &u[u_offset],
                                ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need 2*N, prefer N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, n, n, &u[u_offset], ldu, &work[itau], &
                                work[iwork], &i__2, &ierr);

/*                    Copy R to VT, zeroing out below it */

                        dlacpy_((char *)"U", n, n, &a[a_offset], lda, &vt[vt_offset],
                                ldvt, (ftnlen)1);
                        if (*n > 1) {
                            i__2 = *n - 1;
                            i__3 = *n - 1;
                            dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &vt[
                                    vt_dim1 + 2], ldvt, (ftnlen)1);
                        }
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;

/*                    Bidiagonalize R in VT */
/*                    (Workspace: need 4*N, prefer 3*N + 2*N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie],
                                &work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in VT */
/*                    (Workspace: need 3*N + M, prefer 3*N + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"Q", (char *)"R", (char *)"N", m, n, n, &vt[vt_offset], ldvt,
                                &work[itauq], &u[u_offset], ldu, &work[iwork],
                                 &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
                                1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", n, n, n, &vt[vt_offset], ldvt, &work[
                                itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
                                ;
                        iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

                        dbdsqr_((char *)"U", n, n, m, &c__0, &s[1], &work[ie], &vt[
                                vt_offset], ldvt, &u[u_offset], ldu, dum, &
                                c__1, &work[iwork], info, (ftnlen)1);

                    }

                }

            } else if (wntua) {

                if (wntvn) {

/*                 Path 7 (M much larger than N, JOBU='A', JOBVT='N') */
/*                 M left singular vectors to be computed in U and */
/*                 no right singular vectors to be computed */

/* Computing MAX */
                    i__2 = *n + *m, i__3 = *n << 2, i__2 = max(i__2,i__3);
                    if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

                        ir = 1;
                        if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IR) is LDA by N */

                            ldwrkr = *lda;
                        } else {

/*                       WORK(IR) is N by N */

                            ldwrkr = *n;
                        }
                        itau = ir + ldwrkr * *n;
                        iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, n, &a[a_offset], lda, &u[u_offset],
                                ldu, (ftnlen)1);

/*                    Copy R to WORK(IR), zeroing out below it */

                        dlacpy_((char *)"U", n, n, &a[a_offset], lda, &work[ir], &
                                ldwrkr, (ftnlen)1);
                        i__2 = *n - 1;
                        i__3 = *n - 1;
                        dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &work[ir +
                                1], &ldwrkr, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N*N + N + M, prefer N*N + N + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
                                work[iwork], &i__2, &ierr);
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IR) */
/*                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &work[ir], &ldwrkr, &s[1], &work[ie], &
                                work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", n, n, n, &work[ir], &ldwrkr, &work[itauq]
                                , &work[iwork], &i__2, &ierr, (ftnlen)1);
                        iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IR) */
/*                    (Workspace: need N*N + BDSPAC) */

                        dbdsqr_((char *)"U", n, &c__0, n, &c__0, &s[1], &work[ie],
                                dum, &c__1, &work[ir], &ldwrkr, dum, &c__1, &
                                work[iwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IR), storing result in A */
/*                    (Workspace: need N*N) */

                        dgemm_((char *)"N", (char *)"N", m, n, n, &c_b79, &u[u_offset], ldu, &
                                work[ir], &ldwrkr, &c_b57, &a[a_offset], lda,
                                (ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

                        dlacpy_((char *)"F", m, n, &a[a_offset], lda, &u[u_offset],
                                ldu, (ftnlen)1);

                    } else {

/*                    Insufficient workspace for a fast algorithm */

                        itau = 1;
                        iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, n, &a[a_offset], lda, &u[u_offset],
                                ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N + M, prefer N + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
                                work[iwork], &i__2, &ierr);
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;

/*                    Zero out below R in A */

                        if (*n > 1) {
                            i__2 = *n - 1;
                            i__3 = *n - 1;
                            dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &a[
                                    a_dim1 + 2], lda, (ftnlen)1);
                        }

/*                    Bidiagonalize R in A */
/*                    (Workspace: need 4*N, prefer 3*N + 2*N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &
                                work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in A */
/*                    (Workspace: need 3*N + M, prefer 3*N + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"Q", (char *)"R", (char *)"N", m, n, n, &a[a_offset], lda, &
                                work[itauq], &u[u_offset], ldu, &work[iwork],
                                &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
                                ;
                        iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U */
/*                    (Workspace: need BDSPAC) */

                        dbdsqr_((char *)"U", n, &c__0, m, &c__0, &s[1], &work[ie],
                                dum, &c__1, &u[u_offset], ldu, dum, &c__1, &
                                work[iwork], info, (ftnlen)1);

                    }

                } else if (wntvo) {

/*                 Path 8 (M much larger than N, JOBU='A', JOBVT='O') */
/*                 M left singular vectors to be computed in U and */
/*                 N right singular vectors to be overwritten on A */

/* Computing MAX */
                    i__2 = *n + *m, i__3 = *n << 2, i__2 = max(i__2,i__3);
                    if (*lwork >= (*n << 1) * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

                        iu = 1;
                        if (*lwork >= wrkbl + (*lda << 1) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N */

                            ldwrku = *lda;
                            ir = iu + ldwrku * *n;
                            ldwrkr = *lda;
                        } else if (*lwork >= wrkbl + (*lda + *n) * *n) {

/*                       WORK(IU) is LDA by N and WORK(IR) is N by N */

                            ldwrku = *lda;
                            ir = iu + ldwrku * *n;
                            ldwrkr = *n;
                        } else {

/*                       WORK(IU) is N by N and WORK(IR) is N by N */

                            ldwrku = *n;
                            ir = iu + ldwrku * *n;
                            ldwrkr = *n;
                        }
                        itau = ir + ldwrkr * *n;
                        iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N*N + 2*N, prefer 2*N*N + N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, n, &a[a_offset], lda, &u[u_offset],
                                ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need 2*N*N + N + M, prefer 2*N*N + N + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
                                work[iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

                        dlacpy_((char *)"U", n, n, &a[a_offset], lda, &work[iu], &
                                ldwrku, (ftnlen)1);
                        i__2 = *n - 1;
                        i__3 = *n - 1;
                        dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &work[iu +
                                1], &ldwrku, (ftnlen)1);
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (Workspace: need 2*N*N + 4*N, */
/*                                prefer 2*N*N+3*N+2*N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &
                                work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);
                        dlacpy_((char *)"U", n, n, &work[iu], &ldwrku, &work[ir], &
                                ldwrkr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need 2*N*N + 4*N, prefer 2*N*N + 3*N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
                                , &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need 2*N*N + 4*N-1, */
/*                                prefer 2*N*N+3*N+(N-1)*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", n, n, n, &work[ir], &ldwrkr, &work[itaup]
                                , &work[iwork], &i__2, &ierr, (ftnlen)1);
                        iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in WORK(IR) */
/*                    (Workspace: need 2*N*N + BDSPAC) */

                        dbdsqr_((char *)"U", n, n, n, &c__0, &s[1], &work[ie], &work[
                                ir], &ldwrkr, &work[iu], &ldwrku, dum, &c__1,
                                &work[iwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IU), storing result in A */
/*                    (Workspace: need N*N) */

                        dgemm_((char *)"N", (char *)"N", m, n, n, &c_b79, &u[u_offset], ldu, &
                                work[iu], &ldwrku, &c_b57, &a[a_offset], lda,
                                (ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

                        dlacpy_((char *)"F", m, n, &a[a_offset], lda, &u[u_offset],
                                ldu, (ftnlen)1);

/*                    Copy right singular vectors of R from WORK(IR) to A */

                        dlacpy_((char *)"F", n, n, &work[ir], &ldwrkr, &a[a_offset],
                                lda, (ftnlen)1);

                    } else {

/*                    Insufficient workspace for a fast algorithm */

                        itau = 1;
                        iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, n, &a[a_offset], lda, &u[u_offset],
                                ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N + M, prefer N + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
                                work[iwork], &i__2, &ierr);
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;

/*                    Zero out below R in A */

                        if (*n > 1) {
                            i__2 = *n - 1;
                            i__3 = *n - 1;
                            dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &a[
                                    a_dim1 + 2], lda, (ftnlen)1);
                        }

/*                    Bidiagonalize R in A */
/*                    (Workspace: need 4*N, prefer 3*N + 2*N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &a[a_offset], lda, &s[1], &work[ie], &
                                work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in A */
/*                    (Workspace: need 3*N + M, prefer 3*N + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"Q", (char *)"R", (char *)"N", m, n, n, &a[a_offset], lda, &
                                work[itauq], &u[u_offset], ldu, &work[iwork],
                                &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1)
                                ;

/*                    Generate right bidiagonalizing vectors in A */
/*                    (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", n, n, n, &a[a_offset], lda, &work[itaup],
                                 &work[iwork], &i__2, &ierr, (ftnlen)1);
                        iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in A */
/*                    (Workspace: need BDSPAC) */

                        dbdsqr_((char *)"U", n, n, m, &c__0, &s[1], &work[ie], &a[
                                a_offset], lda, &u[u_offset], ldu, dum, &c__1,
                                 &work[iwork], info, (ftnlen)1);

                    }

                } else if (wntvas) {

/*                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' */
/*                         or 'A') */
/*                 M left singular vectors to be computed in U and */
/*                 N right singular vectors to be computed in VT */

/* Computing MAX */
                    i__2 = *n + *m, i__3 = *n << 2, i__2 = max(i__2,i__3);
                    if (*lwork >= *n * *n + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

                        iu = 1;
                        if (*lwork >= wrkbl + *lda * *n) {

/*                       WORK(IU) is LDA by N */

                            ldwrku = *lda;
                        } else {

/*                       WORK(IU) is N by N */

                            ldwrku = *n;
                        }
                        itau = iu + ldwrku * *n;
                        iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, n, &a[a_offset], lda, &u[u_offset],
                                ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N*N + N + M, prefer N*N + N + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
                                work[iwork], &i__2, &ierr);

/*                    Copy R to WORK(IU), zeroing out below it */

                        dlacpy_((char *)"U", n, n, &a[a_offset], lda, &work[iu], &
                                ldwrku, (ftnlen)1);
                        i__2 = *n - 1;
                        i__3 = *n - 1;
                        dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &work[iu +
                                1], &ldwrku, (ftnlen)1);
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;

/*                    Bidiagonalize R in WORK(IU), copying result to VT */
/*                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &work[iu], &ldwrku, &s[1], &work[ie], &
                                work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);
                        dlacpy_((char *)"U", n, n, &work[iu], &ldwrku, &vt[vt_offset],
                                 ldvt, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", n, n, n, &work[iu], &ldwrku, &work[itauq]
                                , &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (Workspace: need N*N + 4*N-1, */
/*                                prefer N*N+3*N+(N-1)*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", n, n, n, &vt[vt_offset], ldvt, &work[
                                itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
                                ;
                        iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of R in WORK(IU) and computing */
/*                    right singular vectors of R in VT */
/*                    (Workspace: need N*N + BDSPAC) */

                        dbdsqr_((char *)"U", n, n, n, &c__0, &s[1], &work[ie], &vt[
                                vt_offset], ldvt, &work[iu], &ldwrku, dum, &
                                c__1, &work[iwork], info, (ftnlen)1);

/*                    Multiply Q in U by left singular vectors of R in */
/*                    WORK(IU), storing result in A */
/*                    (Workspace: need N*N) */

                        dgemm_((char *)"N", (char *)"N", m, n, n, &c_b79, &u[u_offset], ldu, &
                                work[iu], &ldwrku, &c_b57, &a[a_offset], lda,
                                (ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of A from A to U */

                        dlacpy_((char *)"F", m, n, &a[a_offset], lda, &u[u_offset],
                                ldu, (ftnlen)1);

                    } else {

/*                    Insufficient workspace for a fast algorithm */

                        itau = 1;
                        iwork = itau + *n;

/*                    Compute A=Q*R, copying result to U */
/*                    (Workspace: need 2*N, prefer N + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);
                        dlacpy_((char *)"L", m, n, &a[a_offset], lda, &u[u_offset],
                                ldu, (ftnlen)1);

/*                    Generate Q in U */
/*                    (Workspace: need N + M, prefer N + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgqr_(m, m, n, &u[u_offset], ldu, &work[itau], &
                                work[iwork], &i__2, &ierr);

/*                    Copy R from A to VT, zeroing out below it */

                        dlacpy_((char *)"U", n, n, &a[a_offset], lda, &vt[vt_offset],
                                ldvt, (ftnlen)1);
                        if (*n > 1) {
                            i__2 = *n - 1;
                            i__3 = *n - 1;
                            dlaset_((char *)"L", &i__2, &i__3, &c_b57, &c_b57, &vt[
                                    vt_dim1 + 2], ldvt, (ftnlen)1);
                        }
                        ie = itau;
                        itauq = ie + *n;
                        itaup = itauq + *n;
                        iwork = itaup + *n;

/*                    Bidiagonalize R in VT */
/*                    (Workspace: need 4*N, prefer 3*N + 2*N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(n, n, &vt[vt_offset], ldvt, &s[1], &work[ie],
                                &work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);

/*                    Multiply Q in U by left bidiagonalizing vectors */
/*                    in VT */
/*                    (Workspace: need 3*N + M, prefer 3*N + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"Q", (char *)"R", (char *)"N", m, n, n, &vt[vt_offset], ldvt,
                                &work[itauq], &u[u_offset], ldu, &work[iwork],
                                 &i__2, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)
                                1);

/*                    Generate right bidiagonalizing vectors in VT */
/*                    (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", n, n, n, &vt[vt_offset], ldvt, &work[
                                itaup], &work[iwork], &i__2, &ierr, (ftnlen)1)
                                ;
                        iwork = ie + *n;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

                        dbdsqr_((char *)"U", n, n, m, &c__0, &s[1], &work[ie], &vt[
                                vt_offset], ldvt, &u[u_offset], ldu, dum, &
                                c__1, &work[iwork], info, (ftnlen)1);

                    }

                }

            }

        } else {

/*           M .LT. MNTHR */

/*           Path 10 (M at least N, but not much larger) */
/*           Reduce to bidiagonal form without QR decomposition */

            ie = 1;
            itauq = ie + *n;
            itaup = itauq + *n;
            iwork = itaup + *n;

/*           Bidiagonalize A */
/*           (Workspace: need 3*N + M, prefer 3*N + (M + N)*NB) */

            i__2 = *lwork - iwork + 1;
            dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
                    work[itaup], &work[iwork], &i__2, &ierr);
            if (wntuas) {

/*              If left singular vectors desired in U, copy result to U */
/*              and generate left bidiagonalizing vectors in U */
/*              (Workspace: need 3*N + NCU, prefer 3*N + NCU*NB) */

                dlacpy_((char *)"L", m, n, &a[a_offset], lda, &u[u_offset], ldu, (
                        ftnlen)1);
                if (wntus) {
                    ncu = *n;
                }
                if (wntua) {
                    ncu = *m;
                }
                i__2 = *lwork - iwork + 1;
                dorgbr_((char *)"Q", m, &ncu, n, &u[u_offset], ldu, &work[itauq], &
                        work[iwork], &i__2, &ierr, (ftnlen)1);
            }
            if (wntvas) {

/*              If right singular vectors desired in VT, copy result to */
/*              VT and generate right bidiagonalizing vectors in VT */
/*              (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB) */

                dlacpy_((char *)"U", n, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
                        ftnlen)1);
                i__2 = *lwork - iwork + 1;
                dorgbr_((char *)"P", n, n, n, &vt[vt_offset], ldvt, &work[itaup], &
                        work[iwork], &i__2, &ierr, (ftnlen)1);
            }
            if (wntuo) {

/*              If left singular vectors desired in A, generate left */
/*              bidiagonalizing vectors in A */
/*              (Workspace: need 4*N, prefer 3*N + N*NB) */

                i__2 = *lwork - iwork + 1;
                dorgbr_((char *)"Q", m, n, n, &a[a_offset], lda, &work[itauq], &work[
                        iwork], &i__2, &ierr, (ftnlen)1);
            }
            if (wntvo) {

/*              If right singular vectors desired in A, generate right */
/*              bidiagonalizing vectors in A */
/*              (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB) */

                i__2 = *lwork - iwork + 1;
                dorgbr_((char *)"P", n, n, n, &a[a_offset], lda, &work[itaup], &work[
                        iwork], &i__2, &ierr, (ftnlen)1);
            }
            iwork = ie + *n;
            if (wntuas || wntuo) {
                nru = *m;
            }
            if (wntun) {
                nru = 0;
            }
            if (wntvas || wntvo) {
                ncvt = *n;
            }
            if (wntvn) {
                ncvt = 0;
            }
            if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in VT */
/*              (Workspace: need BDSPAC) */

                dbdsqr_((char *)"U", n, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[
                        vt_offset], ldvt, &u[u_offset], ldu, dum, &c__1, &
                        work[iwork], info, (ftnlen)1);
            } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in A */
/*              (Workspace: need BDSPAC) */

                dbdsqr_((char *)"U", n, &ncvt, &nru, &c__0, &s[1], &work[ie], &a[
                        a_offset], lda, &u[u_offset], ldu, dum, &c__1, &work[
                        iwork], info, (ftnlen)1);
            } else {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in A and computing right singular */
/*              vectors in VT */
/*              (Workspace: need BDSPAC) */

                dbdsqr_((char *)"U", n, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[
                        vt_offset], ldvt, &a[a_offset], lda, dum, &c__1, &
                        work[iwork], info, (ftnlen)1);
            }

        }

    } else {

/*        A has more columns than rows. If A has sufficiently more */
/*        columns than rows, first reduce using the LQ decomposition (if */
/*        sufficient workspace available) */

        if (*n >= mnthr) {

            if (wntvn) {

/*              Path 1t(N much larger than M, JOBVT='N') */
/*              No right singular vectors to be computed */

                itau = 1;
                iwork = itau + *m;

/*              Compute A=L*Q */
/*              (Workspace: need 2*M, prefer M + M*NB) */

                i__2 = *lwork - iwork + 1;
                dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &
                        i__2, &ierr);

/*              Zero out above L */

                i__2 = *m - 1;
                i__3 = *m - 1;
                dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1 << 1) +
                        1], lda, (ftnlen)1);
                ie = 1;
                itauq = ie + *m;
                itaup = itauq + *m;
                iwork = itaup + *m;

/*              Bidiagonalize L in A */
/*              (Workspace: need 4*M, prefer 3*M + 2*M*NB) */

                i__2 = *lwork - iwork + 1;
                dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &work[
                        itauq], &work[itaup], &work[iwork], &i__2, &ierr);
                if (wntuo || wntuas) {

/*                 If left singular vectors desired, generate Q */
/*                 (Workspace: need 4*M, prefer 3*M + M*NB) */

                    i__2 = *lwork - iwork + 1;
                    dorgbr_((char *)"Q", m, m, m, &a[a_offset], lda, &work[itauq], &
                            work[iwork], &i__2, &ierr, (ftnlen)1);
                }
                iwork = ie + *m;
                nru = 0;
                if (wntuo || wntuas) {
                    nru = *m;
                }

/*              Perform bidiagonal QR iteration, computing left singular */
/*              vectors of A in A if desired */
/*              (Workspace: need BDSPAC) */

                dbdsqr_((char *)"U", m, &c__0, &nru, &c__0, &s[1], &work[ie], dum, &
                        c__1, &a[a_offset], lda, dum, &c__1, &work[iwork],
                        info, (ftnlen)1);

/*              If left singular vectors desired in U, copy them there */

                if (wntuas) {
                    dlacpy_((char *)"F", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
                            ftnlen)1);
                }

            } else if (wntvo && wntun) {

/*              Path 2t(N much larger than M, JOBU='N', JOBVT='O') */
/*              M right singular vectors to be overwritten on A and */
/*              no left singular vectors to be computed */

/* Computing MAX */
                i__2 = *m << 2;
                if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                 Sufficient workspace for a fast algorithm */

                    ir = 1;
/* Computing MAX */
                    i__2 = wrkbl, i__3 = *lda * *n + *m;
                    if (*lwork >= max(i__2,i__3) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M */

                        ldwrku = *lda;
                        chunk = *n;
                        ldwrkr = *lda;
                    } else /* if(complicated condition) */ {
/* Computing MAX */
                        i__2 = wrkbl, i__3 = *lda * *n + *m;
                        if (*lwork >= max(i__2,i__3) + *m * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is M by M */

                            ldwrku = *lda;
                            chunk = *n;
                            ldwrkr = *m;
                        } else {

/*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M */

                            ldwrku = *m;
                            chunk = (*lwork - *m * *m - *m) / *m;
                            ldwrkr = *m;
                        }
                    }
                    itau = ir + ldwrkr * *m;
                    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (Workspace: need M*M + 2*M, prefer M*M + M + M*NB) */

                    i__2 = *lwork - iwork + 1;
                    dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
                            , &i__2, &ierr);

/*                 Copy L to WORK(IR) and zero out above it */

                    dlacpy_((char *)"L", m, m, &a[a_offset], lda, &work[ir], &ldwrkr,
                            (ftnlen)1);
                    i__2 = *m - 1;
                    i__3 = *m - 1;
                    dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &work[ir +
                            ldwrkr], &ldwrkr, (ftnlen)1);

/*                 Generate Q in A */
/*                 (Workspace: need M*M + 2*M, prefer M*M + M + M*NB) */

                    i__2 = *lwork - iwork + 1;
                    dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
                            iwork], &i__2, &ierr);
                    ie = itau;
                    itauq = ie + *m;
                    itaup = itauq + *m;
                    iwork = itaup + *m;

/*                 Bidiagonalize L in WORK(IR) */
/*                 (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB) */

                    i__2 = *lwork - iwork + 1;
                    dgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &work[ie], &work[
                            itauq], &work[itaup], &work[iwork], &i__2, &ierr);

/*                 Generate right vectors bidiagonalizing L */
/*                 (Workspace: need M*M + 4*M-1, prefer M*M + 3*M + (M-1)*NB) */

                    i__2 = *lwork - iwork + 1;
                    dorgbr_((char *)"P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &
                            work[iwork], &i__2, &ierr, (ftnlen)1);
                    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing right */
/*                 singular vectors of L in WORK(IR) */
/*                 (Workspace: need M*M + BDSPAC) */

                    dbdsqr_((char *)"U", m, m, &c__0, &c__0, &s[1], &work[ie], &work[
                            ir], &ldwrkr, dum, &c__1, dum, &c__1, &work[iwork]
                            , info, (ftnlen)1);
                    iu = ie + *m;

/*                 Multiply right singular vectors of L in WORK(IR) by Q */
/*                 in A, storing result in WORK(IU) and copying to A */
/*                 (Workspace: need M*M + 2*M, prefer M*M + M*N + M) */

                    i__2 = *n;
                    i__3 = chunk;
                    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
                             i__3) {
/* Computing MIN */
                        i__4 = *n - i__ + 1;
                        blk = min(i__4,chunk);
                        dgemm_((char *)"N", (char *)"N", m, &blk, m, &c_b79, &work[ir], &
                                ldwrkr, &a[i__ * a_dim1 + 1], lda, &c_b57, &
                                work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
                        dlacpy_((char *)"F", m, &blk, &work[iu], &ldwrku, &a[i__ *
                                a_dim1 + 1], lda, (ftnlen)1);
/* L30: */
                    }

                } else {

/*                 Insufficient workspace for a fast algorithm */

                    ie = 1;
                    itauq = ie + *m;
                    itaup = itauq + *m;
                    iwork = itaup + *m;

/*                 Bidiagonalize A */
/*                 (Workspace: need 3*M + N, prefer 3*M + (M + N)*NB) */

                    i__3 = *lwork - iwork + 1;
                    dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[
                            itauq], &work[itaup], &work[iwork], &i__3, &ierr);

/*                 Generate right vectors bidiagonalizing A */
/*                 (Workspace: need 4*M, prefer 3*M + M*NB) */

                    i__3 = *lwork - iwork + 1;
                    dorgbr_((char *)"P", m, n, m, &a[a_offset], lda, &work[itaup], &
                            work[iwork], &i__3, &ierr, (ftnlen)1);
                    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing right */
/*                 singular vectors of A in A */
/*                 (Workspace: need BDSPAC) */

                    dbdsqr_((char *)"L", m, n, &c__0, &c__0, &s[1], &work[ie], &a[
                            a_offset], lda, dum, &c__1, dum, &c__1, &work[
                            iwork], info, (ftnlen)1);

                }

            } else if (wntvo && wntuas) {

/*              Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O') */
/*              M right singular vectors to be overwritten on A and */
/*              M left singular vectors to be computed in U */

/* Computing MAX */
                i__3 = *m << 2;
                if (*lwork >= *m * *m + max(i__3,bdspac)) {

/*                 Sufficient workspace for a fast algorithm */

                    ir = 1;
/* Computing MAX */
                    i__3 = wrkbl, i__2 = *lda * *n + *m;
                    if (*lwork >= max(i__3,i__2) + *lda * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M */

                        ldwrku = *lda;
                        chunk = *n;
                        ldwrkr = *lda;
                    } else /* if(complicated condition) */ {
/* Computing MAX */
                        i__3 = wrkbl, i__2 = *lda * *n + *m;
                        if (*lwork >= max(i__3,i__2) + *m * *m) {

/*                    WORK(IU) is LDA by N and WORK(IR) is M by M */

                            ldwrku = *lda;
                            chunk = *n;
                            ldwrkr = *m;
                        } else {

/*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M */

                            ldwrku = *m;
                            chunk = (*lwork - *m * *m - *m) / *m;
                            ldwrkr = *m;
                        }
                    }
                    itau = ir + ldwrkr * *m;
                    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (Workspace: need M*M + 2*M, prefer M*M + M + M*NB) */

                    i__3 = *lwork - iwork + 1;
                    dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
                            , &i__3, &ierr);

/*                 Copy L to U, zeroing about above it */

                    dlacpy_((char *)"L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
                            ftnlen)1);
                    i__3 = *m - 1;
                    i__2 = *m - 1;
                    dlaset_((char *)"U", &i__3, &i__2, &c_b57, &c_b57, &u[(u_dim1 <<
                            1) + 1], ldu, (ftnlen)1);

/*                 Generate Q in A */
/*                 (Workspace: need M*M + 2*M, prefer M*M + M + M*NB) */

                    i__3 = *lwork - iwork + 1;
                    dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
                            iwork], &i__3, &ierr);
                    ie = itau;
                    itauq = ie + *m;
                    itaup = itauq + *m;
                    iwork = itaup + *m;

/*                 Bidiagonalize L in U, copying result to WORK(IR) */
/*                 (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB) */

                    i__3 = *lwork - iwork + 1;
                    dgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &work[
                            itauq], &work[itaup], &work[iwork], &i__3, &ierr);
                    dlacpy_((char *)"U", m, m, &u[u_offset], ldu, &work[ir], &ldwrkr,
                            (ftnlen)1);

/*                 Generate right vectors bidiagonalizing L in WORK(IR) */
/*                 (Workspace: need M*M + 4*M-1, prefer M*M + 3*M + (M-1)*NB) */

                    i__3 = *lwork - iwork + 1;
                    dorgbr_((char *)"P", m, m, m, &work[ir], &ldwrkr, &work[itaup], &
                            work[iwork], &i__3, &ierr, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing L in U */
/*                 (Workspace: need M*M + 4*M, prefer M*M + 3*M + M*NB) */

                    i__3 = *lwork - iwork + 1;
                    dorgbr_((char *)"Q", m, m, m, &u[u_offset], ldu, &work[itauq], &
                            work[iwork], &i__3, &ierr, (ftnlen)1);
                    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of L in U, and computing right */
/*                 singular vectors of L in WORK(IR) */
/*                 (Workspace: need M*M + BDSPAC) */

                    dbdsqr_((char *)"U", m, m, m, &c__0, &s[1], &work[ie], &work[ir],
                            &ldwrkr, &u[u_offset], ldu, dum, &c__1, &work[
                            iwork], info, (ftnlen)1);
                    iu = ie + *m;

/*                 Multiply right singular vectors of L in WORK(IR) by Q */
/*                 in A, storing result in WORK(IU) and copying to A */
/*                 (Workspace: need M*M + 2*M, prefer M*M + M*N + M)) */

                    i__3 = *n;
                    i__2 = chunk;
                    for (i__ = 1; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ +=
                             i__2) {
/* Computing MIN */
                        i__4 = *n - i__ + 1;
                        blk = min(i__4,chunk);
                        dgemm_((char *)"N", (char *)"N", m, &blk, m, &c_b79, &work[ir], &
                                ldwrkr, &a[i__ * a_dim1 + 1], lda, &c_b57, &
                                work[iu], &ldwrku, (ftnlen)1, (ftnlen)1);
                        dlacpy_((char *)"F", m, &blk, &work[iu], &ldwrku, &a[i__ *
                                a_dim1 + 1], lda, (ftnlen)1);
/* L40: */
                    }

                } else {

/*                 Insufficient workspace for a fast algorithm */

                    itau = 1;
                    iwork = itau + *m;

/*                 Compute A=L*Q */
/*                 (Workspace: need 2*M, prefer M + M*NB) */

                    i__2 = *lwork - iwork + 1;
                    dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork]
                            , &i__2, &ierr);

/*                 Copy L to U, zeroing out above it */

                    dlacpy_((char *)"L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
                            ftnlen)1);
                    i__2 = *m - 1;
                    i__3 = *m - 1;
                    dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &u[(u_dim1 <<
                            1) + 1], ldu, (ftnlen)1);

/*                 Generate Q in A */
/*                 (Workspace: need 2*M, prefer M + M*NB) */

                    i__2 = *lwork - iwork + 1;
                    dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &work[
                            iwork], &i__2, &ierr);
                    ie = itau;
                    itauq = ie + *m;
                    itaup = itauq + *m;
                    iwork = itaup + *m;

/*                 Bidiagonalize L in U */
/*                 (Workspace: need 4*M, prefer 3*M + 2*M*NB) */

                    i__2 = *lwork - iwork + 1;
                    dgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &work[
                            itauq], &work[itaup], &work[iwork], &i__2, &ierr);

/*                 Multiply right vectors bidiagonalizing L by Q in A */
/*                 (Workspace: need 3*M + N, prefer 3*M + N*NB) */

                    i__2 = *lwork - iwork + 1;
                    dormbr_((char *)"P", (char *)"L", (char *)"T", m, n, m, &u[u_offset], ldu, &work[
                            itaup], &a[a_offset], lda, &work[iwork], &i__2, &
                            ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*                 Generate left vectors bidiagonalizing L in U */
/*                 (Workspace: need 4*M, prefer 3*M + M*NB) */

                    i__2 = *lwork - iwork + 1;
                    dorgbr_((char *)"Q", m, m, m, &u[u_offset], ldu, &work[itauq], &
                            work[iwork], &i__2, &ierr, (ftnlen)1);
                    iwork = ie + *m;

/*                 Perform bidiagonal QR iteration, computing left */
/*                 singular vectors of A in U and computing right */
/*                 singular vectors of A in A */
/*                 (Workspace: need BDSPAC) */

                    dbdsqr_((char *)"U", m, n, m, &c__0, &s[1], &work[ie], &a[
                            a_offset], lda, &u[u_offset], ldu, dum, &c__1, &
                            work[iwork], info, (ftnlen)1);

                }

            } else if (wntvs) {

                if (wntun) {

/*                 Path 4t(N much larger than M, JOBU='N', JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 no left singular vectors to be computed */

/* Computing MAX */
                    i__2 = *m << 2;
                    if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

                        ir = 1;
                        if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

                            ldwrkr = *lda;
                        } else {

/*                       WORK(IR) is M by M */

                            ldwrkr = *m;
                        }
                        itau = ir + ldwrkr * *m;
                        iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);

/*                    Copy L to WORK(IR), zeroing out above it */

                        dlacpy_((char *)"L", m, m, &a[a_offset], lda, &work[ir], &
                                ldwrkr, (ftnlen)1);
                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &work[ir +
                                ldwrkr], &ldwrkr, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &
                                work[iwork], &i__2, &ierr);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR) */
/*                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &work[ie], &
                                work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);

/*                    Generate right vectors bidiagonalizing L in */
/*                    WORK(IR) */
/*                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + (M-1)*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", m, m, m, &work[ir], &ldwrkr, &work[itaup]
                                , &work[iwork], &i__2, &ierr, (ftnlen)1);
                        iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of L in WORK(IR) */
/*                    (Workspace: need M*M + BDSPAC) */

                        dbdsqr_((char *)"U", m, m, &c__0, &c__0, &s[1], &work[ie], &
                                work[ir], &ldwrkr, dum, &c__1, dum, &c__1, &
                                work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IR) by */
/*                    Q in A, storing result in VT */
/*                    (Workspace: need M*M) */

                        dgemm_((char *)"N", (char *)"N", m, n, m, &c_b79, &work[ir], &ldwrkr,
                                &a[a_offset], lda, &c_b57, &vt[vt_offset],
                                ldvt, (ftnlen)1, (ftnlen)1);

                    } else {

/*                    Insufficient workspace for a fast algorithm */

                        itau = 1;
                        iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (Workspace: need 2*M, prefer M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);

/*                    Copy result to VT */

                        dlacpy_((char *)"U", m, n, &a[a_offset], lda, &vt[vt_offset],
                                ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need 2*M, prefer M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
                                work[iwork], &i__2, &ierr);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;

/*                    Zero out above L in A */

                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1
                                << 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (Workspace: need 4*M, prefer 3*M + 2*M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &
                                work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);

/*                    Multiply right vectors bidiagonalizing L by Q in VT */
/*                    (Workspace: need 3*M + N, prefer 3*M + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"P", (char *)"L", (char *)"T", m, n, m, &a[a_offset], lda, &
                                work[itaup], &vt[vt_offset], ldvt, &work[
                                iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
                                ftnlen)1);
                        iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

                        dbdsqr_((char *)"U", m, n, &c__0, &c__0, &s[1], &work[ie], &
                                vt[vt_offset], ldvt, dum, &c__1, dum, &c__1, &
                                work[iwork], info, (ftnlen)1);

                    }

                } else if (wntuo) {

/*                 Path 5t(N much larger than M, JOBU='O', JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 M left singular vectors to be overwritten on A */

/* Computing MAX */
                    i__2 = *m << 2;
                    if (*lwork >= (*m << 1) * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

                        iu = 1;
                        if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M */

                            ldwrku = *lda;
                            ir = iu + ldwrku * *m;
                            ldwrkr = *lda;
                        } else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is M by M */

                            ldwrku = *lda;
                            ir = iu + ldwrku * *m;
                            ldwrkr = *m;
                        } else {

/*                       WORK(IU) is M by M and WORK(IR) is M by M */

                            ldwrku = *m;
                            ir = iu + ldwrku * *m;
                            ldwrkr = *m;
                        }
                        itau = ir + ldwrkr * *m;
                        iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (Workspace: need 2*M*M + 2*M, prefer 2*M*M + M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out below it */

                        dlacpy_((char *)"L", m, m, &a[a_offset], lda, &work[iu], &
                                ldwrku, (ftnlen)1);
                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &work[iu +
                                ldwrku], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need 2*M*M + 2*M, prefer 2*M*M + M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &
                                work[iwork], &i__2, &ierr);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (Workspace: need 2*M*M + 4*M, */
/*                                prefer 2*M*M+3*M+2*M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &
                                work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);
                        dlacpy_((char *)"L", m, m, &work[iu], &ldwrku, &work[ir], &
                                ldwrkr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need 2*M*M + 4*M-1, */
/*                                prefer 2*M*M+3*M+(M-1)*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", m, m, m, &work[iu], &ldwrku, &work[itaup]
                                , &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need 2*M*M + 4*M, prefer 2*M*M + 3*M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", m, m, m, &work[ir], &ldwrkr, &work[itauq]
                                , &work[iwork], &i__2, &ierr, (ftnlen)1);
                        iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in WORK(IR) and computing */
/*                    right singular vectors of L in WORK(IU) */
/*                    (Workspace: need 2*M*M + BDSPAC) */

                        dbdsqr_((char *)"U", m, m, m, &c__0, &s[1], &work[ie], &work[
                                iu], &ldwrku, &work[ir], &ldwrkr, dum, &c__1,
                                &work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in A, storing result in VT */
/*                    (Workspace: need M*M) */

                        dgemm_((char *)"N", (char *)"N", m, n, m, &c_b79, &work[iu], &ldwrku,
                                &a[a_offset], lda, &c_b57, &vt[vt_offset],
                                ldvt, (ftnlen)1, (ftnlen)1);

/*                    Copy left singular vectors of L to A */
/*                    (Workspace: need M*M) */

                        dlacpy_((char *)"F", m, m, &work[ir], &ldwrkr, &a[a_offset],
                                lda, (ftnlen)1);

                    } else {

/*                    Insufficient workspace for a fast algorithm */

                        itau = 1;
                        iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", m, n, &a[a_offset], lda, &vt[vt_offset],
                                ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need 2*M, prefer M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
                                work[iwork], &i__2, &ierr);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;

/*                    Zero out above L in A */

                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1
                                << 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (Workspace: need 4*M, prefer 3*M + 2*M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &
                                work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);

/*                    Multiply right vectors bidiagonalizing L by Q in VT */
/*                    (Workspace: need 3*M + N, prefer 3*M + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"P", (char *)"L", (char *)"T", m, n, m, &a[a_offset], lda, &
                                work[itaup], &vt[vt_offset], ldvt, &work[
                                iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
                                ftnlen)1);

/*                    Generate left bidiagonalizing vectors of L in A */
/*                    (Workspace: need 4*M, prefer 3*M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", m, m, m, &a[a_offset], lda, &work[itauq],
                                 &work[iwork], &i__2, &ierr, (ftnlen)1);
                        iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, compute left */
/*                    singular vectors of A in A and compute right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

                        dbdsqr_((char *)"U", m, n, m, &c__0, &s[1], &work[ie], &vt[
                                vt_offset], ldvt, &a[a_offset], lda, dum, &
                                c__1, &work[iwork], info, (ftnlen)1);

                    }

                } else if (wntuas) {

/*                 Path 6t(N much larger than M, JOBU='S' or 'A', */
/*                         JOBVT='S') */
/*                 M right singular vectors to be computed in VT and */
/*                 M left singular vectors to be computed in U */

/* Computing MAX */
                    i__2 = *m << 2;
                    if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

                        iu = 1;
                        if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by N */

                            ldwrku = *lda;
                        } else {

/*                       WORK(IU) is LDA by M */

                            ldwrku = *m;
                        }
                        itau = iu + ldwrku * *m;
                        iwork = itau + *m;

/*                    Compute A=L*Q */
/*                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

                        dlacpy_((char *)"L", m, m, &a[a_offset], lda, &work[iu], &
                                ldwrku, (ftnlen)1);
                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &work[iu +
                                ldwrku], &ldwrku, (ftnlen)1);

/*                    Generate Q in A */
/*                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorglq_(m, n, m, &a[a_offset], lda, &work[itau], &
                                work[iwork], &i__2, &ierr);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to U */
/*                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &
                                work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);
                        dlacpy_((char *)"L", m, m, &work[iu], &ldwrku, &u[u_offset],
                                ldu, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need M*M + 4*M-1, */
/*                                prefer M*M+3*M+(M-1)*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", m, m, m, &work[iu], &ldwrku, &work[itaup]
                                , &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", m, m, m, &u[u_offset], ldu, &work[itauq],
                                 &work[iwork], &i__2, &ierr, (ftnlen)1);
                        iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in U and computing right */
/*                    singular vectors of L in WORK(IU) */
/*                    (Workspace: need M*M + BDSPAC) */

                        dbdsqr_((char *)"U", m, m, m, &c__0, &s[1], &work[ie], &work[
                                iu], &ldwrku, &u[u_offset], ldu, dum, &c__1, &
                                work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in A, storing result in VT */
/*                    (Workspace: need M*M) */

                        dgemm_((char *)"N", (char *)"N", m, n, m, &c_b79, &work[iu], &ldwrku,
                                &a[a_offset], lda, &c_b57, &vt[vt_offset],
                                ldvt, (ftnlen)1, (ftnlen)1);

                    } else {

/*                    Insufficient workspace for a fast algorithm */

                        itau = 1;
                        iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", m, n, &a[a_offset], lda, &vt[vt_offset],
                                ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need 2*M, prefer M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorglq_(m, n, m, &vt[vt_offset], ldvt, &work[itau], &
                                work[iwork], &i__2, &ierr);

/*                    Copy L to U, zeroing out above it */

                        dlacpy_((char *)"L", m, m, &a[a_offset], lda, &u[u_offset],
                                ldu, (ftnlen)1);
                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &u[(u_dim1
                                << 1) + 1], ldu, (ftnlen)1);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;

/*                    Bidiagonalize L in U */
/*                    (Workspace: need 4*M, prefer 3*M + 2*M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &
                                work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in U by Q */
/*                    in VT */
/*                    (Workspace: need 3*M + N, prefer 3*M + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"P", (char *)"L", (char *)"T", m, n, m, &u[u_offset], ldu, &
                                work[itaup], &vt[vt_offset], ldvt, &work[
                                iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
                                ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (Workspace: need 4*M, prefer 3*M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", m, m, m, &u[u_offset], ldu, &work[itauq],
                                 &work[iwork], &i__2, &ierr, (ftnlen)1);
                        iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

                        dbdsqr_((char *)"U", m, n, m, &c__0, &s[1], &work[ie], &vt[
                                vt_offset], ldvt, &u[u_offset], ldu, dum, &
                                c__1, &work[iwork], info, (ftnlen)1);

                    }

                }

            } else if (wntva) {

                if (wntun) {

/*                 Path 7t(N much larger than M, JOBU='N', JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 no left singular vectors to be computed */

/* Computing MAX */
                    i__2 = *n + *m, i__3 = *m << 2, i__2 = max(i__2,i__3);
                    if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

                        ir = 1;
                        if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IR) is LDA by M */

                            ldwrkr = *lda;
                        } else {

/*                       WORK(IR) is M by M */

                            ldwrkr = *m;
                        }
                        itau = ir + ldwrkr * *m;
                        iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", m, n, &a[a_offset], lda, &vt[vt_offset],
                                ldvt, (ftnlen)1);

/*                    Copy L to WORK(IR), zeroing out above it */

                        dlacpy_((char *)"L", m, m, &a[a_offset], lda, &work[ir], &
                                ldwrkr, (ftnlen)1);
                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &work[ir +
                                ldwrkr], &ldwrkr, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M*M + M + N, prefer M*M + M + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
                                work[iwork], &i__2, &ierr);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IR) */
/*                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &work[ir], &ldwrkr, &s[1], &work[ie], &
                                work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);

/*                    Generate right bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need M*M + 4*M-1, */
/*                                prefer M*M+3*M+(M-1)*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", m, m, m, &work[ir], &ldwrkr, &work[itaup]
                                , &work[iwork], &i__2, &ierr, (ftnlen)1);
                        iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of L in WORK(IR) */
/*                    (Workspace: need M*M + BDSPAC) */

                        dbdsqr_((char *)"U", m, m, &c__0, &c__0, &s[1], &work[ie], &
                                work[ir], &ldwrkr, dum, &c__1, dum, &c__1, &
                                work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IR) by */
/*                    Q in VT, storing result in A */
/*                    (Workspace: need M*M) */

                        dgemm_((char *)"N", (char *)"N", m, n, m, &c_b79, &work[ir], &ldwrkr,
                                &vt[vt_offset], ldvt, &c_b57, &a[a_offset],
                                lda, (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

                        dlacpy_((char *)"F", m, n, &a[a_offset], lda, &vt[vt_offset],
                                ldvt, (ftnlen)1);

                    } else {

/*                    Insufficient workspace for a fast algorithm */

                        itau = 1;
                        iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", m, n, &a[a_offset], lda, &vt[vt_offset],
                                ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M + N, prefer M + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
                                work[iwork], &i__2, &ierr);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;

/*                    Zero out above L in A */

                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1
                                << 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (Workspace: need 4*M, prefer 3*M + 2*M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &
                                work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in A by Q */
/*                    in VT */
/*                    (Workspace: need 3*M + N, prefer 3*M + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"P", (char *)"L", (char *)"T", m, n, m, &a[a_offset], lda, &
                                work[itaup], &vt[vt_offset], ldvt, &work[
                                iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
                                ftnlen)1);
                        iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

                        dbdsqr_((char *)"U", m, n, &c__0, &c__0, &s[1], &work[ie], &
                                vt[vt_offset], ldvt, dum, &c__1, dum, &c__1, &
                                work[iwork], info, (ftnlen)1);

                    }

                } else if (wntuo) {

/*                 Path 8t(N much larger than M, JOBU='O', JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 M left singular vectors to be overwritten on A */

/* Computing MAX */
                    i__2 = *n + *m, i__3 = *m << 2, i__2 = max(i__2,i__3);
                    if (*lwork >= (*m << 1) * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

                        iu = 1;
                        if (*lwork >= wrkbl + (*lda << 1) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M */

                            ldwrku = *lda;
                            ir = iu + ldwrku * *m;
                            ldwrkr = *lda;
                        } else if (*lwork >= wrkbl + (*lda + *m) * *m) {

/*                       WORK(IU) is LDA by M and WORK(IR) is M by M */

                            ldwrku = *lda;
                            ir = iu + ldwrku * *m;
                            ldwrkr = *m;
                        } else {

/*                       WORK(IU) is M by M and WORK(IR) is M by M */

                            ldwrku = *m;
                            ir = iu + ldwrku * *m;
                            ldwrkr = *m;
                        }
                        itau = ir + ldwrkr * *m;
                        iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M*M + 2*M, prefer 2*M*M + M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", m, n, &a[a_offset], lda, &vt[vt_offset],
                                ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need 2*M*M + M + N, prefer 2*M*M + M + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
                                work[iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

                        dlacpy_((char *)"L", m, m, &a[a_offset], lda, &work[iu], &
                                ldwrku, (ftnlen)1);
                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &work[iu +
                                ldwrku], &ldwrku, (ftnlen)1);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to */
/*                    WORK(IR) */
/*                    (Workspace: need 2*M*M + 4*M, */
/*                                prefer 2*M*M+3*M+2*M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &
                                work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);
                        dlacpy_((char *)"L", m, m, &work[iu], &ldwrku, &work[ir], &
                                ldwrkr, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need 2*M*M + 4*M-1, */
/*                                prefer 2*M*M+3*M+(M-1)*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", m, m, m, &work[iu], &ldwrku, &work[itaup]
                                , &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in WORK(IR) */
/*                    (Workspace: need 2*M*M + 4*M, prefer 2*M*M + 3*M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", m, m, m, &work[ir], &ldwrkr, &work[itauq]
                                , &work[iwork], &i__2, &ierr, (ftnlen)1);
                        iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in WORK(IR) and computing */
/*                    right singular vectors of L in WORK(IU) */
/*                    (Workspace: need 2*M*M + BDSPAC) */

                        dbdsqr_((char *)"U", m, m, m, &c__0, &s[1], &work[ie], &work[
                                iu], &ldwrku, &work[ir], &ldwrkr, dum, &c__1,
                                &work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in VT, storing result in A */
/*                    (Workspace: need M*M) */

                        dgemm_((char *)"N", (char *)"N", m, n, m, &c_b79, &work[iu], &ldwrku,
                                &vt[vt_offset], ldvt, &c_b57, &a[a_offset],
                                lda, (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

                        dlacpy_((char *)"F", m, n, &a[a_offset], lda, &vt[vt_offset],
                                ldvt, (ftnlen)1);

/*                    Copy left singular vectors of A from WORK(IR) to A */

                        dlacpy_((char *)"F", m, m, &work[ir], &ldwrkr, &a[a_offset],
                                lda, (ftnlen)1);

                    } else {

/*                    Insufficient workspace for a fast algorithm */

                        itau = 1;
                        iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", m, n, &a[a_offset], lda, &vt[vt_offset],
                                ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M + N, prefer M + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
                                work[iwork], &i__2, &ierr);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;

/*                    Zero out above L in A */

                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &a[(a_dim1
                                << 1) + 1], lda, (ftnlen)1);

/*                    Bidiagonalize L in A */
/*                    (Workspace: need 4*M, prefer 3*M + 2*M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &a[a_offset], lda, &s[1], &work[ie], &
                                work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in A by Q */
/*                    in VT */
/*                    (Workspace: need 3*M + N, prefer 3*M + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"P", (char *)"L", (char *)"T", m, n, m, &a[a_offset], lda, &
                                work[itaup], &vt[vt_offset], ldvt, &work[
                                iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
                                ftnlen)1);

/*                    Generate left bidiagonalizing vectors in A */
/*                    (Workspace: need 4*M, prefer 3*M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", m, m, m, &a[a_offset], lda, &work[itauq],
                                 &work[iwork], &i__2, &ierr, (ftnlen)1);
                        iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in A and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

                        dbdsqr_((char *)"U", m, n, m, &c__0, &s[1], &work[ie], &vt[
                                vt_offset], ldvt, &a[a_offset], lda, dum, &
                                c__1, &work[iwork], info, (ftnlen)1);

                    }

                } else if (wntuas) {

/*                 Path 9t(N much larger than M, JOBU='S' or 'A', */
/*                         JOBVT='A') */
/*                 N right singular vectors to be computed in VT and */
/*                 M left singular vectors to be computed in U */

/* Computing MAX */
                    i__2 = *n + *m, i__3 = *m << 2, i__2 = max(i__2,i__3);
                    if (*lwork >= *m * *m + max(i__2,bdspac)) {

/*                    Sufficient workspace for a fast algorithm */

                        iu = 1;
                        if (*lwork >= wrkbl + *lda * *m) {

/*                       WORK(IU) is LDA by M */

                            ldwrku = *lda;
                        } else {

/*                       WORK(IU) is M by M */

                            ldwrku = *m;
                        }
                        itau = iu + ldwrku * *m;
                        iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", m, n, &a[a_offset], lda, &vt[vt_offset],
                                ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M*M + M + N, prefer M*M + M + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
                                work[iwork], &i__2, &ierr);

/*                    Copy L to WORK(IU), zeroing out above it */

                        dlacpy_((char *)"L", m, m, &a[a_offset], lda, &work[iu], &
                                ldwrku, (ftnlen)1);
                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &work[iu +
                                ldwrku], &ldwrku, (ftnlen)1);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;

/*                    Bidiagonalize L in WORK(IU), copying result to U */
/*                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &work[iu], &ldwrku, &s[1], &work[ie], &
                                work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);
                        dlacpy_((char *)"L", m, m, &work[iu], &ldwrku, &u[u_offset],
                                ldu, (ftnlen)1);

/*                    Generate right bidiagonalizing vectors in WORK(IU) */
/*                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + (M-1)*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"P", m, m, m, &work[iu], &ldwrku, &work[itaup]
                                , &work[iwork], &i__2, &ierr, (ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", m, m, m, &u[u_offset], ldu, &work[itauq],
                                 &work[iwork], &i__2, &ierr, (ftnlen)1);
                        iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of L in U and computing right */
/*                    singular vectors of L in WORK(IU) */
/*                    (Workspace: need M*M + BDSPAC) */

                        dbdsqr_((char *)"U", m, m, m, &c__0, &s[1], &work[ie], &work[
                                iu], &ldwrku, &u[u_offset], ldu, dum, &c__1, &
                                work[iwork], info, (ftnlen)1);

/*                    Multiply right singular vectors of L in WORK(IU) by */
/*                    Q in VT, storing result in A */
/*                    (Workspace: need M*M) */

                        dgemm_((char *)"N", (char *)"N", m, n, m, &c_b79, &work[iu], &ldwrku,
                                &vt[vt_offset], ldvt, &c_b57, &a[a_offset],
                                lda, (ftnlen)1, (ftnlen)1);

/*                    Copy right singular vectors of A from A to VT */

                        dlacpy_((char *)"F", m, n, &a[a_offset], lda, &vt[vt_offset],
                                ldvt, (ftnlen)1);

                    } else {

/*                    Insufficient workspace for a fast algorithm */

                        itau = 1;
                        iwork = itau + *m;

/*                    Compute A=L*Q, copying result to VT */
/*                    (Workspace: need 2*M, prefer M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[
                                iwork], &i__2, &ierr);
                        dlacpy_((char *)"U", m, n, &a[a_offset], lda, &vt[vt_offset],
                                ldvt, (ftnlen)1);

/*                    Generate Q in VT */
/*                    (Workspace: need M + N, prefer M + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorglq_(n, n, m, &vt[vt_offset], ldvt, &work[itau], &
                                work[iwork], &i__2, &ierr);

/*                    Copy L to U, zeroing out above it */

                        dlacpy_((char *)"L", m, m, &a[a_offset], lda, &u[u_offset],
                                ldu, (ftnlen)1);
                        i__2 = *m - 1;
                        i__3 = *m - 1;
                        dlaset_((char *)"U", &i__2, &i__3, &c_b57, &c_b57, &u[(u_dim1
                                << 1) + 1], ldu, (ftnlen)1);
                        ie = itau;
                        itauq = ie + *m;
                        itaup = itauq + *m;
                        iwork = itaup + *m;

/*                    Bidiagonalize L in U */
/*                    (Workspace: need 4*M, prefer 3*M + 2*M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dgebrd_(m, m, &u[u_offset], ldu, &s[1], &work[ie], &
                                work[itauq], &work[itaup], &work[iwork], &
                                i__2, &ierr);

/*                    Multiply right bidiagonalizing vectors in U by Q */
/*                    in VT */
/*                    (Workspace: need 3*M + N, prefer 3*M + N*NB) */

                        i__2 = *lwork - iwork + 1;
                        dormbr_((char *)"P", (char *)"L", (char *)"T", m, n, m, &u[u_offset], ldu, &
                                work[itaup], &vt[vt_offset], ldvt, &work[
                                iwork], &i__2, &ierr, (ftnlen)1, (ftnlen)1, (
                                ftnlen)1);

/*                    Generate left bidiagonalizing vectors in U */
/*                    (Workspace: need 4*M, prefer 3*M + M*NB) */

                        i__2 = *lwork - iwork + 1;
                        dorgbr_((char *)"Q", m, m, m, &u[u_offset], ldu, &work[itauq],
                                 &work[iwork], &i__2, &ierr, (ftnlen)1);
                        iwork = ie + *m;

/*                    Perform bidiagonal QR iteration, computing left */
/*                    singular vectors of A in U and computing right */
/*                    singular vectors of A in VT */
/*                    (Workspace: need BDSPAC) */

                        dbdsqr_((char *)"U", m, n, m, &c__0, &s[1], &work[ie], &vt[
                                vt_offset], ldvt, &u[u_offset], ldu, dum, &
                                c__1, &work[iwork], info, (ftnlen)1);

                    }

                }

            }

        } else {

/*           N .LT. MNTHR */

/*           Path 10t(N greater than M, but not much larger) */
/*           Reduce to bidiagonal form without LQ decomposition */

            ie = 1;
            itauq = ie + *m;
            itaup = itauq + *m;
            iwork = itaup + *m;

/*           Bidiagonalize A */
/*           (Workspace: need 3*M + N, prefer 3*M + (M + N)*NB) */

            i__2 = *lwork - iwork + 1;
            dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
                    work[itaup], &work[iwork], &i__2, &ierr);
            if (wntuas) {

/*              If left singular vectors desired in U, copy result to U */
/*              and generate left bidiagonalizing vectors in U */
/*              (Workspace: need 4*M-1, prefer 3*M + (M-1)*NB) */

                dlacpy_((char *)"L", m, m, &a[a_offset], lda, &u[u_offset], ldu, (
                        ftnlen)1);
                i__2 = *lwork - iwork + 1;
                dorgbr_((char *)"Q", m, m, n, &u[u_offset], ldu, &work[itauq], &work[
                        iwork], &i__2, &ierr, (ftnlen)1);
            }
            if (wntvas) {

/*              If right singular vectors desired in VT, copy result to */
/*              VT and generate right bidiagonalizing vectors in VT */
/*              (Workspace: need 3*M + NRVT, prefer 3*M + NRVT*NB) */

                dlacpy_((char *)"U", m, n, &a[a_offset], lda, &vt[vt_offset], ldvt, (
                        ftnlen)1);
                if (wntva) {
                    nrvt = *n;
                }
                if (wntvs) {
                    nrvt = *m;
                }
                i__2 = *lwork - iwork + 1;
                dorgbr_((char *)"P", &nrvt, n, m, &vt[vt_offset], ldvt, &work[itaup],
                        &work[iwork], &i__2, &ierr, (ftnlen)1);
            }
            if (wntuo) {

/*              If left singular vectors desired in A, generate left */
/*              bidiagonalizing vectors in A */
/*              (Workspace: need 4*M-1, prefer 3*M + (M-1)*NB) */

                i__2 = *lwork - iwork + 1;
                dorgbr_((char *)"Q", m, m, n, &a[a_offset], lda, &work[itauq], &work[
                        iwork], &i__2, &ierr, (ftnlen)1);
            }
            if (wntvo) {

/*              If right singular vectors desired in A, generate right */
/*              bidiagonalizing vectors in A */
/*              (Workspace: need 4*M, prefer 3*M + M*NB) */

                i__2 = *lwork - iwork + 1;
                dorgbr_((char *)"P", m, n, m, &a[a_offset], lda, &work[itaup], &work[
                        iwork], &i__2, &ierr, (ftnlen)1);
            }
            iwork = ie + *m;
            if (wntuas || wntuo) {
                nru = *m;
            }
            if (wntun) {
                nru = 0;
            }
            if (wntvas || wntvo) {
                ncvt = *n;
            }
            if (wntvn) {
                ncvt = 0;
            }
            if (! wntuo && ! wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in VT */
/*              (Workspace: need BDSPAC) */

                dbdsqr_((char *)"L", m, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[
                        vt_offset], ldvt, &u[u_offset], ldu, dum, &c__1, &
                        work[iwork], info, (ftnlen)1);
            } else if (! wntuo && wntvo) {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in U and computing right singular */
/*              vectors in A */
/*              (Workspace: need BDSPAC) */

                dbdsqr_((char *)"L", m, &ncvt, &nru, &c__0, &s[1], &work[ie], &a[
                        a_offset], lda, &u[u_offset], ldu, dum, &c__1, &work[
                        iwork], info, (ftnlen)1);
            } else {

/*              Perform bidiagonal QR iteration, if desired, computing */
/*              left singular vectors in A and computing right singular */
/*              vectors in VT */
/*              (Workspace: need BDSPAC) */

                dbdsqr_((char *)"L", m, &ncvt, &nru, &c__0, &s[1], &work[ie], &vt[
                        vt_offset], ldvt, &a[a_offset], lda, dum, &c__1, &
                        work[iwork], info, (ftnlen)1);
            }

        }

    }

/*     If DBDSQR failed to converge, copy unconverged superdiagonals */
/*     to WORK( 2:MINMN ) */

    if (*info != 0) {
        if (ie > 2) {
            i__2 = minmn - 1;
            for (i__ = 1; i__ <= i__2; ++i__) {
                work[i__ + 1] = work[i__ + ie - 1];
/* L50: */
            }
        }
        if (ie < 2) {
            for (i__ = minmn - 1; i__ >= 1; --i__) {
                work[i__ + 1] = work[i__ + ie - 1];
/* L60: */
            }
        }
    }

/*     Undo scaling if necessary */

    if (iscl == 1) {
        if (anrm > bignum) {
            dlascl_((char *)"G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
                    minmn, &ierr, (ftnlen)1);
        }
        if (*info != 0 && anrm > bignum) {
            i__2 = minmn - 1;
            dlascl_((char *)"G", &c__0, &c__0, &bignum, &anrm, &i__2, &c__1, &work[2],
                     &minmn, &ierr, (ftnlen)1);
        }
        if (anrm < smlnum) {
            dlascl_((char *)"G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
                    minmn, &ierr, (ftnlen)1);
        }
        if (*info != 0 && anrm < smlnum) {
            i__2 = minmn - 1;
            dlascl_((char *)"G", &c__0, &c__0, &smlnum, &anrm, &i__2, &c__1, &work[2],
                     &minmn, &ierr, (ftnlen)1);
        }
    }

/*     Return optimal workspace in WORK(1) */

    work[1] = (doublereal) maxwrk;

    return 0;

/*     End of DGESVD */

} /* dgesvd_ */

#ifdef __cplusplus
        }
#endif
