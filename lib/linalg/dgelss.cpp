/* fortran/dgelss.f -- translated by f2c (version 20200916).
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
static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b46 = 0.;
static integer c__1 = 1;
static doublereal c_b79 = 1.;

/* > \brief <b> DGELSS solves overdetermined or underdetermined systems for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGELSS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgelss.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgelss.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgelss.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, */
/*                          WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK */
/*       DOUBLE PRECISION   RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGELSS computes the minimum norm solution to a real linear least */
/* > squares problem: */
/* > */
/* > Minimize 2-norm(| b - A*x |). */
/* > */
/* > using the singular value decomposition (SVD) of A. A is an M-by-N */
/* > matrix which may be rank-deficient. */
/* > */
/* > Several right hand side vectors b and solution vectors x can be */
/* > handled in a single call; they are stored as the columns of the */
/* > M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix */
/* > X. */
/* > */
/* > The effective rank of A is determined by treating as zero those */
/* > singular values which are less than RCOND times the largest singular */
/* > value. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrices B and X. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, the first min(m,n) rows of A are overwritten with */
/* >          its right singular vectors, stored rowwise. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
/* >          On entry, the M-by-NRHS right hand side matrix B. */
/* >          On exit, B is overwritten by the N-by-NRHS solution */
/* >          matrix X.  If m >= n and RANK = n, the residual */
/* >          sum-of-squares for the solution in the i-th column is given */
/* >          by the sum of squares of elements n+1:m in that column. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. LDB >= max(1,max(M,N)). */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is DOUBLE PRECISION array, dimension (min(M,N)) */
/* >          The singular values of A in decreasing order. */
/* >          The condition number of A in the 2-norm = S(1)/S(min(m,n)). */
/* > \endverbatim */
/* > */
/* > \param[in] RCOND */
/* > \verbatim */
/* >          RCOND is DOUBLE PRECISION */
/* >          RCOND is used to determine the effective rank of A. */
/* >          Singular values S(i) <= RCOND*S(1) are treated as zero. */
/* >          If RCOND < 0, machine precision is used instead. */
/* > \endverbatim */
/* > */
/* > \param[out] RANK */
/* > \verbatim */
/* >          RANK is INTEGER */
/* >          The effective rank of A, i.e., the number of singular values */
/* >          which are greater than RCOND*S(1). */
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
/* >          The dimension of the array WORK. LWORK >= 1, and also: */
/* >          LWORK >= 3*min(M,N) + max( 2*min(M,N), max(M,N), NRHS ) */
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
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  the algorithm for computing the SVD failed to converge; */
/* >                if INFO = i, i off-diagonal elements of an intermediate */
/* >                bidiagonal form did not converge to zero. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doubleGEsolve */

/*  ===================================================================== */
/* Subroutine */ int dgelss_(integer *m, integer *n, integer *nrhs,
        doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
        s, doublereal *rcond, integer *rank, doublereal *work, integer *lwork,
         integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Local variables */
    integer i__, bl, ie, il, mm;
    doublereal dum[1], eps, thr, anrm, bnrm;
    integer itau, lwork_dgebrd__, lwork_dgelqf__, lwork_dgeqrf__,
            lwork_dorgbr__, lwork_dormbr__, lwork_dormlq__, lwork_dormqr__;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *,
            integer *, doublereal *, doublereal *, integer *, doublereal *,
            integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    integer iascl, ibscl;
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *,
            doublereal *, doublereal *, integer *, doublereal *, integer *,
            doublereal *, doublereal *, integer *, ftnlen), drscl_(integer *,
            doublereal *, doublereal *, integer *);
    integer chunk;
    doublereal sfmin;
    integer minmn;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *,
            doublereal *, integer *);
    integer maxmn, itaup, itauq, mnthr, iwork;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *), dgebrd_(
            integer *, integer *, doublereal *, integer *, doublereal *,
            doublereal *, doublereal *, doublereal *, doublereal *, integer *,
             integer *);
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
            xerbla_(char *, integer *, ftnlen), dbdsqr_(char *, integer *,
            integer *, integer *, integer *, doublereal *, doublereal *,
            doublereal *, integer *, doublereal *, integer *, doublereal *,
            integer *, doublereal *, integer *, ftnlen), dorgbr_(char *,
            integer *, integer *, integer *, doublereal *, integer *,
            doublereal *, doublereal *, integer *, integer *, ftnlen);
    doublereal bignum;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *,
            integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dormbr_(char *, char *, char *, integer *,
            integer *, integer *, doublereal *, integer *, doublereal *,
            doublereal *, integer *, doublereal *, integer *, integer *,
            ftnlen, ftnlen, ftnlen), dormlq_(char *, char *, integer *,
            integer *, integer *, doublereal *, integer *, doublereal *,
            doublereal *, integer *, doublereal *, integer *, integer *,
            ftnlen, ftnlen);
    integer ldwork;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *,
            integer *, doublereal *, integer *, doublereal *, doublereal *,
            integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    integer minwrk, maxwrk;
    doublereal smlnum;
    logical lquery;


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
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --s;
    --work;

    /* Function Body */
    *info = 0;
    minmn = min(*m,*n);
    maxmn = max(*m,*n);
    lquery = *lwork == -1;
    if (*m < 0) {
        *info = -1;
    } else if (*n < 0) {
        *info = -2;
    } else if (*nrhs < 0) {
        *info = -3;
    } else if (*lda < max(1,*m)) {
        *info = -5;
    } else if (*ldb < max(1,maxmn)) {
        *info = -7;
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
        if (minmn > 0) {
            mm = *m;
            mnthr = ilaenv_(&c__6, (char *)"DGELSS", (char *)" ", m, n, nrhs, &c_n1, (ftnlen)
                    6, (ftnlen)1);
            if (*m >= *n && *m >= mnthr) {

/*              Path 1a - overdetermined, with many more rows than */
/*                        columns */

/*              Compute space needed for DGEQRF */
                dgeqrf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, info);
                lwork_dgeqrf__ = (integer) dum[0];
/*              Compute space needed for DORMQR */
                dormqr_((char *)"L", (char *)"T", m, nrhs, n, &a[a_offset], lda, dum, &b[
                        b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (ftnlen)
                        1);
                lwork_dormqr__ = (integer) dum[0];
                mm = *n;
/* Computing MAX */
                i__1 = maxwrk, i__2 = *n + lwork_dgeqrf__;
                maxwrk = max(i__1,i__2);
/* Computing MAX */
                i__1 = maxwrk, i__2 = *n + lwork_dormqr__;
                maxwrk = max(i__1,i__2);
            }
            if (*m >= *n) {

/*              Path 1 - overdetermined or exactly determined */

/*              Compute workspace needed for DBDSQR */

/* Computing MAX */
                i__1 = 1, i__2 = *n * 5;
                bdspac = max(i__1,i__2);
/*              Compute space needed for DGEBRD */
                dgebrd_(&mm, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum,
                        &c_n1, info);
                lwork_dgebrd__ = (integer) dum[0];
/*              Compute space needed for DORMBR */
                dormbr_((char *)"Q", (char *)"L", (char *)"T", &mm, nrhs, n, &a[a_offset], lda, dum, &
                        b[b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (
                        ftnlen)1, (ftnlen)1);
                lwork_dormbr__ = (integer) dum[0];
/*              Compute space needed for DORGBR */
                dorgbr_((char *)"P", n, n, n, &a[a_offset], lda, dum, dum, &c_n1,
                        info, (ftnlen)1);
                lwork_dorgbr__ = (integer) dum[0];
/*              Compute total workspace needed */
/* Computing MAX */
                i__1 = maxwrk, i__2 = *n * 3 + lwork_dgebrd__;
                maxwrk = max(i__1,i__2);
/* Computing MAX */
                i__1 = maxwrk, i__2 = *n * 3 + lwork_dormbr__;
                maxwrk = max(i__1,i__2);
/* Computing MAX */
                i__1 = maxwrk, i__2 = *n * 3 + lwork_dorgbr__;
                maxwrk = max(i__1,i__2);
                maxwrk = max(maxwrk,bdspac);
/* Computing MAX */
                i__1 = maxwrk, i__2 = *n * *nrhs;
                maxwrk = max(i__1,i__2);
/* Computing MAX */
                i__1 = *n * 3 + mm, i__2 = *n * 3 + *nrhs, i__1 = max(i__1,
                        i__2);
                minwrk = max(i__1,bdspac);
                maxwrk = max(minwrk,maxwrk);
            }
            if (*n > *m) {

/*              Compute workspace needed for DBDSQR */

/* Computing MAX */
                i__1 = 1, i__2 = *m * 5;
                bdspac = max(i__1,i__2);
/* Computing MAX */
                i__1 = *m * 3 + *nrhs, i__2 = *m * 3 + *n, i__1 = max(i__1,
                        i__2);
                minwrk = max(i__1,bdspac);
                if (*n >= mnthr) {

/*                 Path 2a - underdetermined, with many more columns */
/*                 than rows */

/*                 Compute space needed for DGELQF */
                    dgelqf_(m, n, &a[a_offset], lda, dum, dum, &c_n1, info);
                    lwork_dgelqf__ = (integer) dum[0];
/*                 Compute space needed for DGEBRD */
                    dgebrd_(m, m, &a[a_offset], lda, &s[1], dum, dum, dum,
                            dum, &c_n1, info);
                    lwork_dgebrd__ = (integer) dum[0];
/*                 Compute space needed for DORMBR */
                    dormbr_((char *)"Q", (char *)"L", (char *)"T", m, nrhs, n, &a[a_offset], lda, dum,
                             &b[b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (
                            ftnlen)1, (ftnlen)1);
                    lwork_dormbr__ = (integer) dum[0];
/*                 Compute space needed for DORGBR */
                    dorgbr_((char *)"P", m, m, m, &a[a_offset], lda, dum, dum, &c_n1,
                            info, (ftnlen)1);
                    lwork_dorgbr__ = (integer) dum[0];
/*                 Compute space needed for DORMLQ */
                    dormlq_((char *)"L", (char *)"T", n, nrhs, m, &a[a_offset], lda, dum, &b[
                            b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (
                            ftnlen)1);
                    lwork_dormlq__ = (integer) dum[0];
/*                 Compute total workspace needed */
                    maxwrk = *m + lwork_dgelqf__;
/* Computing MAX */
                    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) +
                            lwork_dgebrd__;
                    maxwrk = max(i__1,i__2);
/* Computing MAX */
                    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) +
                            lwork_dormbr__;
                    maxwrk = max(i__1,i__2);
/* Computing MAX */
                    i__1 = maxwrk, i__2 = *m * *m + (*m << 2) +
                            lwork_dorgbr__;
                    maxwrk = max(i__1,i__2);
/* Computing MAX */
                    i__1 = maxwrk, i__2 = *m * *m + *m + bdspac;
                    maxwrk = max(i__1,i__2);
                    if (*nrhs > 1) {
/* Computing MAX */
                        i__1 = maxwrk, i__2 = *m * *m + *m + *m * *nrhs;
                        maxwrk = max(i__1,i__2);
                    } else {
/* Computing MAX */
                        i__1 = maxwrk, i__2 = *m * *m + (*m << 1);
                        maxwrk = max(i__1,i__2);
                    }
/* Computing MAX */
                    i__1 = maxwrk, i__2 = *m + lwork_dormlq__;
                    maxwrk = max(i__1,i__2);
                } else {

/*                 Path 2 - underdetermined */

/*                 Compute space needed for DGEBRD */
                    dgebrd_(m, n, &a[a_offset], lda, &s[1], dum, dum, dum,
                            dum, &c_n1, info);
                    lwork_dgebrd__ = (integer) dum[0];
/*                 Compute space needed for DORMBR */
                    dormbr_((char *)"Q", (char *)"L", (char *)"T", m, nrhs, m, &a[a_offset], lda, dum,
                             &b[b_offset], ldb, dum, &c_n1, info, (ftnlen)1, (
                            ftnlen)1, (ftnlen)1);
                    lwork_dormbr__ = (integer) dum[0];
/*                 Compute space needed for DORGBR */
                    dorgbr_((char *)"P", m, n, m, &a[a_offset], lda, dum, dum, &c_n1,
                            info, (ftnlen)1);
                    lwork_dorgbr__ = (integer) dum[0];
                    maxwrk = *m * 3 + lwork_dgebrd__;
/* Computing MAX */
                    i__1 = maxwrk, i__2 = *m * 3 + lwork_dormbr__;
                    maxwrk = max(i__1,i__2);
/* Computing MAX */
                    i__1 = maxwrk, i__2 = *m * 3 + lwork_dorgbr__;
                    maxwrk = max(i__1,i__2);
                    maxwrk = max(maxwrk,bdspac);
/* Computing MAX */
                    i__1 = maxwrk, i__2 = *n * *nrhs;
                    maxwrk = max(i__1,i__2);
                }
            }
            maxwrk = max(minwrk,maxwrk);
        }
        work[1] = (doublereal) maxwrk;

        if (*lwork < minwrk && ! lquery) {
            *info = -12;
        }
    }

    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"DGELSS", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
        *rank = 0;
        return 0;
    }

/*     Get machine parameters */

    eps = dlamch_((char *)"P", (ftnlen)1);
    sfmin = dlamch_((char *)"S", (ftnlen)1);
    smlnum = sfmin / eps;
    bignum = 1. / smlnum;
    dlabad_(&smlnum, &bignum);

/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

    anrm = dlange_((char *)"M", m, n, &a[a_offset], lda, &work[1], (ftnlen)1);
    iascl = 0;
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

        dlascl_((char *)"G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda,
                info, (ftnlen)1);
        iascl = 1;
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

        dlascl_((char *)"G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda,
                info, (ftnlen)1);
        iascl = 2;
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

        i__1 = max(*m,*n);
        dlaset_((char *)"F", &i__1, nrhs, &c_b46, &c_b46, &b[b_offset], ldb, (ftnlen)
                1);
        dlaset_((char *)"F", &minmn, &c__1, &c_b46, &c_b46, &s[1], &minmn, (ftnlen)1);
        *rank = 0;
        goto L70;
    }

/*     Scale B if max element outside range [SMLNUM,BIGNUM] */

    bnrm = dlange_((char *)"M", m, nrhs, &b[b_offset], ldb, &work[1], (ftnlen)1);
    ibscl = 0;
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

        dlascl_((char *)"G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &b[b_offset], ldb,
                 info, (ftnlen)1);
        ibscl = 1;
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

        dlascl_((char *)"G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &b[b_offset], ldb,
                 info, (ftnlen)1);
        ibscl = 2;
    }

/*     Overdetermined case */

    if (*m >= *n) {

/*        Path 1 - overdetermined or exactly determined */

        mm = *m;
        if (*m >= mnthr) {

/*           Path 1a - overdetermined, with many more rows than columns */

            mm = *n;
            itau = 1;
            iwork = itau + *n;

/*           Compute A=Q*R */
/*           (Workspace: need 2*N, prefer N+N*NB) */

            i__1 = *lwork - iwork + 1;
            dgeqrf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__1,
                     info);

/*           Multiply B by transpose(Q) */
/*           (Workspace: need N+NRHS, prefer N+NRHS*NB) */

            i__1 = *lwork - iwork + 1;
            dormqr_((char *)"L", (char *)"T", m, nrhs, n, &a[a_offset], lda, &work[itau], &b[
                    b_offset], ldb, &work[iwork], &i__1, info, (ftnlen)1, (
                    ftnlen)1);

/*           Zero out below R */

            if (*n > 1) {
                i__1 = *n - 1;
                i__2 = *n - 1;
                dlaset_((char *)"L", &i__1, &i__2, &c_b46, &c_b46, &a[a_dim1 + 2],
                        lda, (ftnlen)1);
            }
        }

        ie = 1;
        itauq = ie + *n;
        itaup = itauq + *n;
        iwork = itaup + *n;

/*        Bidiagonalize R in A */
/*        (Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB) */

        i__1 = *lwork - iwork + 1;
        dgebrd_(&mm, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
                work[itaup], &work[iwork], &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors of R */
/*        (Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB) */

        i__1 = *lwork - iwork + 1;
        dormbr_((char *)"Q", (char *)"L", (char *)"T", &mm, nrhs, n, &a[a_offset], lda, &work[itauq],
                &b[b_offset], ldb, &work[iwork], &i__1, info, (ftnlen)1, (
                ftnlen)1, (ftnlen)1);

/*        Generate right bidiagonalizing vectors of R in A */
/*        (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */

        i__1 = *lwork - iwork + 1;
        dorgbr_((char *)"P", n, n, n, &a[a_offset], lda, &work[itaup], &work[iwork], &
                i__1, info, (ftnlen)1);
        iwork = ie + *n;

/*        Perform bidiagonal QR iteration */
/*          multiply B by transpose of left singular vectors */
/*          compute right singular vectors in A */
/*        (Workspace: need BDSPAC) */

        dbdsqr_((char *)"U", n, n, &c__0, nrhs, &s[1], &work[ie], &a[a_offset], lda,
                dum, &c__1, &b[b_offset], ldb, &work[iwork], info, (ftnlen)1);
        if (*info != 0) {
            goto L70;
        }

/*        Multiply B by reciprocals of singular values */

/* Computing MAX */
        d__1 = *rcond * s[1];
        thr = max(d__1,sfmin);
        if (*rcond < 0.) {
/* Computing MAX */
            d__1 = eps * s[1];
            thr = max(d__1,sfmin);
        }
        *rank = 0;
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            if (s[i__] > thr) {
                drscl_(nrhs, &s[i__], &b[i__ + b_dim1], ldb);
                ++(*rank);
            } else {
                dlaset_((char *)"F", &c__1, nrhs, &c_b46, &c_b46, &b[i__ + b_dim1],
                        ldb, (ftnlen)1);
            }
/* L10: */
        }

/*        Multiply B by right singular vectors */
/*        (Workspace: need N, prefer N*NRHS) */

        if (*lwork >= *ldb * *nrhs && *nrhs > 1) {
            dgemm_((char *)"T", (char *)"N", n, nrhs, n, &c_b79, &a[a_offset], lda, &b[
                    b_offset], ldb, &c_b46, &work[1], ldb, (ftnlen)1, (ftnlen)
                    1);
            dlacpy_((char *)"G", n, nrhs, &work[1], ldb, &b[b_offset], ldb, (ftnlen)1)
                    ;
        } else if (*nrhs > 1) {
            chunk = *lwork / *n;
            i__1 = *nrhs;
            i__2 = chunk;
            for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
                i__3 = *nrhs - i__ + 1;
                bl = min(i__3,chunk);
                dgemm_((char *)"T", (char *)"N", n, &bl, n, &c_b79, &a[a_offset], lda, &b[i__
                        * b_dim1 + 1], ldb, &c_b46, &work[1], n, (ftnlen)1, (
                        ftnlen)1);
                dlacpy_((char *)"G", n, &bl, &work[1], n, &b[i__ * b_dim1 + 1], ldb, (
                        ftnlen)1);
/* L20: */
            }
        } else {
            dgemv_((char *)"T", n, n, &c_b79, &a[a_offset], lda, &b[b_offset], &c__1,
                    &c_b46, &work[1], &c__1, (ftnlen)1);
            dcopy_(n, &work[1], &c__1, &b[b_offset], &c__1);
        }

    } else /* if(complicated condition) */ {
/* Computing MAX */
        i__2 = *m, i__1 = (*m << 1) - 4, i__2 = max(i__2,i__1), i__2 = max(
                i__2,*nrhs), i__1 = *n - *m * 3;
        if (*n >= mnthr && *lwork >= (*m << 2) + *m * *m + max(i__2,i__1)) {

/*        Path 2a - underdetermined, with many more columns than rows */
/*        and sufficient workspace for an efficient algorithm */

            ldwork = *m;
/* Computing MAX */
/* Computing MAX */
            i__3 = *m, i__4 = (*m << 1) - 4, i__3 = max(i__3,i__4), i__3 =
                    max(i__3,*nrhs), i__4 = *n - *m * 3;
            i__2 = (*m << 2) + *m * *lda + max(i__3,i__4), i__1 = *m * *lda +
                    *m + *m * *nrhs;
            if (*lwork >= max(i__2,i__1)) {
                ldwork = *lda;
            }
            itau = 1;
            iwork = *m + 1;

/*        Compute A=L*Q */
/*        (Workspace: need 2*M, prefer M+M*NB) */

            i__2 = *lwork - iwork + 1;
            dgelqf_(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2,
                     info);
            il = iwork;

/*        Copy L to WORK(IL), zeroing out above it */

            dlacpy_((char *)"L", m, m, &a[a_offset], lda, &work[il], &ldwork, (ftnlen)
                    1);
            i__2 = *m - 1;
            i__1 = *m - 1;
            dlaset_((char *)"U", &i__2, &i__1, &c_b46, &c_b46, &work[il + ldwork], &
                    ldwork, (ftnlen)1);
            ie = il + ldwork * *m;
            itauq = ie + *m;
            itaup = itauq + *m;
            iwork = itaup + *m;

/*        Bidiagonalize L in WORK(IL) */
/*        (Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB) */

            i__2 = *lwork - iwork + 1;
            dgebrd_(m, m, &work[il], &ldwork, &s[1], &work[ie], &work[itauq],
                    &work[itaup], &work[iwork], &i__2, info);

/*        Multiply B by transpose of left bidiagonalizing vectors of L */
/*        (Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB) */

            i__2 = *lwork - iwork + 1;
            dormbr_((char *)"Q", (char *)"L", (char *)"T", m, nrhs, m, &work[il], &ldwork, &work[
                    itauq], &b[b_offset], ldb, &work[iwork], &i__2, info, (
                    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*        Generate right bidiagonalizing vectors of R in WORK(IL) */
/*        (Workspace: need M*M+5*M-1, prefer M*M+4*M+(M-1)*NB) */

            i__2 = *lwork - iwork + 1;
            dorgbr_((char *)"P", m, m, m, &work[il], &ldwork, &work[itaup], &work[
                    iwork], &i__2, info, (ftnlen)1);
            iwork = ie + *m;

/*        Perform bidiagonal QR iteration, */
/*           computing right singular vectors of L in WORK(IL) and */
/*           multiplying B by transpose of left singular vectors */
/*        (Workspace: need M*M+M+BDSPAC) */

            dbdsqr_((char *)"U", m, m, &c__0, nrhs, &s[1], &work[ie], &work[il], &
                    ldwork, &a[a_offset], lda, &b[b_offset], ldb, &work[iwork]
                    , info, (ftnlen)1);
            if (*info != 0) {
                goto L70;
            }

/*        Multiply B by reciprocals of singular values */

/* Computing MAX */
            d__1 = *rcond * s[1];
            thr = max(d__1,sfmin);
            if (*rcond < 0.) {
/* Computing MAX */
                d__1 = eps * s[1];
                thr = max(d__1,sfmin);
            }
            *rank = 0;
            i__2 = *m;
            for (i__ = 1; i__ <= i__2; ++i__) {
                if (s[i__] > thr) {
                    drscl_(nrhs, &s[i__], &b[i__ + b_dim1], ldb);
                    ++(*rank);
                } else {
                    dlaset_((char *)"F", &c__1, nrhs, &c_b46, &c_b46, &b[i__ + b_dim1]
                            , ldb, (ftnlen)1);
                }
/* L30: */
            }
            iwork = ie;

/*        Multiply B by right singular vectors of L in WORK(IL) */
/*        (Workspace: need M*M+2*M, prefer M*M+M+M*NRHS) */

            if (*lwork >= *ldb * *nrhs + iwork - 1 && *nrhs > 1) {
                dgemm_((char *)"T", (char *)"N", m, nrhs, m, &c_b79, &work[il], &ldwork, &b[
                        b_offset], ldb, &c_b46, &work[iwork], ldb, (ftnlen)1,
                        (ftnlen)1);
                dlacpy_((char *)"G", m, nrhs, &work[iwork], ldb, &b[b_offset], ldb, (
                        ftnlen)1);
            } else if (*nrhs > 1) {
                chunk = (*lwork - iwork + 1) / *m;
                i__2 = *nrhs;
                i__1 = chunk;
                for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ +=
                        i__1) {
/* Computing MIN */
                    i__3 = *nrhs - i__ + 1;
                    bl = min(i__3,chunk);
                    dgemm_((char *)"T", (char *)"N", m, &bl, m, &c_b79, &work[il], &ldwork, &
                            b[i__ * b_dim1 + 1], ldb, &c_b46, &work[iwork], m,
                             (ftnlen)1, (ftnlen)1);
                    dlacpy_((char *)"G", m, &bl, &work[iwork], m, &b[i__ * b_dim1 + 1]
                            , ldb, (ftnlen)1);
/* L40: */
                }
            } else {
                dgemv_((char *)"T", m, m, &c_b79, &work[il], &ldwork, &b[b_dim1 + 1],
                        &c__1, &c_b46, &work[iwork], &c__1, (ftnlen)1);
                dcopy_(m, &work[iwork], &c__1, &b[b_dim1 + 1], &c__1);
            }

/*        Zero out below first M rows of B */

            i__1 = *n - *m;
            dlaset_((char *)"F", &i__1, nrhs, &c_b46, &c_b46, &b[*m + 1 + b_dim1],
                    ldb, (ftnlen)1);
            iwork = itau + *m;

/*        Multiply transpose(Q) by B */
/*        (Workspace: need M+NRHS, prefer M+NRHS*NB) */

            i__1 = *lwork - iwork + 1;
            dormlq_((char *)"L", (char *)"T", n, nrhs, m, &a[a_offset], lda, &work[itau], &b[
                    b_offset], ldb, &work[iwork], &i__1, info, (ftnlen)1, (
                    ftnlen)1);

        } else {

/*        Path 2 - remaining underdetermined cases */

            ie = 1;
            itauq = ie + *m;
            itaup = itauq + *m;
            iwork = itaup + *m;

/*        Bidiagonalize A */
/*        (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB) */

            i__1 = *lwork - iwork + 1;
            dgebrd_(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &
                    work[itaup], &work[iwork], &i__1, info);

/*        Multiply B by transpose of left bidiagonalizing vectors */
/*        (Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB) */

            i__1 = *lwork - iwork + 1;
            dormbr_((char *)"Q", (char *)"L", (char *)"T", m, nrhs, n, &a[a_offset], lda, &work[itauq]
                    , &b[b_offset], ldb, &work[iwork], &i__1, info, (ftnlen)1,
                     (ftnlen)1, (ftnlen)1);

/*        Generate right bidiagonalizing vectors in A */
/*        (Workspace: need 4*M, prefer 3*M+M*NB) */

            i__1 = *lwork - iwork + 1;
            dorgbr_((char *)"P", m, n, m, &a[a_offset], lda, &work[itaup], &work[
                    iwork], &i__1, info, (ftnlen)1);
            iwork = ie + *m;

/*        Perform bidiagonal QR iteration, */
/*           computing right singular vectors of A in A and */
/*           multiplying B by transpose of left singular vectors */
/*        (Workspace: need BDSPAC) */

            dbdsqr_((char *)"L", m, n, &c__0, nrhs, &s[1], &work[ie], &a[a_offset],
                    lda, dum, &c__1, &b[b_offset], ldb, &work[iwork], info, (
                    ftnlen)1);
            if (*info != 0) {
                goto L70;
            }

/*        Multiply B by reciprocals of singular values */

/* Computing MAX */
            d__1 = *rcond * s[1];
            thr = max(d__1,sfmin);
            if (*rcond < 0.) {
/* Computing MAX */
                d__1 = eps * s[1];
                thr = max(d__1,sfmin);
            }
            *rank = 0;
            i__1 = *m;
            for (i__ = 1; i__ <= i__1; ++i__) {
                if (s[i__] > thr) {
                    drscl_(nrhs, &s[i__], &b[i__ + b_dim1], ldb);
                    ++(*rank);
                } else {
                    dlaset_((char *)"F", &c__1, nrhs, &c_b46, &c_b46, &b[i__ + b_dim1]
                            , ldb, (ftnlen)1);
                }
/* L50: */
            }

/*        Multiply B by right singular vectors of A */
/*        (Workspace: need N, prefer N*NRHS) */

            if (*lwork >= *ldb * *nrhs && *nrhs > 1) {
                dgemm_((char *)"T", (char *)"N", n, nrhs, m, &c_b79, &a[a_offset], lda, &b[
                        b_offset], ldb, &c_b46, &work[1], ldb, (ftnlen)1, (
                        ftnlen)1);
                dlacpy_((char *)"F", n, nrhs, &work[1], ldb, &b[b_offset], ldb, (
                        ftnlen)1);
            } else if (*nrhs > 1) {
                chunk = *lwork / *n;
                i__1 = *nrhs;
                i__2 = chunk;
                for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ +=
                        i__2) {
/* Computing MIN */
                    i__3 = *nrhs - i__ + 1;
                    bl = min(i__3,chunk);
                    dgemm_((char *)"T", (char *)"N", n, &bl, m, &c_b79, &a[a_offset], lda, &b[
                            i__ * b_dim1 + 1], ldb, &c_b46, &work[1], n, (
                            ftnlen)1, (ftnlen)1);
                    dlacpy_((char *)"F", n, &bl, &work[1], n, &b[i__ * b_dim1 + 1],
                            ldb, (ftnlen)1);
/* L60: */
                }
            } else {
                dgemv_((char *)"T", m, n, &c_b79, &a[a_offset], lda, &b[b_offset], &
                        c__1, &c_b46, &work[1], &c__1, (ftnlen)1);
                dcopy_(n, &work[1], &c__1, &b[b_offset], &c__1);
            }
        }
    }

/*     Undo scaling */

    if (iascl == 1) {
        dlascl_((char *)"G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &b[b_offset], ldb,
                 info, (ftnlen)1);
        dlascl_((char *)"G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &
                minmn, info, (ftnlen)1);
    } else if (iascl == 2) {
        dlascl_((char *)"G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &b[b_offset], ldb,
                 info, (ftnlen)1);
        dlascl_((char *)"G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &
                minmn, info, (ftnlen)1);
    }
    if (ibscl == 1) {
        dlascl_((char *)"G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &b[b_offset], ldb,
                 info, (ftnlen)1);
    } else if (ibscl == 2) {
        dlascl_((char *)"G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &b[b_offset], ldb,
                 info, (ftnlen)1);
    }

L70:
    work[1] = (doublereal) maxwrk;
    return 0;

/*     End of DGELSS */

} /* dgelss_ */

#ifdef __cplusplus
        }
#endif
