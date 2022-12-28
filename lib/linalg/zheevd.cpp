/* fortran/zheevd.f -- translated by f2c (version 20200916).
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
static integer c__0 = 0;
static doublereal c_b18 = 1.;

/* > \brief <b> ZHEEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE mat
rices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHEEVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zheevd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zheevd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zheevd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, */
/*                          LRWORK, IWORK, LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, LDA, LIWORK, LRWORK, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   RWORK( * ), W( * ) */
/*       COMPLEX*16         A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHEEVD computes all eigenvalues and, optionally, eigenvectors of a */
/* > complex Hermitian matrix A.  If eigenvectors are desired, it uses a */
/* > divide and conquer algorithm. */
/* > */
/* > The divide and conquer algorithm makes very mild assumptions about */
/* > floating point arithmetic. It will work on machines with a guard */
/* > digit in add/subtract, or on those binary machines without guard */
/* > digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or */
/* > Cray-2. It could conceivably fail on hexadecimal or decimal machines */
/* > without guard digits, but we know of none. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBZ */
/* > \verbatim */
/* >          JOBZ is CHARACTER*1 */
/* >          = 'N':  Compute eigenvalues only; */
/* >          = 'V':  Compute eigenvalues and eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored; */
/* >          = 'L':  Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA, N) */
/* >          On entry, the Hermitian matrix A.  If UPLO = 'U', the */
/* >          leading N-by-N upper triangular part of A contains the */
/* >          upper triangular part of the matrix A.  If UPLO = 'L', */
/* >          the leading N-by-N lower triangular part of A contains */
/* >          the lower triangular part of the matrix A. */
/* >          On exit, if JOBZ = 'V', then if INFO = 0, A contains the */
/* >          orthonormal eigenvectors of the matrix A. */
/* >          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L') */
/* >          or the upper triangle (if UPLO='U') of A, including the */
/* >          diagonal, is destroyed. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >          If INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The length of the array WORK. */
/* >          If N <= 1,                LWORK must be at least 1. */
/* >          If JOBZ  = 'N' and N > 1, LWORK must be at least N + 1. */
/* >          If JOBZ  = 'V' and N > 1, LWORK must be at least 2*N + N**2. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal sizes of the WORK, RWORK and */
/* >          IWORK arrays, returns these values as the first entries of */
/* >          the WORK, RWORK and IWORK arrays, and no error message */
/* >          related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, */
/* >                                         dimension (LRWORK) */
/* >          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* >          LRWORK is INTEGER */
/* >          The dimension of the array RWORK. */
/* >          If N <= 1,                LRWORK must be at least 1. */
/* >          If JOBZ  = 'N' and N > 1, LRWORK must be at least N. */
/* >          If JOBZ  = 'V' and N > 1, LRWORK must be at least */
/* >                         1 + 5*N + 2*N**2. */
/* > */
/* >          If LRWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the optimal sizes of the WORK, RWORK */
/* >          and IWORK arrays, returns these values as the first entries */
/* >          of the WORK, RWORK and IWORK arrays, and no error message */
/* >          related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* >          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* >          LIWORK is INTEGER */
/* >          The dimension of the array IWORK. */
/* >          If N <= 1,                LIWORK must be at least 1. */
/* >          If JOBZ  = 'N' and N > 1, LIWORK must be at least 1. */
/* >          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N. */
/* > */
/* >          If LIWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the optimal sizes of the WORK, RWORK */
/* >          and IWORK arrays, returns these values as the first entries */
/* >          of the WORK, RWORK and IWORK arrays, and no error message */
/* >          related to LWORK or LRWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed */
/* >                to converge; i off-diagonal elements of an intermediate */
/* >                tridiagonal form did not converge to zero; */
/* >                if INFO = i and JOBZ = 'V', then the algorithm failed */
/* >                to compute an eigenvalue while working on the submatrix */
/* >                lying in rows and columns INFO/(N+1) through */
/* >                mod(INFO,N+1). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complex16HEeigen */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* >  Modified description of INFO. Sven, 16 Feb 05. */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Jeff Rutter, Computer Science Division, University of California */
/* > at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zheevd_(char *jobz, char *uplo, integer *n,
        doublecomplex *a, integer *lda, doublereal *w, doublecomplex *work,
        integer *lwork, doublereal *rwork, integer *lrwork, integer *iwork,
        integer *liwork, integer *info, ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal eps;
    integer inde;
    doublereal anrm;
    integer imax;
    doublereal rmin, rmax;
    integer lopt;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *,
            integer *);
    doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer iinfo, lwmin, liopt;
    logical lower;
    integer llrwk, lropt;
    logical wantz;
    integer indwk2, llwrk2;
    extern doublereal dlamch_(char *, ftnlen);
    integer iscale;
    doublereal safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *,
            integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    doublereal bignum;
    extern doublereal zlanhe_(char *, char *, integer *, doublecomplex *,
            integer *, doublereal *, ftnlen, ftnlen);
    integer indtau;
    extern /* Subroutine */ int dsterf_(integer *, doublereal *, doublereal *,
             integer *), zlascl_(char *, integer *, integer *, doublereal *,
            doublereal *, integer *, integer *, doublecomplex *, integer *,
            integer *, ftnlen), zstedc_(char *, integer *, doublereal *,
            doublereal *, doublecomplex *, integer *, doublecomplex *,
            integer *, doublereal *, integer *, integer *, integer *, integer
            *, ftnlen);
    integer indrwk, indwrk, liwmin;
    extern /* Subroutine */ int zhetrd_(char *, integer *, doublecomplex *,
            integer *, doublereal *, doublereal *, doublecomplex *,
            doublecomplex *, integer *, integer *, ftnlen), zlacpy_(char *,
            integer *, integer *, doublecomplex *, integer *, doublecomplex *,
             integer *, ftnlen);
    integer lrwmin, llwork;
    doublereal smlnum;
    logical lquery;
    extern /* Subroutine */ int zunmtr_(char *, char *, char *, integer *,
            integer *, doublecomplex *, integer *, doublecomplex *,
            doublecomplex *, integer *, doublecomplex *, integer *, integer *,
             ftnlen, ftnlen, ftnlen);


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
    --w;
    --work;
    --rwork;
    --iwork;

    /* Function Body */
    wantz = lsame_(jobz, (char *)"V", (ftnlen)1, (ftnlen)1);
    lower = lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1);
    lquery = *lwork == -1 || *lrwork == -1 || *liwork == -1;

    *info = 0;
    if (! (wantz || lsame_(jobz, (char *)"N", (ftnlen)1, (ftnlen)1))) {
        *info = -1;
    } else if (! (lower || lsame_(uplo, (char *)"U", (ftnlen)1, (ftnlen)1))) {
        *info = -2;
    } else if (*n < 0) {
        *info = -3;
    } else if (*lda < max(1,*n)) {
        *info = -5;
    }

    if (*info == 0) {
        if (*n <= 1) {
            lwmin = 1;
            lrwmin = 1;
            liwmin = 1;
            lopt = lwmin;
            lropt = lrwmin;
            liopt = liwmin;
        } else {
            if (wantz) {
                lwmin = (*n << 1) + *n * *n;
/* Computing 2nd power */
                i__1 = *n;
                lrwmin = *n * 5 + 1 + (i__1 * i__1 << 1);
                liwmin = *n * 5 + 3;
            } else {
                lwmin = *n + 1;
                lrwmin = *n;
                liwmin = 1;
            }
/* Computing MAX */
            i__1 = lwmin, i__2 = *n + *n * ilaenv_(&c__1, (char *)"ZHETRD", uplo, n, &
                    c_n1, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
            lopt = max(i__1,i__2);
            lropt = lrwmin;
            liopt = liwmin;
        }
        work[1].r = (doublereal) lopt, work[1].i = 0.;
        rwork[1] = (doublereal) lropt;
        iwork[1] = liopt;

        if (*lwork < lwmin && ! lquery) {
            *info = -8;
        } else if (*lrwork < lrwmin && ! lquery) {
            *info = -10;
        } else if (*liwork < liwmin && ! lquery) {
            *info = -12;
        }
    }

    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"ZHEEVD", &i__1, (ftnlen)6);
        return 0;
    } else if (lquery) {
        return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
        return 0;
    }

    if (*n == 1) {
        i__1 = a_dim1 + 1;
        w[1] = a[i__1].r;
        if (wantz) {
            i__1 = a_dim1 + 1;
            a[i__1].r = 1., a[i__1].i = 0.;
        }
        return 0;
    }

/*     Get machine constants. */

    safmin = dlamch_((char *)"Safe minimum", (ftnlen)12);
    eps = dlamch_((char *)"Precision", (ftnlen)9);
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = sqrt(smlnum);
    rmax = sqrt(bignum);

/*     Scale matrix to allowable range, if necessary. */

    anrm = zlanhe_((char *)"M", uplo, n, &a[a_offset], lda, &rwork[1], (ftnlen)1, (
            ftnlen)1);
    iscale = 0;
    if (anrm > 0. && anrm < rmin) {
        iscale = 1;
        sigma = rmin / anrm;
    } else if (anrm > rmax) {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if (iscale == 1) {
        zlascl_(uplo, &c__0, &c__0, &c_b18, &sigma, n, n, &a[a_offset], lda,
                info, (ftnlen)1);
    }

/*     Call ZHETRD to reduce Hermitian matrix to tridiagonal form. */

    inde = 1;
    indtau = 1;
    indwrk = indtau + *n;
    indrwk = inde + *n;
    indwk2 = indwrk + *n * *n;
    llwork = *lwork - indwrk + 1;
    llwrk2 = *lwork - indwk2 + 1;
    llrwk = *lrwork - indrwk + 1;
    zhetrd_(uplo, n, &a[a_offset], lda, &w[1], &rwork[inde], &work[indtau], &
            work[indwrk], &llwork, &iinfo, (ftnlen)1);

/*     For eigenvalues only, call DSTERF.  For eigenvectors, first call */
/*     ZSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the */
/*     tridiagonal matrix, then call ZUNMTR to multiply it to the */
/*     Householder transformations represented as Householder vectors in */
/*     A. */

    if (! wantz) {
        dsterf_(n, &w[1], &rwork[inde], info);
    } else {
        zstedc_((char *)"I", n, &w[1], &rwork[inde], &work[indwrk], n, &work[indwk2],
                &llwrk2, &rwork[indrwk], &llrwk, &iwork[1], liwork, info, (
                ftnlen)1);
        zunmtr_((char *)"L", uplo, (char *)"N", n, n, &a[a_offset], lda, &work[indtau], &work[
                indwrk], n, &work[indwk2], &llwrk2, &iinfo, (ftnlen)1, (
                ftnlen)1, (ftnlen)1);
        zlacpy_((char *)"A", n, n, &work[indwrk], n, &a[a_offset], lda, (ftnlen)1);
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

    if (iscale == 1) {
        if (*info == 0) {
            imax = *n;
        } else {
            imax = *info - 1;
        }
        d__1 = 1. / sigma;
        dscal_(&imax, &d__1, &w[1], &c__1);
    }

    work[1].r = (doublereal) lopt, work[1].i = 0.;
    rwork[1] = (doublereal) lropt;
    iwork[1] = liopt;

    return 0;

/*     End of ZHEEVD */

} /* zheevd_ */

#ifdef __cplusplus
        }
#endif
