/* fortran/zheev.f -- translated by f2c (version 20200916).
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

/* > \brief <b> ZHEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matr
ices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHEEV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zheev.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zheev.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zheev.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, */
/*                         INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBZ, UPLO */
/*       INTEGER            INFO, LDA, LWORK, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   RWORK( * ), W( * ) */
/*       COMPLEX*16         A( LDA, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHEEV computes all eigenvalues and, optionally, eigenvectors of a */
/* > complex Hermitian matrix A. */
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
/* >          The length of the array WORK.  LWORK >= max(1,2*N-1). */
/* >          For optimal efficiency, LWORK >= (NB+1)*N, */
/* >          where NB is the blocksize for ZHETRD returned by ILAENV. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (max(1, 3*N-2)) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, the algorithm failed to converge; i */
/* >                off-diagonal elements of an intermediate tridiagonal */
/* >                form did not converge to zero. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complex16HEeigen */

/*  ===================================================================== */
/* Subroutine */ int zheev_(char *jobz, char *uplo, integer *n, doublecomplex
        *a, integer *lda, doublereal *w, doublecomplex *work, integer *lwork,
        doublereal *rwork, integer *info, ftnlen jobz_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer nb;
    doublereal eps;
    integer inde;
    doublereal anrm;
    integer imax;
    doublereal rmin, rmax;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *,
            integer *);
    doublereal sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    integer iinfo;
    logical lower, wantz;
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
            integer *, ftnlen);
    integer indwrk;
    extern /* Subroutine */ int zhetrd_(char *, integer *, doublecomplex *,
            integer *, doublereal *, doublereal *, doublecomplex *,
            doublecomplex *, integer *, integer *, ftnlen);
    integer llwork;
    doublereal smlnum;
    integer lwkopt;
    logical lquery;
    extern /* Subroutine */ int zsteqr_(char *, integer *, doublereal *,
            doublereal *, doublecomplex *, integer *, doublereal *, integer *,
             ftnlen), zungtr_(char *, integer *, doublecomplex *, integer *,
            doublecomplex *, doublecomplex *, integer *, integer *, ftnlen);


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

    /* Function Body */
    wantz = lsame_(jobz, (char *)"V", (ftnlen)1, (ftnlen)1);
    lower = lsame_(uplo, (char *)"L", (ftnlen)1, (ftnlen)1);
    lquery = *lwork == -1;

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
        nb = ilaenv_(&c__1, (char *)"ZHETRD", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6,
                 (ftnlen)1);
/* Computing MAX */
        i__1 = 1, i__2 = (nb + 1) * *n;
        lwkopt = max(i__1,i__2);
        work[1].r = (doublereal) lwkopt, work[1].i = 0.;

/* Computing MAX */
        i__1 = 1, i__2 = (*n << 1) - 1;
        if (*lwork < max(i__1,i__2) && ! lquery) {
            *info = -8;
        }
    }

    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"ZHEEV ", &i__1, (ftnlen)6);
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
        work[1].r = 1., work[1].i = 0.;
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
    llwork = *lwork - indwrk + 1;
    zhetrd_(uplo, n, &a[a_offset], lda, &w[1], &rwork[inde], &work[indtau], &
            work[indwrk], &llwork, &iinfo, (ftnlen)1);

/*     For eigenvalues only, call DSTERF.  For eigenvectors, first call */
/*     ZUNGTR to generate the unitary matrix, then call ZSTEQR. */

    if (! wantz) {
        dsterf_(n, &w[1], &rwork[inde], info);
    } else {
        zungtr_(uplo, n, &a[a_offset], lda, &work[indtau], &work[indwrk], &
                llwork, &iinfo, (ftnlen)1);
        indwrk = inde + *n;
        zsteqr_(jobz, n, &w[1], &rwork[inde], &a[a_offset], lda, &rwork[
                indwrk], info, (ftnlen)1);
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

/*     Set WORK(1) to optimal complex workspace size. */

    work[1].r = (doublereal) lwkopt, work[1].i = 0.;

    return 0;

/*     End of ZHEEV */

} /* zheev_ */

#ifdef __cplusplus
        }
#endif
