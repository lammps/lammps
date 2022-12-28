/* fortran/zlaed7.f -- translated by f2c (version 20200916).
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

static integer c__2 = 2;
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b ZLAED7 used by ZSTEDC. Computes the updated eigensystem of a diagonal matrix after modification
 by a rank-one symmetric matrix. Used when the original matrix is dense. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAED7 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaed7.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaed7.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaed7.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAED7( N, CUTPNT, QSIZ, TLVLS, CURLVL, CURPBM, D, Q, */
/*                          LDQ, RHO, INDXQ, QSTORE, QPTR, PRMPTR, PERM, */
/*                          GIVPTR, GIVCOL, GIVNUM, WORK, RWORK, IWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            CURLVL, CURPBM, CUTPNT, INFO, LDQ, N, QSIZ, */
/*      $                   TLVLS */
/*       DOUBLE PRECISION   RHO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            GIVCOL( 2, * ), GIVPTR( * ), INDXQ( * ), */
/*      $                   IWORK( * ), PERM( * ), PRMPTR( * ), QPTR( * ) */
/*       DOUBLE PRECISION   D( * ), GIVNUM( 2, * ), QSTORE( * ), RWORK( * ) */
/*       COMPLEX*16         Q( LDQ, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAED7 computes the updated eigensystem of a diagonal */
/* > matrix after modification by a rank-one symmetric matrix. This */
/* > routine is used only for the eigenproblem which requires all */
/* > eigenvalues and optionally eigenvectors of a dense or banded */
/* > Hermitian matrix that has been reduced to tridiagonal form. */
/* > */
/* >   T = Q(in) ( D(in) + RHO * Z*Z**H ) Q**H(in) = Q(out) * D(out) * Q**H(out) */
/* > */
/* >   where Z = Q**Hu, u is a vector of length N with ones in the */
/* >   CUTPNT and CUTPNT + 1 th elements and zeros elsewhere. */
/* > */
/* >    The eigenvectors of the original matrix are stored in Q, and the */
/* >    eigenvalues are in D.  The algorithm consists of three stages: */
/* > */
/* >       The first stage consists of deflating the size of the problem */
/* >       when there are multiple eigenvalues or if there is a zero in */
/* >       the Z vector.  For each such occurrence the dimension of the */
/* >       secular equation problem is reduced by one.  This stage is */
/* >       performed by the routine DLAED2. */
/* > */
/* >       The second stage consists of calculating the updated */
/* >       eigenvalues. This is done by finding the roots of the secular */
/* >       equation via the routine DLAED4 (as called by SLAED3). */
/* >       This routine also calculates the eigenvectors of the current */
/* >       problem. */
/* > */
/* >       The final stage consists of computing the updated eigenvectors */
/* >       directly using the updated eigenvalues.  The eigenvectors for */
/* >       the current problem are multiplied with the eigenvectors from */
/* >       the overall problem. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         The dimension of the symmetric tridiagonal matrix.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] CUTPNT */
/* > \verbatim */
/* >          CUTPNT is INTEGER */
/* >         Contains the location of the last eigenvalue in the leading */
/* >         sub-matrix.  min(1,N) <= CUTPNT <= N. */
/* > \endverbatim */
/* > */
/* > \param[in] QSIZ */
/* > \verbatim */
/* >          QSIZ is INTEGER */
/* >         The dimension of the unitary matrix used to reduce */
/* >         the full matrix to tridiagonal form.  QSIZ >= N. */
/* > \endverbatim */
/* > */
/* > \param[in] TLVLS */
/* > \verbatim */
/* >          TLVLS is INTEGER */
/* >         The total number of merging levels in the overall divide and */
/* >         conquer tree. */
/* > \endverbatim */
/* > */
/* > \param[in] CURLVL */
/* > \verbatim */
/* >          CURLVL is INTEGER */
/* >         The current level in the overall merge routine, */
/* >         0 <= curlvl <= tlvls. */
/* > \endverbatim */
/* > */
/* > \param[in] CURPBM */
/* > \verbatim */
/* >          CURPBM is INTEGER */
/* >         The current problem in the current level in the overall */
/* >         merge routine (counting from upper left to lower right). */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >         On entry, the eigenvalues of the rank-1-perturbed matrix. */
/* >         On exit, the eigenvalues of the repaired matrix. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* >          Q is COMPLEX*16 array, dimension (LDQ,N) */
/* >         On entry, the eigenvectors of the rank-1-perturbed matrix. */
/* >         On exit, the eigenvectors of the repaired tridiagonal matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >         The leading dimension of the array Q.  LDQ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] RHO */
/* > \verbatim */
/* >          RHO is DOUBLE PRECISION */
/* >         Contains the subdiagonal element used to create the rank-1 */
/* >         modification. */
/* > \endverbatim */
/* > */
/* > \param[out] INDXQ */
/* > \verbatim */
/* >          INDXQ is INTEGER array, dimension (N) */
/* >         This contains the permutation which will reintegrate the */
/* >         subproblem just solved back into sorted order, */
/* >         ie. D( INDXQ( I = 1, N ) ) will be in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (4*N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, */
/* >                                 dimension (3*N+2*QSIZ*N) */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (QSIZ*N) */
/* > \endverbatim */
/* > */
/* > \param[in,out] QSTORE */
/* > \verbatim */
/* >          QSTORE is DOUBLE PRECISION array, dimension (N**2+1) */
/* >         Stores eigenvectors of submatrices encountered during */
/* >         divide and conquer, packed together. QPTR points to */
/* >         beginning of the submatrices. */
/* > \endverbatim */
/* > */
/* > \param[in,out] QPTR */
/* > \verbatim */
/* >          QPTR is INTEGER array, dimension (N+2) */
/* >         List of indices pointing to beginning of submatrices stored */
/* >         in QSTORE. The submatrices are numbered starting at the */
/* >         bottom left of the divide and conquer tree, from left to */
/* >         right and bottom to top. */
/* > \endverbatim */
/* > */
/* > \param[in] PRMPTR */
/* > \verbatim */
/* >          PRMPTR is INTEGER array, dimension (N lg N) */
/* >         Contains a list of pointers which indicate where in PERM a */
/* >         level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i) */
/* >         indicates the size of the permutation and also the size of */
/* >         the full, non-deflated problem. */
/* > \endverbatim */
/* > */
/* > \param[in] PERM */
/* > \verbatim */
/* >          PERM is INTEGER array, dimension (N lg N) */
/* >         Contains the permutations (from deflation and sorting) to be */
/* >         applied to each eigenblock. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVPTR */
/* > \verbatim */
/* >          GIVPTR is INTEGER array, dimension (N lg N) */
/* >         Contains a list of pointers which indicate where in GIVCOL a */
/* >         level's Givens rotations are stored.  GIVPTR(i+1) - GIVPTR(i) */
/* >         indicates the number of Givens rotations. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVCOL */
/* > \verbatim */
/* >          GIVCOL is INTEGER array, dimension (2, N lg N) */
/* >         Each pair of numbers indicates a pair of columns to take place */
/* >         in a Givens rotation. */
/* > \endverbatim */
/* > */
/* > \param[in] GIVNUM */
/* > \verbatim */
/* >          GIVNUM is DOUBLE PRECISION array, dimension (2, N lg N) */
/* >         Each number indicates the S value to be used in the */
/* >         corresponding Givens rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  if INFO = 1, an eigenvalue did not converge */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zlaed7_(integer *n, integer *cutpnt, integer *qsiz,
        integer *tlvls, integer *curlvl, integer *curpbm, doublereal *d__,
        doublecomplex *q, integer *ldq, doublereal *rho, integer *indxq,
        doublereal *qstore, integer *qptr, integer *prmptr, integer *perm,
        integer *givptr, integer *givcol, doublereal *givnum, doublecomplex *
        work, doublereal *rwork, integer *iwork, integer *info)
{
    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2;

    /* Builtin functions */
    integer pow_ii(integer *, integer *);

    /* Local variables */
    integer i__, k, n1, n2, iq, iw, iz, ptr, indx, curr, indxc, indxp;
    extern /* Subroutine */ int dlaed9_(integer *, integer *, integer *,
            integer *, doublereal *, doublereal *, integer *, doublereal *,
            doublereal *, doublereal *, doublereal *, integer *, integer *),
            zlaed8_(integer *, integer *, integer *, doublecomplex *, integer
            *, doublereal *, doublereal *, integer *, doublereal *,
            doublereal *, doublecomplex *, integer *, doublereal *, integer *,
             integer *, integer *, integer *, integer *, integer *,
            doublereal *, integer *), dlaeda_(integer *, integer *, integer *,
             integer *, integer *, integer *, integer *, integer *,
            doublereal *, doublereal *, integer *, doublereal *, doublereal *,
             integer *);
    integer idlmda;
    extern /* Subroutine */ int dlamrg_(integer *, integer *, doublereal *,
            integer *, integer *, integer *), xerbla_(char *, integer *,
            ftnlen), zlacrm_(integer *, integer *, doublecomplex *, integer *,
             doublereal *, integer *, doublecomplex *, integer *, doublereal *
            );
    integer coltyp;


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    --d__;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --indxq;
    --qstore;
    --qptr;
    --prmptr;
    --perm;
    --givptr;
    givcol -= 3;
    givnum -= 3;
    --work;
    --rwork;
    --iwork;

    /* Function Body */
    *info = 0;

/*     IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.1 ) THEN */
/*        INFO = -1 */
/*     ELSE IF( N.LT.0 ) THEN */
    if (*n < 0) {
        *info = -1;
    } else if (min(1,*n) > *cutpnt || *n < *cutpnt) {
        *info = -2;
    } else if (*qsiz < *n) {
        *info = -3;
    } else if (*ldq < max(1,*n)) {
        *info = -9;
    }
    if (*info != 0) {
        i__1 = -(*info);
        xerbla_((char *)"ZLAED7", &i__1, (ftnlen)6);
        return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
        return 0;
    }

/*     The following values are for bookkeeping purposes only.  They are */
/*     integer pointers which indicate the portion of the workspace */
/*     used by a particular array in DLAED2 and SLAED3. */

    iz = 1;
    idlmda = iz + *n;
    iw = idlmda + *n;
    iq = iw + *n;

    indx = 1;
    indxc = indx + *n;
    coltyp = indxc + *n;
    indxp = coltyp + *n;

/*     Form the z-vector which consists of the last row of Q_1 and the */
/*     first row of Q_2. */

    ptr = pow_ii(&c__2, tlvls) + 1;
    i__1 = *curlvl - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
        i__2 = *tlvls - i__;
        ptr += pow_ii(&c__2, &i__2);
/* L10: */
    }
    curr = ptr + *curpbm;
    dlaeda_(n, tlvls, curlvl, curpbm, &prmptr[1], &perm[1], &givptr[1], &
            givcol[3], &givnum[3], &qstore[1], &qptr[1], &rwork[iz], &rwork[
            iz + *n], info);

/*     When solving the final problem, we no longer need the stored data, */
/*     so we will overwrite the data from this level onto the previously */
/*     used storage space. */

    if (*curlvl == *tlvls) {
        qptr[curr] = 1;
        prmptr[curr] = 1;
        givptr[curr] = 1;
    }

/*     Sort and Deflate eigenvalues. */

    zlaed8_(&k, n, qsiz, &q[q_offset], ldq, &d__[1], rho, cutpnt, &rwork[iz],
            &rwork[idlmda], &work[1], qsiz, &rwork[iw], &iwork[indxp], &iwork[
            indx], &indxq[1], &perm[prmptr[curr]], &givptr[curr + 1], &givcol[
            (givptr[curr] << 1) + 1], &givnum[(givptr[curr] << 1) + 1], info);
    prmptr[curr + 1] = prmptr[curr] + *n;
    givptr[curr + 1] += givptr[curr];

/*     Solve Secular Equation. */

    if (k != 0) {
        dlaed9_(&k, &c__1, &k, n, &d__[1], &rwork[iq], &k, rho, &rwork[idlmda]
                , &rwork[iw], &qstore[qptr[curr]], &k, info);
        zlacrm_(qsiz, &k, &work[1], qsiz, &qstore[qptr[curr]], &k, &q[
                q_offset], ldq, &rwork[iq]);
/* Computing 2nd power */
        i__1 = k;
        qptr[curr + 1] = qptr[curr] + i__1 * i__1;
        if (*info != 0) {
            return 0;
        }

/*     Prepare the INDXQ sorting premutation. */

        n1 = k;
        n2 = *n - k;
        dlamrg_(&n1, &n2, &d__[1], &c__1, &c_n1, &indxq[1]);
    } else {
        qptr[curr + 1] = qptr[curr];
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            indxq[i__] = i__;
/* L20: */
        }
    }

    return 0;

/*     End of ZLAED7 */

} /* zlaed7_ */

#ifdef __cplusplus
        }
#endif
