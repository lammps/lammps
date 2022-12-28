/* fortran/dlarfb.f -- translated by f2c (version 20200916).
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
static doublereal c_b14 = 1.;
static doublereal c_b25 = -1.;

/* > \brief \b DLARFB applies a block reflector or its transpose to a general rectangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLARFB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarfb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarfb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, */
/*                          T, LDT, C, LDC, WORK, LDWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIRECT, SIDE, STOREV, TRANS */
/*       INTEGER            K, LDC, LDT, LDV, LDWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ), */
/*      $                   WORK( LDWORK, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARFB applies a real block reflector H or its transpose H**T to a */
/* > real m by n matrix C, from either the left or the right. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply H or H**T from the Left */
/* >          = 'R': apply H or H**T from the Right */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N': apply H (No transpose) */
/* >          = 'T': apply H**T (Transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] DIRECT */
/* > \verbatim */
/* >          DIRECT is CHARACTER*1 */
/* >          Indicates how H is formed from a product of elementary */
/* >          reflectors */
/* >          = 'F': H = H(1) H(2) . . . H(k) (Forward) */
/* >          = 'B': H = H(k) . . . H(2) H(1) (Backward) */
/* > \endverbatim */
/* > */
/* > \param[in] STOREV */
/* > \verbatim */
/* >          STOREV is CHARACTER*1 */
/* >          Indicates how the vectors which define the elementary */
/* >          reflectors are stored: */
/* >          = 'C': Columnwise */
/* >          = 'R': Rowwise */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The order of the matrix T (= the number of elementary */
/* >          reflectors whose product defines the block reflector). */
/* >          If SIDE = 'L', M >= K >= 0; */
/* >          if SIDE = 'R', N >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* >          V is DOUBLE PRECISION array, dimension */
/* >                                (LDV,K) if STOREV = 'C' */
/* >                                (LDV,M) if STOREV = 'R' and SIDE = 'L' */
/* >                                (LDV,N) if STOREV = 'R' and SIDE = 'R' */
/* >          The matrix V. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >          The leading dimension of the array V. */
/* >          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M); */
/* >          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N); */
/* >          if STOREV = 'R', LDV >= K. */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* >          T is DOUBLE PRECISION array, dimension (LDT,K) */
/* >          The triangular k by k matrix T in the representation of the */
/* >          block reflector. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T. LDT >= K. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (LDC,N) */
/* >          On entry, the m by n matrix C. */
/* >          On exit, C is overwritten by H*C or H**T*C or C*H or C*H**T. */
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
/* >          WORK is DOUBLE PRECISION array, dimension (LDWORK,K) */
/* > \endverbatim */
/* > */
/* > \param[in] LDWORK */
/* > \verbatim */
/* >          LDWORK is INTEGER */
/* >          The leading dimension of the array WORK. */
/* >          If SIDE = 'L', LDWORK >= max(1,N); */
/* >          if SIDE = 'R', LDWORK >= max(1,M). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doubleOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The shape of the matrix V and the storage of the vectors which define */
/* >  the H(i) is best illustrated by the following example with n = 5 and */
/* >  k = 3. The elements equal to 1 are not stored; the corresponding */
/* >  array elements are modified but restored on exit. The rest of the */
/* >  array is not used. */
/* > */
/* >  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R': */
/* > */
/* >               V = (  1       )                 V = (  1 v1 v1 v1 v1 ) */
/* >                   ( v1  1    )                     (     1 v2 v2 v2 ) */
/* >                   ( v1 v2  1 )                     (        1 v3 v3 ) */
/* >                   ( v1 v2 v3 ) */
/* >                   ( v1 v2 v3 ) */
/* > */
/* >  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R': */
/* > */
/* >               V = ( v1 v2 v3 )                 V = ( v1 v1  1       ) */
/* >                   ( v1 v2 v3 )                     ( v2 v2 v2  1    ) */
/* >                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 ) */
/* >                   (     1 v3 ) */
/* >                   (        1 ) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlarfb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, doublereal *v, integer *
	ldv, doublereal *t, integer *ldt, doublereal *c__, integer *ldc, 
	doublereal *work, integer *ldwork, ftnlen side_len, ftnlen trans_len, 
	ftnlen direct_len, ftnlen storev_len)
{
    /* System generated locals */
    integer c_dim1, c_offset, t_dim1, t_offset, v_dim1, v_offset, work_dim1, 
	    work_offset, i__1, i__2;

    /* Local variables */
    integer i__, j;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    char transt[1];


/*  -- LAPACK auxiliary routine -- */
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
/*     .. Executable Statements .. */

/*     Quick return if possible */

    /* Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    work_dim1 = *ldwork;
    work_offset = 1 + work_dim1;
    work -= work_offset;

    /* Function Body */
    if (*m <= 0 || *n <= 0) {
	return 0;
    }

    if (lsame_(trans, (char *)"N", (ftnlen)1, (ftnlen)1)) {
	*(unsigned char *)transt = 'T';
    } else {
	*(unsigned char *)transt = 'N';
    }

    if (lsame_(storev, (char *)"C", (ftnlen)1, (ftnlen)1)) {

	if (lsame_(direct, (char *)"F", (ftnlen)1, (ftnlen)1)) {

/*           Let  V =  ( V1 )    (first K rows) */
/*                     ( V2 ) */
/*           where  V1  is unit lower triangular. */

	    if (lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1)) {

/*              Form  H * C  or  H**T * C  where  C = ( C1 ) */
/*                                                    ( C2 ) */

/*              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK) */

/*              W := C1**T */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    dcopy_(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1],
			     &c__1);
/* L10: */
		}

/*              W := W * V1 */

		dtrmm_((char *)"Right", (char *)"Lower", (char *)"No transpose", (char *)"Unit", n, k, &c_b14,
			 &v[v_offset], ldv, &work[work_offset], ldwork, (
			ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
		if (*m > *k) {

/*                 W := W + C2**T * V2 */

		    i__1 = *m - *k;
		    dgemm_((char *)"Transpose", (char *)"No transpose", n, k, &i__1, &c_b14, &
			    c__[*k + 1 + c_dim1], ldc, &v[*k + 1 + v_dim1], 
			    ldv, &c_b14, &work[work_offset], ldwork, (ftnlen)
			    9, (ftnlen)12);
		}

/*              W := W * T**T  or  W * T */

		dtrmm_((char *)"Right", (char *)"Upper", transt, (char *)"Non-unit", n, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - V * W**T */

		if (*m > *k) {

/*                 C2 := C2 - V2 * W**T */

		    i__1 = *m - *k;
		    dgemm_((char *)"No transpose", (char *)"Transpose", &i__1, n, k, &c_b25, &
			    v[*k + 1 + v_dim1], ldv, &work[work_offset], 
			    ldwork, &c_b14, &c__[*k + 1 + c_dim1], ldc, (
			    ftnlen)12, (ftnlen)9);
		}

/*              W := W * V1**T */

		dtrmm_((char *)"Right", (char *)"Lower", (char *)"Transpose", (char *)"Unit", n, k, &c_b14, &
			v[v_offset], ldv, &work[work_offset], ldwork, (ftnlen)
			5, (ftnlen)5, (ftnlen)9, (ftnlen)4);

/*              C1 := C1 - W**T */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[j + i__ * c_dim1] -= work[i__ + j * work_dim1];
/* L20: */
		    }
/* L30: */
		}

	    } else if (lsame_(side, (char *)"R", (ftnlen)1, (ftnlen)1)) {

/*              Form  C * H  or  C * H**T  where  C = ( C1  C2 ) */

/*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK) */

/*              W := C1 */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    dcopy_(m, &c__[j * c_dim1 + 1], &c__1, &work[j * 
			    work_dim1 + 1], &c__1);
/* L40: */
		}

/*              W := W * V1 */

		dtrmm_((char *)"Right", (char *)"Lower", (char *)"No transpose", (char *)"Unit", m, k, &c_b14,
			 &v[v_offset], ldv, &work[work_offset], ldwork, (
			ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
		if (*n > *k) {

/*                 W := W + C2 * V2 */

		    i__1 = *n - *k;
		    dgemm_((char *)"No transpose", (char *)"No transpose", m, k, &i__1, &
			    c_b14, &c__[(*k + 1) * c_dim1 + 1], ldc, &v[*k + 
			    1 + v_dim1], ldv, &c_b14, &work[work_offset], 
			    ldwork, (ftnlen)12, (ftnlen)12);
		}

/*              W := W * T  or  W * T**T */

		dtrmm_((char *)"Right", (char *)"Upper", trans, (char *)"Non-unit", m, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - W * V**T */

		if (*n > *k) {

/*                 C2 := C2 - W * V2**T */

		    i__1 = *n - *k;
		    dgemm_((char *)"No transpose", (char *)"Transpose", m, &i__1, k, &c_b25, &
			    work[work_offset], ldwork, &v[*k + 1 + v_dim1], 
			    ldv, &c_b14, &c__[(*k + 1) * c_dim1 + 1], ldc, (
			    ftnlen)12, (ftnlen)9);
		}

/*              W := W * V1**T */

		dtrmm_((char *)"Right", (char *)"Lower", (char *)"Transpose", (char *)"Unit", m, k, &c_b14, &
			v[v_offset], ldv, &work[work_offset], ldwork, (ftnlen)
			5, (ftnlen)5, (ftnlen)9, (ftnlen)4);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] -= work[i__ + j * work_dim1];
/* L50: */
		    }
/* L60: */
		}
	    }

	} else {

/*           Let  V =  ( V1 ) */
/*                     ( V2 )    (last K rows) */
/*           where  V2  is unit upper triangular. */

	    if (lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1)) {

/*              Form  H * C  or  H**T * C  where  C = ( C1 ) */
/*                                                    ( C2 ) */

/*              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK) */

/*              W := C2**T */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    dcopy_(n, &c__[*m - *k + j + c_dim1], ldc, &work[j * 
			    work_dim1 + 1], &c__1);
/* L70: */
		}

/*              W := W * V2 */

		dtrmm_((char *)"Right", (char *)"Upper", (char *)"No transpose", (char *)"Unit", n, k, &c_b14,
			 &v[*m - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
		if (*m > *k) {

/*                 W := W + C1**T * V1 */

		    i__1 = *m - *k;
		    dgemm_((char *)"Transpose", (char *)"No transpose", n, k, &i__1, &c_b14, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &c_b14, &
			    work[work_offset], ldwork, (ftnlen)9, (ftnlen)12);
		}

/*              W := W * T**T  or  W * T */

		dtrmm_((char *)"Right", (char *)"Lower", transt, (char *)"Non-unit", n, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - V * W**T */

		if (*m > *k) {

/*                 C1 := C1 - V1 * W**T */

		    i__1 = *m - *k;
		    dgemm_((char *)"No transpose", (char *)"Transpose", &i__1, n, k, &c_b25, &
			    v[v_offset], ldv, &work[work_offset], ldwork, &
			    c_b14, &c__[c_offset], ldc, (ftnlen)12, (ftnlen)9)
			    ;
		}

/*              W := W * V2**T */

		dtrmm_((char *)"Right", (char *)"Upper", (char *)"Transpose", (char *)"Unit", n, k, &c_b14, &
			v[*m - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)4);

/*              C2 := C2 - W**T */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[*m - *k + j + i__ * c_dim1] -= work[i__ + j * 
				work_dim1];
/* L80: */
		    }
/* L90: */
		}

	    } else if (lsame_(side, (char *)"R", (ftnlen)1, (ftnlen)1)) {

/*              Form  C * H  or  C * H**T  where  C = ( C1  C2 ) */

/*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK) */

/*              W := C2 */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    dcopy_(m, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &work[
			    j * work_dim1 + 1], &c__1);
/* L100: */
		}

/*              W := W * V2 */

		dtrmm_((char *)"Right", (char *)"Upper", (char *)"No transpose", (char *)"Unit", m, k, &c_b14,
			 &v[*n - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
		if (*n > *k) {

/*                 W := W + C1 * V1 */

		    i__1 = *n - *k;
		    dgemm_((char *)"No transpose", (char *)"No transpose", m, k, &i__1, &
			    c_b14, &c__[c_offset], ldc, &v[v_offset], ldv, &
			    c_b14, &work[work_offset], ldwork, (ftnlen)12, (
			    ftnlen)12);
		}

/*              W := W * T  or  W * T**T */

		dtrmm_((char *)"Right", (char *)"Lower", trans, (char *)"Non-unit", m, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - W * V**T */

		if (*n > *k) {

/*                 C1 := C1 - W * V1**T */

		    i__1 = *n - *k;
		    dgemm_((char *)"No transpose", (char *)"Transpose", m, &i__1, k, &c_b25, &
			    work[work_offset], ldwork, &v[v_offset], ldv, &
			    c_b14, &c__[c_offset], ldc, (ftnlen)12, (ftnlen)9)
			    ;
		}

/*              W := W * V2**T */

		dtrmm_((char *)"Right", (char *)"Upper", (char *)"Transpose", (char *)"Unit", m, k, &c_b14, &
			v[*n - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)4);

/*              C2 := C2 - W */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + (*n - *k + j) * c_dim1] -= work[i__ + j * 
				work_dim1];
/* L110: */
		    }
/* L120: */
		}
	    }
	}

    } else if (lsame_(storev, (char *)"R", (ftnlen)1, (ftnlen)1)) {

	if (lsame_(direct, (char *)"F", (ftnlen)1, (ftnlen)1)) {

/*           Let  V =  ( V1  V2 )    (V1: first K columns) */
/*           where  V1  is unit upper triangular. */

	    if (lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1)) {

/*              Form  H * C  or  H**T * C  where  C = ( C1 ) */
/*                                                    ( C2 ) */

/*              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK) */

/*              W := C1**T */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    dcopy_(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1],
			     &c__1);
/* L130: */
		}

/*              W := W * V1**T */

		dtrmm_((char *)"Right", (char *)"Upper", (char *)"Transpose", (char *)"Unit", n, k, &c_b14, &
			v[v_offset], ldv, &work[work_offset], ldwork, (ftnlen)
			5, (ftnlen)5, (ftnlen)9, (ftnlen)4);
		if (*m > *k) {

/*                 W := W + C2**T * V2**T */

		    i__1 = *m - *k;
		    dgemm_((char *)"Transpose", (char *)"Transpose", n, k, &i__1, &c_b14, &
			    c__[*k + 1 + c_dim1], ldc, &v[(*k + 1) * v_dim1 + 
			    1], ldv, &c_b14, &work[work_offset], ldwork, (
			    ftnlen)9, (ftnlen)9);
		}

/*              W := W * T**T  or  W * T */

		dtrmm_((char *)"Right", (char *)"Upper", transt, (char *)"Non-unit", n, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - V**T * W**T */

		if (*m > *k) {

/*                 C2 := C2 - V2**T * W**T */

		    i__1 = *m - *k;
		    dgemm_((char *)"Transpose", (char *)"Transpose", &i__1, n, k, &c_b25, &v[(
			    *k + 1) * v_dim1 + 1], ldv, &work[work_offset], 
			    ldwork, &c_b14, &c__[*k + 1 + c_dim1], ldc, (
			    ftnlen)9, (ftnlen)9);
		}

/*              W := W * V1 */

		dtrmm_((char *)"Right", (char *)"Upper", (char *)"No transpose", (char *)"Unit", n, k, &c_b14,
			 &v[v_offset], ldv, &work[work_offset], ldwork, (
			ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);

/*              C1 := C1 - W**T */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[j + i__ * c_dim1] -= work[i__ + j * work_dim1];
/* L140: */
		    }
/* L150: */
		}

	    } else if (lsame_(side, (char *)"R", (ftnlen)1, (ftnlen)1)) {

/*              Form  C * H  or  C * H**T  where  C = ( C1  C2 ) */

/*              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK) */

/*              W := C1 */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    dcopy_(m, &c__[j * c_dim1 + 1], &c__1, &work[j * 
			    work_dim1 + 1], &c__1);
/* L160: */
		}

/*              W := W * V1**T */

		dtrmm_((char *)"Right", (char *)"Upper", (char *)"Transpose", (char *)"Unit", m, k, &c_b14, &
			v[v_offset], ldv, &work[work_offset], ldwork, (ftnlen)
			5, (ftnlen)5, (ftnlen)9, (ftnlen)4);
		if (*n > *k) {

/*                 W := W + C2 * V2**T */

		    i__1 = *n - *k;
		    dgemm_((char *)"No transpose", (char *)"Transpose", m, k, &i__1, &c_b14, &
			    c__[(*k + 1) * c_dim1 + 1], ldc, &v[(*k + 1) * 
			    v_dim1 + 1], ldv, &c_b14, &work[work_offset], 
			    ldwork, (ftnlen)12, (ftnlen)9);
		}

/*              W := W * T  or  W * T**T */

		dtrmm_((char *)"Right", (char *)"Upper", trans, (char *)"Non-unit", m, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - W * V */

		if (*n > *k) {

/*                 C2 := C2 - W * V2 */

		    i__1 = *n - *k;
		    dgemm_((char *)"No transpose", (char *)"No transpose", m, &i__1, k, &
			    c_b25, &work[work_offset], ldwork, &v[(*k + 1) * 
			    v_dim1 + 1], ldv, &c_b14, &c__[(*k + 1) * c_dim1 
			    + 1], ldc, (ftnlen)12, (ftnlen)12);
		}

/*              W := W * V1 */

		dtrmm_((char *)"Right", (char *)"Upper", (char *)"No transpose", (char *)"Unit", m, k, &c_b14,
			 &v[v_offset], ldv, &work[work_offset], ldwork, (
			ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + j * c_dim1] -= work[i__ + j * work_dim1];
/* L170: */
		    }
/* L180: */
		}

	    }

	} else {

/*           Let  V =  ( V1  V2 )    (V2: last K columns) */
/*           where  V2  is unit lower triangular. */

	    if (lsame_(side, (char *)"L", (ftnlen)1, (ftnlen)1)) {

/*              Form  H * C  or  H**T * C  where  C = ( C1 ) */
/*                                                    ( C2 ) */

/*              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK) */

/*              W := C2**T */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    dcopy_(n, &c__[*m - *k + j + c_dim1], ldc, &work[j * 
			    work_dim1 + 1], &c__1);
/* L190: */
		}

/*              W := W * V2**T */

		dtrmm_((char *)"Right", (char *)"Lower", (char *)"Transpose", (char *)"Unit", n, k, &c_b14, &
			v[(*m - *k + 1) * v_dim1 + 1], ldv, &work[work_offset]
			, ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)4);
		if (*m > *k) {

/*                 W := W + C1**T * V1**T */

		    i__1 = *m - *k;
		    dgemm_((char *)"Transpose", (char *)"Transpose", n, k, &i__1, &c_b14, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &c_b14, &
			    work[work_offset], ldwork, (ftnlen)9, (ftnlen)9);
		}

/*              W := W * T**T  or  W * T */

		dtrmm_((char *)"Right", (char *)"Lower", transt, (char *)"Non-unit", n, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - V**T * W**T */

		if (*m > *k) {

/*                 C1 := C1 - V1**T * W**T */

		    i__1 = *m - *k;
		    dgemm_((char *)"Transpose", (char *)"Transpose", &i__1, n, k, &c_b25, &v[
			    v_offset], ldv, &work[work_offset], ldwork, &
			    c_b14, &c__[c_offset], ldc, (ftnlen)9, (ftnlen)9);
		}

/*              W := W * V2 */

		dtrmm_((char *)"Right", (char *)"Lower", (char *)"No transpose", (char *)"Unit", n, k, &c_b14,
			 &v[(*m - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)
			12, (ftnlen)4);

/*              C2 := C2 - W**T */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[*m - *k + j + i__ * c_dim1] -= work[i__ + j * 
				work_dim1];
/* L200: */
		    }
/* L210: */
		}

	    } else if (lsame_(side, (char *)"R", (ftnlen)1, (ftnlen)1)) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 ) */

/*              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK) */

/*              W := C2 */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    dcopy_(m, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &work[
			    j * work_dim1 + 1], &c__1);
/* L220: */
		}

/*              W := W * V2**T */

		dtrmm_((char *)"Right", (char *)"Lower", (char *)"Transpose", (char *)"Unit", m, k, &c_b14, &
			v[(*n - *k + 1) * v_dim1 + 1], ldv, &work[work_offset]
			, ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)4);
		if (*n > *k) {

/*                 W := W + C1 * V1**T */

		    i__1 = *n - *k;
		    dgemm_((char *)"No transpose", (char *)"Transpose", m, k, &i__1, &c_b14, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &c_b14, &
			    work[work_offset], ldwork, (ftnlen)12, (ftnlen)9);
		}

/*              W := W * T  or  W * T**T */

		dtrmm_((char *)"Right", (char *)"Lower", trans, (char *)"Non-unit", m, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - W * V */

		if (*n > *k) {

/*                 C1 := C1 - W * V1 */

		    i__1 = *n - *k;
		    dgemm_((char *)"No transpose", (char *)"No transpose", m, &i__1, k, &
			    c_b25, &work[work_offset], ldwork, &v[v_offset], 
			    ldv, &c_b14, &c__[c_offset], ldc, (ftnlen)12, (
			    ftnlen)12);
		}

/*              W := W * V2 */

		dtrmm_((char *)"Right", (char *)"Lower", (char *)"No transpose", (char *)"Unit", m, k, &c_b14,
			 &v[(*n - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)
			12, (ftnlen)4);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c__[i__ + (*n - *k + j) * c_dim1] -= work[i__ + j * 
				work_dim1];
/* L230: */
		    }
/* L240: */
		}

	    }

	}
    }

    return 0;

/*     End of DLARFB */

} /* dlarfb_ */

#ifdef __cplusplus
	}
#endif
