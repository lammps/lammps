*> \brief \b ZDROT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZDROT( N, ZX, INCX, ZY, INCY, C, S )
*
*       .. Scalar Arguments ..
*       INTEGER            INCX, INCY, N
*       DOUBLE PRECISION   C, S
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         ZX( * ), ZY( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> Applies a plane rotation, where the cos and sin (c and s) are real
*> and the vectors cx and cy are complex.
*> jack dongarra, linpack, 3/11/78.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the order of the vectors cx and cy.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in,out] ZX
*> \verbatim
*>          ZX is COMPLEX*16 array, dimension at least
*>           ( 1 + ( N - 1 )*abs( INCX ) ).
*>           Before entry, the incremented array ZX must contain the n
*>           element vector cx. On exit, ZX is overwritten by the updated
*>           vector cx.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           ZX. INCX must not be zero.
*> \endverbatim
*>
*> \param[in,out] ZY
*> \verbatim
*>          ZY is COMPLEX*16 array, dimension at least
*>           ( 1 + ( N - 1 )*abs( INCY ) ).
*>           Before entry, the incremented array ZY must contain the n
*>           element vector cy. On exit, ZY is overwritten by the updated
*>           vector cy.
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>           On entry, INCY specifies the increment for the elements of
*>           ZY. INCY must not be zero.
*> \endverbatim
*>
*> \param[in] C
*> \verbatim
*>          C is DOUBLE PRECISION
*>           On entry, C specifies the cosine, cos.
*> \endverbatim
*>
*> \param[in] S
*> \verbatim
*>          S is DOUBLE PRECISION
*>           On entry, S specifies the sine, sin.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complex16_blas_level1
*
*  =====================================================================
      SUBROUTINE ZDROT( N, ZX, INCX, ZY, INCY, C, S )
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      DOUBLE PRECISION   C, S
*     ..
*     .. Array Arguments ..
      COMPLEX*16         ZX( * ), ZY( * )
*     ..
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IX, IY
      COMPLEX*16         CTEMP
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 )
     $   RETURN
      IF( INCX.EQ.1 .AND. INCY.EQ.1 ) THEN
*
*        code for both increments equal to 1
*
         DO I = 1, N
            CTEMP = C*ZX( I ) + S*ZY( I )
            ZY( I ) = C*ZY( I ) - S*ZX( I )
            ZX( I ) = CTEMP
         END DO
      ELSE
*
*        code for unequal increments or equal increments not equal
*          to 1
*
         IX = 1
         IY = 1
         IF( INCX.LT.0 )
     $      IX = ( -N+1 )*INCX + 1
         IF( INCY.LT.0 )
     $      IY = ( -N+1 )*INCY + 1
         DO I = 1, N
            CTEMP = C*ZX( IX ) + S*ZY( IY )
            ZY( IY ) = C*ZY( IY ) - S*ZX( IX )
            ZX( IX ) = CTEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
*
*     End of ZDROT
*
      END
