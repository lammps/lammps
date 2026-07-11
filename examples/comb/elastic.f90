! This program calculates elastic constants
! Must be used with in.comb-elastic script,
! and reads a LAMMPS log file from standard input:
!
! f95 -O -o elastic.x elastic.f90
! ./elastic.x < log.lammps
!
! Written by T-R Shan, MSE, UF, Feb 2010

PROGRAM main
  IMPLICIT NONE

  INTEGER, PARAMETER :: dbl=8
  INTEGER :: i,j,k,l
  REAL(kind=dbl) :: box(6,11),force(6,11),vol(11),eps(11),strbox(6,11)
  REAL(kind=dbl) :: sumx,sumx2,sumy,sumxy
  REAL(kind=dbl) :: c11,c12,c13,c14,c15,c16,c33,c44,c66
  REAL(kind=dbl) :: bulk,shear
  CHARACTER*7  :: header1

!  open(5,status='old',file='log.lammps')

  DO i=1,500
      READ(5,'(a7)') header1
      IF(header1.EQ.'Step Lx')THEN
          DO j=1,11
              READ(5,*) k,(box(l,j),l=1,6),(force(l,j),l=1,6)
              vol(j)=box(1,j)*box(2,j)*box(3,j) * dsqrt(1-(box(4,j)/box(1,j)) &
                  *(box(4,j)/box(1,j))-(box(6,j)/box(3,j))*(box(6,j)/box(3,j) &
                  -(box(5,j)/box(3,j))*(box(5,j)/box(3,j))))
              eps(j)=(box(1,j)-box(1,1))/box(1,1)
              DO l=1,6
                  strbox(l,j)=force(l,j)/vol(j)/1.0d4
              ENDDO
          ENDDO
          EXIT
      ENDIF
  ENDDO

  !!    C11
  sumx = 0.0_dbl; sumx2 = 0.0_dbl; sumy = 0.0_dbl; sumxy = 0.0_dbl
  DO j=2,11
      sumx = sumx + eps(j)
      sumx2= sumx2+ eps(j)*eps(j)
      sumy = sumy + strbox(1,j)
      sumxy= sumxy+ strbox(1,j)*eps(j)
  ENDDO
  c11 = (sumxy-sumx*sumy/10.0_dbl)/(sumx2-sumx*sumx/10.0_dbl)
  !!    C12
  sumy = 0.0_dbl; sumxy = 0.0_dbl
  DO j=2,11
      sumy = sumy + strbox(2,j)
      sumxy= sumxy+ strbox(2,j)*eps(j)
  ENDDO
  c12 = (sumxy-sumx*sumy/10.0_dbl)/(sumx2-sumx*sumx/10.0_dbl)
  !!    C13
  sumy = 0.0_dbl; sumxy = 0.0_dbl
  DO j=2,11
      sumy = sumy + strbox(3,j)
      sumxy= sumxy+ strbox(3,j)*eps(j)
  ENDDO
  c13 = (sumxy-sumx*sumy/10.0_dbl)/(sumx2-sumx*sumx/10.0_dbl)
  !!    C14
  sumy = 0.0_dbl; sumxy = 0.0_dbl
  DO j=2,11
      sumy = sumy + strbox(4,j)
      sumxy= sumxy+ strbox(4,j)*eps(j)
  ENDDO
  c14 = (sumxy-sumx*sumy/10.0_dbl)/(sumx2-sumx*sumx/10.0_dbl)
  !!    C15
  sumy = 0.0_dbl; sumxy = 0.0_dbl
  DO j=2,11
      sumy = sumy + strbox(5,j)
      sumxy= sumxy+ strbox(5,j)*eps(j)
  ENDDO
  c15 = (sumxy-sumx*sumy/10.0_dbl)/(sumx2-sumx*sumx/10.0_dbl)
  !!    C16
  sumy = 0.0_dbl; sumxy = 0.0_dbl
  DO j=2,11
      sumy = sumy + strbox(6,j)
      sumxy= sumxy+ strbox(6,j)*eps(j)
  ENDDO
  c16 = (sumxy-sumx*sumy/10.0_dbl)/(sumx2-sumx*sumx/10.0_dbl)

  !!    C33
  DO i=1,500
      READ(5,'(a7)') header1
      IF(header1.EQ.'Step Lx')THEN
          DO j=1,11
              READ(5,*) k,(box(l,j),l=1,6),(force(l,j),l=1,6)
              vol(j)=box(1,j)*box(2,j)*box(3,j) * dsqrt(1-(box(4,j)/box(1,j)) &
                  *(box(4,j)/box(1,j))-(box(6,j)/box(3,j))*(box(6,j)/box(3,j) &
                  -(box(5,j)/box(3,j))*(box(5,j)/box(3,j))))
              eps(j)=(box(3,j)-box(3,1))/box(3,1)
              DO l=1,6
                  strbox(l,j)=force(l,j)/vol(j)/1.0d4
              ENDDO
          ENDDO
          EXIT
      ENDIF
  ENDDO

  sumx = 0.0_dbl; sumx2 = 0.0_dbl; sumy = 0.0_dbl; sumxy = 0.0_dbl
  DO j=2,11
      sumx = sumx + eps(j)
      sumx2= sumx2+ eps(j)*eps(j)
      sumy = sumy + strbox(3,j)
      sumxy= sumxy+ strbox(3,j)*eps(j)
  ENDDO
  c33 = (sumxy-sumx*sumy/10.0_dbl)/(sumx2-sumx*sumx/10.0_dbl)

  !!    C44
  DO i=1,500
      READ(5,'(a7)') header1
      IF(header1.EQ.'Step Lx')THEN
          DO j=1,11
              READ(5,*) k,(box(l,j),l=1,6),(force(l,j),l=1,6)
              vol(j)=box(1,j)*box(2,j)*box(3,j)*dsqrt(1.0-(box(4,j)/box(1,j)) &
                  *(box(4,j)/box(1,j))-(box(6,j)/box(3,j))*(box(6,j)/box(3,j) &
                  -(box(5,j)/box(3,j))*(box(5,j)/box(3,j))))
              DO l=1,6
                  strbox(l,j)=force(l,j)/vol(j)/1.0d4
              ENDDO
          ENDDO
          EXIT
      ENDIF
  ENDDO

  sumx = 0.0_dbl; sumx2 = 0.0_dbl; sumy = 0.0_dbl; sumxy = 0.0_dbl
  DO j=2,11
      eps(j)=ASIN(box(6,j)/box(3,j))
      sumx = sumx + eps(j)
      sumx2= sumx2+ eps(j)*eps(j)
      sumy = sumy + strbox(6,j)
      sumxy= sumxy+ strbox(6,j)*eps(j)
  enddo
  c44 = (sumxy-sumx*sumy/10.0_dbl)/(sumx2-sumx*sumx/10.0_dbl)

  !!    C66
  DO i=1,500
      READ(5,'(a7)') header1
      IF(header1.EQ.'Step Lx')THEN
          DO j=1,11
              READ(5,*) k,(box(l,j),l=1,6),(force(l,j),l=1,6)
              vol(j)=box(1,j)*box(2,j)*box(3,j)*dsqrt(1.0-(box(4,j)/box(1,j)) &
                  *(box(4,j)/box(1,j))-(box(6,j)/box(3,j))*(box(6,j)/box(3,j) &
                  -(box(5,j)/box(3,j))*(box(5,j)/box(3,j))))
              DO l=1,6
                  strbox(l,j)=force(l,j)/vol(j)/1.0d4
              ENDDO
          ENDDO
          EXIT
      ENDIF
  ENDDO

  sumx = 0.0_dbl; sumx2 = 0.0_dbl; sumy = 0.0_dbl; sumxy = 0.0_dbl
  DO j=2,11
      eps(j)=ASIN(box(4,j)/box(1,j))
      sumx = sumx + eps(j)
      sumx2= sumx2+ eps(j)*eps(j)
      sumy = sumy + strbox(4,j)
      sumxy= sumxy+ strbox(4,j)*eps(j)
  ENDDO
  c66 = (sumxy-sumx*sumy/10.0_dbl)/(sumx2-sumx*sumx/10.0_dbl)

  bulk=(c11*2.0_dbl+c33)/9.0_dbl+2.0_dbl*(c12+c13*2.0_dbl)/9.0_dbl
  shear=(2.0_dbl*c11+c33-c12-2.0_dbl*c13)/15.0_dbl+(2.0_dbl*c44+c66)/5.0_dbl

  WRITE(*,*)
  WRITE(*,*)'Elastic constants (GPa):'
  WRITE(*,101) 'C11 = ',c11
  WRITE(*,101) 'C12 = ',c12
  WRITE(*,101) 'C13 = ',c13
  WRITE(*,101) 'C14 = ',c14
  WRITE(*,101) 'C15 = ',c15
  WRITE(*,101) 'C16 = ',c16
  WRITE(*,101) 'C33 = ',c33
  WRITE(*,101) 'C44 = ',c44
  WRITE(*,101) 'C66 = ',c66
  WRITE(*,101) 'B = ',bulk
  WRITE(*,101) 'G = ',shear

  ! WRITE(5,*)
  ! WRITE(5,*)'Elastic constants (GPa):'
  ! WRITE(5,101) 'C11 = ',c11
  ! WRITE(5,101) 'C12 = ',c12
  ! WRITE(5,101) 'C13 = ',c13
  ! WRITE(5,101) 'C14 = ',c14
  ! WRITE(5,101) 'C15 = ',c15
  ! WRITE(5,101) 'C16 = ',c16
  ! WRITE(5,101) 'C33 = ',c33
  ! WRITE(5,101) 'C44 = ',c44
  ! WRITE(5,101) 'C66 = ',c66
  ! WRITE(5,101) 'B = ',bulk
  ! WRITE(5,101) 'G = ',shear

101 FORMAT(a6,f11.3)

  STOP
END PROGRAM main
