!	This program calculates elastic constants
!	Must be used with in.comb-elastic script, and reads output log.lammps
!	Written by T-R Shan, MSE, UF, Feb 2010

	program main
	implicit none
	
	integer :: i,j,k,l
	real*8  :: box(6,11),force(6,11),vol(11),eps(11),strbox(6,11)
	real*8  :: sumx,sumx2,sumy,sumxy,c11,c12,c13,c14,c15,c16,c33,c44,c66,bulk,shear
	character*7  :: header1

	open(5,status='old',file='log.lammps')

	do i=1,500
	 read(5,'(a7)') header1
	 if(header1.eq.'Step Lx')then
	  do j=1,11
	   read(5,*) k,(box(l,j),l=1,6),(force(l,j),l=1,6)
	   vol(j)=box(1,j)*box(2,j)*box(3,j)*dsqrt(1-(box(4,j)/box(1,j))*(box(4,j)/box(1,j))-(box(6,j)/box(3,j))*(box(6,j)/box(3,j)-(box(5,j)/box(3,j))*(box(5,j)/box(3,j))))
	   eps(j)=(box(1,j)-box(1,1))/box(1,1)
	   do l=1,6
	    strbox(l,j)=force(l,j)/vol(j)/1.0d4
	   enddo 
	  enddo
	  exit
	 endif
	enddo

!!	C11
	sumx = 0.0; sumx2 = 0.0; sumy = 0.0; sumxy = 0.0
	do j=2,11
	  sumx = sumx + eps(j)
	  sumx2= sumx2+ eps(j)*eps(j)
	  sumy = sumy + strbox(1,j)
	  sumxy= sumxy+ strbox(1,j)*eps(j)
	enddo
	c11 = (sumxy-sumx*sumy/10.0)/(sumx2-sumx*sumx/10.0)
!!	C12
	sumy = 0.0; sumxy = 0.0
	do j=2,11
	  sumy = sumy + strbox(2,j)
	  sumxy= sumxy+ strbox(2,j)*eps(j)
	enddo
	c12 = (sumxy-sumx*sumy/10.0)/(sumx2-sumx*sumx/10.0)
!!	C13
	sumy = 0.0; sumxy = 0.0
	do j=2,11
	  sumy = sumy + strbox(3,j)
	  sumxy= sumxy+ strbox(3,j)*eps(j)
	enddo
	c13 = (sumxy-sumx*sumy/10.0)/(sumx2-sumx*sumx/10.0)
!!	C14
	sumy = 0.0; sumxy = 0.0
	do j=2,11
	  sumy = sumy + strbox(4,j)
	  sumxy= sumxy+ strbox(4,j)*eps(j)
	enddo
	c14 = (sumxy-sumx*sumy/10.0)/(sumx2-sumx*sumx/10.0)
!!	C15
	sumy = 0.0; sumxy = 0.0
	do j=2,11
	  sumy = sumy + strbox(5,j)
	  sumxy= sumxy+ strbox(5,j)*eps(j)
	enddo
	c15 = (sumxy-sumx*sumy/10.0)/(sumx2-sumx*sumx/10.0)
!!	C16
	sumy = 0.0; sumxy = 0.0
	do j=2,11
	  sumy = sumy + strbox(6,j)
	  sumxy= sumxy+ strbox(6,j)*eps(j)
	enddo
	c16 = (sumxy-sumx*sumy/10.0)/(sumx2-sumx*sumx/10.0)
	
!!	C33
	do i=1,500
	 read(5,'(a7)') header1
	 if(header1.eq.'Step Lx')then
	  do j=1,11
	   read(5,*) k,(box(l,j),l=1,6),(force(l,j),l=1,6)
	   vol(j)=box(1,j)*box(2,j)*box(3,j)*dsqrt(1-(box(4,j)/box(1,j))*(box(4,j)/box(1,j))-(box(6,j)/box(3,j))*(box(6,j)/box(3,j)-(box(5,j)/box(3,j))*(box(5,j)/box(3,j))))
	   eps(j)=(box(3,j)-box(3,1))/box(3,1)
	   do l=1,6
	    strbox(l,j)=force(l,j)/vol(j)/1.0d4
	   enddo 
	  enddo
	  exit
	 endif
	enddo

	sumx = 0.0; sumx2 = 0.0; sumy = 0.0; sumxy = 0.0
	do j=2,11
	  sumx = sumx + eps(j)
	  sumx2= sumx2+ eps(j)*eps(j)
	  sumy = sumy + strbox(3,j)
	  sumxy= sumxy+ strbox(3,j)*eps(j)
	enddo
	c33 = (sumxy-sumx*sumy/10.0)/(sumx2-sumx*sumx/10.0)

!!	C44
	do i=1,500
	 read(5,'(a7)') header1
	 if(header1.eq.'Step Lx')then
	  do j=1,11
	   read(5,*) k,(box(l,j),l=1,6),(force(l,j),l=1,6)
	   vol(j)=box(1,j)*box(2,j)*box(3,j)*dsqrt(1.0-(box(4,j)/box(1,j))*(box(4,j)/box(1,j))-(box(6,j)/box(3,j))*(box(6,j)/box(3,j)-(box(5,j)/box(3,j))*(box(5,j)/box(3,j))))
	   do l=1,6
	    strbox(l,j)=force(l,j)/vol(j)/1.0d4
	   enddo 
	  enddo
	  exit
	 endif
	enddo
	
	sumx = 0.0; sumx2 = 0.0; sumy = 0.0; sumxy = 0.0
	do j=2,11
	  eps(j)=asin(box(6,j)/box(3,j))
	  sumx = sumx + eps(j)
	  sumx2= sumx2+ eps(j)*eps(j)
	  sumy = sumy + strbox(6,j)
	  sumxy= sumxy+ strbox(6,j)*eps(j)
	enddo
	c44 = (sumxy-sumx*sumy/10.0)/(sumx2-sumx*sumx/10.0)

!!	C66
	do i=1,500
	 read(5,'(a7)') header1
	 if(header1.eq.'Step Lx')then
	  do j=1,11
	   read(5,*) k,(box(l,j),l=1,6),(force(l,j),l=1,6)
	   vol(j)=box(1,j)*box(2,j)*box(3,j)*dsqrt(1.0-(box(4,j)/box(1,j))*(box(4,j)/box(1,j))-(box(6,j)/box(3,j))*(box(6,j)/box(3,j)-(box(5,j)/box(3,j))*(box(5,j)/box(3,j))))
	   do l=1,6
	    strbox(l,j)=force(l,j)/vol(j)/1.0d4
	   enddo 
	  enddo
	  exit
	 endif
	enddo

	sumx = 0.0; sumx2 = 0.0; sumy = 0.0; sumxy = 0.0
	do j=2,11
	  eps(j)=asin(box(4,j)/box(1,j))
	  sumx = sumx + eps(j)
	  sumx2= sumx2+ eps(j)*eps(j)
	  sumy = sumy + strbox(4,j)
	  sumxy= sumxy+ strbox(4,j)*eps(j)
	enddo
	c66 = (sumxy-sumx*sumy/10.0)/(sumx2-sumx*sumx/10.0)

	bulk=(c11*2.0+c33)/9.0+2.0*(c12+c13*2.0)/9.0
	shear=(2.*c11+c33-c12-2.*c13)/15.+(2.*c44+c66)/5.

	write(*,*)
	write(*,*)'Elastic constants (GPa):'
	write(*,101) 'C11 = ',c11
	write(*,101) 'C12 = ',c12
	write(*,101) 'C13 = ',c13
	write(*,101) 'C14 = ',c14
	write(*,101) 'C15 = ',c15
	write(*,101) 'C16 = ',c16
	write(*,101) 'C33 = ',c33
	write(*,101) 'C44 = ',c44
	write(*,101) 'C66 = ',c66
	write(*,101) 'B = ',bulk
	write(*,101) 'G = ',shear

	write(5,*)
	write(5,*)'Elastic constants (GPa):'
	write(5,101) 'C11 = ',c11
	write(5,101) 'C12 = ',c12
	write(5,101) 'C13 = ',c13
	write(5,101) 'C14 = ',c14
	write(5,101) 'C15 = ',c15
	write(5,101) 'C16 = ',c16
	write(5,101) 'C33 = ',c33
	write(5,101) 'C44 = ',c44
	write(5,101) 'C66 = ',c66
	write(5,101) 'B = ',bulk
	write(5,101) 'G = ',shear

101	format(a6,f11.3)

	stop
	end
