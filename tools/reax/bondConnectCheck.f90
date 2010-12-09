!# DEC.9, 2010
!# HLL
!# NCSU
!#
!# This is a program to read the output from 'fix reax/bond', TPRD, Lammps
!# The output is saved into file "bonds.reax", where each image is divided
!# into three parts: 
!# 
!# (1) Head, 7 Lines;
!# (2) Body, No._of_atom Lines;
!# (3) Tail, 1 Line
!# 
!# The total number of images is related with the output frequence and number of iterations. 
!# In this case, it is "number of iteration+1".
!#
!# Each line in Body part is made up of the following parameters:
!# id, type, nb, id_1, id_2, ... id_nb, mol, bo_1, bo_2, ... bo_nb, abo, nlp, q
!# abo = atomic bond order
!# nlp = number of lone pairs
!# q = atomic charge
!#
!# PLEASE DOUBLE CHECK YOUR OWN LAMMPS INPUT SCRIPT & OUTPUT AND MAKE CORRESPONDING CHSNGES

program main
implicit none

integer I, J, K, L
integer image, natom
integer headline, tailline
integer id, atype, nb, bd1, bd2, bd3, bd4, mol
double precision bo1, bo2, bo3, bo4, abo, nlp, q

open (unit=10, file='bonds.reax')

open (unit=20, file='N129.txt', status='unknown')
open (unit=21, file='N133.txt', status='unknown')
open (unit=22, file='N137.txt', status='unknown')
open (unit=23, file='N141.txt', status='unknown')
open (unit=24, file='N145.txt', status='unknown')
open (unit=25, file='N149.txt', status='unknown')
open (unit=26, file='N153.txt', status='unknown')
open (unit=27, file='N157.txt', status='unknown')

open (unit=30, file='reactionRecord.txt', status='unknown')

!# Make changes accordingly.
image = 1
headline = 7
tailline = 1
natom = 384

do I = 1, image+1

! Skip the head part
  do J = 1, headline
  read(10,*)
  end do 
  
! Each image has 'natom' lines
    do K = 1, natom
    
! read in the first three number each line to determine:
! (1) what type of atom it is, atype
! the correspondence in Lammps: 1-C, 2-H, 3-O, 4-N, 5-S
! (2) how many bonds it has, nb
! this 'nb' determines the following bond_link information & bond_order paramaters of the same line

  read(10,*) id, atype, nb

! TEST  
! write(*,*) id, atype, nb
  
  if (atype .eq. 4) then 
  
  backspace 10
  
! Should have some easier way to replace this "IF", I am just toooo lazy.
! Thanks to the fact that the maximum number of bonds is 4. ^-^
!??? is it possible that nb = 0 ??? KEEP THAT IN MIND.

  if (nb.eq.0) then 
  
  	read(10,*) id, atype, nb, mol, abo, nlp, q

  	if (id .eq. 129) then 
  	write(20, 200) id, atype, nb, mol, abo, nlp, q
  	elseif (id .eq. 133) then 
  	write(21, 200) id, atype, nb, mol, abo, nlp, q
  	elseif (id .eq. 137) then 
  	write(22, 200) id, atype, nb, mol, abo, nlp, q
  	elseif (id .eq. 141) then 
  	write(23, 200) id, atype, nb, mol, abo, nlp, q
  	elseif (id .eq. 145) then 
  	write(24, 200) id, atype, nb, mol, abo, nlp, q
  	elseif (id .eq. 149) then 
  	write(25, 200) id, atype, nb, mol, abo, nlp, q
  	elseif (id .eq. 153) then 
  	write(26, 200) id, atype, nb, mol, abo, nlp, q
  	elseif (id .eq. 157) then 
  	write(27, 200) id, atype, nb, mol, abo, nlp, q
  	200 format(4I4, 3f14.3)
    endif 
    
! If bd .ne. 3, it measn reaction is happening to Nitrogen atom.
    write (30, 300)  I, id, atype, nb, mol, abo, nlp, q
    300 format(5I4, 3f14.3)

    
  elseif (nb.eq.1) then
  
    read(10,*) id, atype, nb, bd1, mol, bo1, abo, nlp, q

   	if (id .eq. 129) then 
  	write(20, 201) id, atype, nb, bd1, mol, bo1, abo, nlp, q
  	elseif (id .eq. 133) then 
  	write(21, 201) id, atype, nb, bd1, mol, bo1, abo, nlp, q
  	elseif (id .eq. 137) then 
  	write(22, 201) id, atype, nb, bd1, mol, bo1, abo, nlp, q
  	elseif (id .eq. 141) then 
  	write(23, 201) id, atype, nb, bd1, mol, bo1, abo, nlp, q
  	elseif (id .eq. 145) then 
  	write(24, 201) id, atype, nb, bd1, mol, bo1, abo, nlp, q
  	elseif (id .eq. 149) then 
  	write(25, 201) id, atype, nb, bd1, mol, bo1, abo, nlp, q
  	elseif (id .eq. 153) then 
  	write(26, 201) id, atype, nb, bd1, mol, bo1, abo, nlp, q
  	elseif (id .eq. 157) then 
  	write(27, 201) id, atype, nb, bd1, mol, bo1, abo, nlp, q
  	201 format(5I4, 4f14.3) 
  	endif 
  	
 ! If bd .ne. 3, it measn reaction is happening to Nitrogen atom.
    write (30, 301)  I, id, atype, nb, bd1, mol, bo1, abo, nlp, q
    301 format(6I4, 4f14.3)
    
  elseif (nb.eq.2) then
  
    read(10,*) id, atype, nb, bd1, bd2, mol, bo1, bo2, abo, nlp, q
  
   	if (id .eq. 129) then 
  	write(20, 202) id, atype, nb, bd1, bd2, mol, bo1, bo2, abo, nlp, q
  	elseif (id .eq. 133) then 
  	write(21, 202) id, atype, nb, bd1, bd2, mol, bo1, bo2, abo, nlp, q
  	elseif (id .eq. 137) then 
  	write(22, 202) id, atype, nb, bd1, bd2, mol, bo1, bo2, abo, nlp, q
  	elseif (id .eq. 141) then 
  	write(23, 202) id, atype, nb, bd1, bd2, mol, bo1, bo2, abo, nlp, q
  	elseif (id .eq. 145) then 
  	write(24, 202) id, atype, nb, bd1, bd2, mol, bo1, bo2, abo, nlp, q
  	elseif (id .eq. 149) then 
  	write(25, 202) id, atype, nb, bd1, bd2, mol, bo1, bo2, abo, nlp, q
  	elseif (id .eq. 153) then 
  	write(26, 202) id, atype, nb, bd1, bd2, mol, bo1, bo2, abo, nlp, q
  	elseif (id .eq. 157) then 
  	write(27, 202) id, atype, nb, bd1, bd2, mol, bo1, bo2, abo, nlp, q
  	202 format(6I4, 5f14.3) 
  	endif 
  	
! If bd .ne. 3, it measn reaction is happening to Nitrogen atom.
    write (30, 302)  I, id, atype, nb, bd1, bd2, mol, bo1, bo2, abo, nlp, q
    302 format(7I4, 5f14.3)
  	
  elseif (nb.eq.3) then
  
    read(10,*) id, atype, nb, bd1, bd2, bd3, mol, bo1, bo2, bo3, abo, nlp, q

  
   	if (id .eq. 129) then 
  	write(20, 203) id, atype, nb, bd1, bd2, bd3, mol, bo1, bo2, bo3, abo, nlp, q
  	elseif (id .eq. 133) then 
  	write(21, 203) id, atype, nb, bd1, bd2, bd3, mol, bo1, bo2, bo3, abo, nlp, q
  	elseif (id .eq. 137) then 
  	write(22, 203) id, atype, nb, bd1, bd2, bd3, mol, bo1, bo2, bo3, abo, nlp, q
  	elseif (id .eq. 141) then 
  	write(23, 203) id, atype, nb, bd1, bd2, bd3, mol, bo1, bo2, bo3, abo, nlp, q
  	elseif (id .eq. 145) then 
  	write(24, 203) id, atype, nb, bd1, bd2, bd3, mol, bo1, bo2, bo3, abo, nlp, q
  	elseif (id .eq. 149) then 
  	write(25, 203) id, atype, nb, bd1, bd2, bd3, mol, bo1, bo2, bo3, abo, nlp, q
  	elseif (id .eq. 153) then 
  	write(26, 203) id, atype, nb, bd1, bd2, bd3, mol, bo1, bo2, bd3, abo, nlp, q
  	elseif (id .eq. 157) then 
  	write(27, 203) id, atype, nb, bd1, bd2, bd3, mol, bo1, bo2, bo3, abo, nlp, q
  	203 format(7I4, 6f14.3) 
    endif 
  
 elseif (nb.eq.4) then
 
   	read(10,*) id, atype, nb, bd1, bd2, bd3, bd4, mol, bo1, bo2, bo3, bo4, abo, nlp, q

   	if (id .eq. 129) then 
  	write(20, 204) id, atype, nb, bd1, bd2, bd3, bd4, mol, bo1, bo2, bo3, bo4, abo, nlp, q
  	elseif (id .eq. 133) then 
  	write(21, 204) id, atype, nb, bd1, bd2, bd3, bd4, mol, bo1, bo2, bo3, bo4, abo, nlp, q
  	elseif (id .eq. 137) then 
  	write(22, 204) id, atype, nb, bd1, bd2, bd3, bd4, mol, bo1, bo2, bo3, bo4, abo, nlp, q
  	elseif (id .eq. 141) then 
  	write(23, 204) id, atype, nb, bd1, bd2, bd3, bd4, mol, bo1, bo2, bo3, bo4, abo, nlp, q
  	elseif (id .eq. 145) then 
  	write(24, 204) id, atype, nb, bd1, bd2, bd3, bd4, mol, bo1, bo2, bo3, bo4, abo, nlp, q
  	elseif (id .eq. 149) then 
  	write(25, 204) id, atype, nb, bd1, bd2, bd3, bd4, mol, bo1, bo2, bo3, bo4, abo, nlp, q
  	elseif (id .eq. 153) then 
  	write(26, 204) id, atype, nb, bd1, bd2, bd3, bd4, mol, bo1, bo2, bd3, bo4, abo, nlp, q
  	elseif (id .eq. 157) then 
  	write(27, 204) id, atype, nb, bd1, bd2, bd3, bd4, mol, bo1, bo2, bo3, bo4, abo, nlp, q
  	204 format(8I4, 7f14.3) 
    endif
    
! If bd .ne. 3, it measn reaction is happening to Nitrogen atom.
    write (30, 304)  I, id, atype, nb, bd1, bd2, bd3, bd4, mol, bo1, bo2, bo3, bo4, abo, nlp, q
    304 format(9I4, 7f14.3)

! Corresponding to "if (nb.eq.0) then "
 
  endif
  
! Corresponding to "if (atype .eq. 4) then"
  endif
  
   
  enddo
  
  do L =1,tailline
  read(10,*)
  enddo
  
  enddo
  
  end program main
  
  
