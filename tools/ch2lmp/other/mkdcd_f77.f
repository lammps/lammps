c -------------------------------------------------------------------------
c Code converts LAMMPS atom dumps to CHARMM .dcd files
c Paul Crozier, SNL, 2003
c -------------------------------------------------------------------------
c -------------------------------------------------------------------------

      program mkdcd
      
      parameter 		(max_atoms = 100000)
      common/int_variables/	n_con_files, n_con_files_to_skip, 
     1 n_frames_per_con_file, n_frames_per_dump_to_dcd,
     2 max_dcd_dumps_per_dcd_file, n_dcd_files_to_skip, n_frames,
     3 n_frames_to_skip, ith_frame, ith_dump_to_dcd, n_timesteps, 
     4 n_atoms, n_frozen_atoms, n_atoms_dump, n_atoms_recenter,
     5 i_recenter_flag
      common/real_variables/  xprd, yprd, zprd, box(2,3),  
     1 x(3,max_atoms)
	  common/char_variables/ con_file_path, dcd_file_path
      character*76 con_file_path
      character*76 dcd_file_path
 
      call read_in_mkdcd

      do ith_frame = n_frames_to_skip+1, n_frames
         call find_config
         call read_config
         if (i_recenter_flag == 1) call recenter
         if (mod(ith_frame,n_frames_per_dump_to_dcd) == 0) call mk_dcd 
      enddo

      write(6,*) 'Done.'
      stop
      end

c -------------------------------------------------------------------------

      subroutine find_config
      parameter 		(max_atoms = 100000)
      common/int_variables/	n_con_files, n_con_files_to_skip, 
     1 n_frames_per_con_file, n_frames_per_dump_to_dcd,
     2 max_dcd_dumps_per_dcd_file, n_dcd_files_to_skip, n_frames,
     3 n_frames_to_skip, ith_frame, ith_dump_to_dcd, n_timesteps, 
     4 n_atoms, n_frozen_atoms, n_atoms_dump, n_atoms_recenter,
     5 i_recenter_flag
      common/real_variables/  xprd, yprd, zprd, box(2,3),  
     1 x(3,max_atoms)
	  common/char_variables/ con_file_path, dcd_file_path
      character*76 con_file_path
      character*76 dcd_file_path
 
      integer l,m,n,i,j,itag

      real*8  buf(8)

      if (mod((ith_frame-1),n_frames_per_con_file) == 0) then
        n = (ith_frame-1)/n_frames_per_con_file + 1
        write(6,*) 'On config file # ', n
        close(21)
        if (n < 10) then
          open(21,file="conf"//char(48+n),status='old')
        elseif (n < 100) then
          m = n/10
          n = mod(n,10)
          open(21,file="conf"
     $      //char(48+m)//char(48+n),status='old')
        else
          l = n/100
          m = (n - 100*l)/10
          n = mod(n,10)
          open(21,file="conf"
     $      //char(48+l)//char(48+m)//char(48+n),status='old')
        endif
        rewind 21

c skip the first frame of each config file

        read(21,*) 
        read(21,*) n_timesteps
        read(21,*)
        read(21,*) n_atoms
        read(21,*)
        read(21,*) box(1,1),box(2,1)
        read(21,*) box(1,2),box(2,2)
        read(21,*) box(1,3),box(2,3)
        read(21,*)

        xprd = box(2,1) - box(1,1)
        yprd = box(2,2) - box(1,2)
        zprd = box(2,3) - box(1,3)

        if (n_atoms > max_atoms) write(6,*) "n_atoms > max_atoms"

        n_frozen_atoms = 1000000000
        do i = 1, n_atoms
           read (21,*) (buf(j),j=1,5)
           itag = nint(buf(1))
           n_frozen_atoms = min(itag,n_frozen_atoms)
        enddo

        n_frozen_atoms = n_frozen_atoms - 1

        i_recenter_flag = 0        
        if (n_atoms_recenter < n_atoms) i_recenter_flag = 1

      endif

      return
      end

c -------------------------------------------------------------------------

      subroutine mk_dcd
      parameter 		(max_atoms = 100000)
      common/int_variables/	n_con_files, n_con_files_to_skip, 
     1 n_frames_per_con_file, n_frames_per_dump_to_dcd,
     2 max_dcd_dumps_per_dcd_file, n_dcd_files_to_skip, n_frames,
     3 n_frames_to_skip, ith_frame, ith_dump_to_dcd, n_timesteps, 
     4 n_atoms, n_frozen_atoms, n_atoms_dump, n_atoms_recenter,
     5 i_recenter_flag
      common/real_variables/  xprd, yprd, zprd, box(2,3),  
     1 x(3,max_atoms)
	  common/char_variables/ con_file_path, dcd_file_path
      character*76 con_file_path
      character*76 dcd_file_path
      
      real*8 xtlabc(6)

      integer i

      if (mod(ith_dump_to_dcd,max_dcd_dumps_per_dcd_file) == 0) 
     &  call mk_dcd_start

      ith_dump_to_dcd = ith_dump_to_dcd + 1

      write(6,*) 'Frame # ', ith_frame, 
     &           ', Dump to .dcd #', ith_dump_to_dcd

      xtlabc(1) = xprd
      xtlabc(2) = 0.0

      xtlabc(3) = yprd
      xtlabc(4) = 0.0
      xtlabc(5) = 0.0
      xtlabc(6) = zprd

      write(26) xtlabc
      write(26) (x(1,i), i=1,n_atoms_dump)
      write(26) (x(2,i), i=1,n_atoms_dump)
      write(26) (x(3,i), i=1,n_atoms_dump)

      return
      end

c -------------------------------------------------------------------------

      subroutine mk_dcd_start
      parameter 		(max_atoms = 100000)
      common/int_variables/	n_con_files, n_con_files_to_skip, 
     1 n_frames_per_con_file, n_frames_per_dump_to_dcd,
     2 max_dcd_dumps_per_dcd_file, n_dcd_files_to_skip, n_frames,
     3 n_frames_to_skip, ith_frame, ith_dump_to_dcd, n_timesteps, 
     4 n_atoms, n_frozen_atoms, n_atoms_dump, n_atoms_recenter,
     5 i_recenter_flag
      common/real_variables/  xprd, yprd, zprd, box(2,3),  
     1 x(3,max_atoms)
	  common/char_variables/ con_file_path, dcd_file_path
      character*76 con_file_path
      character*76 dcd_file_path
      
      character*4 hdr

      integer icntrl(20), nstr, n_dcd_file, n, m

      write(6,*) 'Creating new .dcd file . . . '

      n_dcd_file = ith_dump_to_dcd/max_dcd_dumps_per_dcd_file +
     &             n_dcd_files_to_skip + 1

      m = n_dcd_file/10
      n = mod(n_dcd_file,10)

      close(26)
      open(26,file="out"
     $   //char(48+m)//char(48+n)//'.dcd',form='unformatted')

      hdr = 'CORD'
      do n = 1,20
            icntrl(n) = 0
      enddo
      nstr = 0					! number of strings in header
      icntrl(1) = max_dcd_dumps_per_dcd_file	! number of frames in traj file
      icntrl(2) = 1				! number of steps in previous run
      icntrl(3) = 1				! frequency of saving
      icntrl(4) = max_dcd_dumps_per_dcd_file	! total number of steps
      icntrl(8) = n_atoms_dump*3 - 6			! number of degrees of freedom
      icntrl(10) = 981668463			! coded time step
      icntrl(11) = 1				! coded crystallographic group (or zero)
      icntrl(20) = 28				! CHARMM version number

      write(26) hdr, icntrl
      write(26) nstr
      write(26) n_atoms_dump

      return
      end

c -------------------------------------------------------------------------
c input data from config file

      subroutine read_config
      parameter 		(max_atoms = 100000)
      common/int_variables/	n_con_files, n_con_files_to_skip, 
     1 n_frames_per_con_file, n_frames_per_dump_to_dcd,
     2 max_dcd_dumps_per_dcd_file, n_dcd_files_to_skip, n_frames,
     3 n_frames_to_skip, ith_frame, ith_dump_to_dcd, n_timesteps, 
     4 n_atoms, n_frozen_atoms, n_atoms_dump, n_atoms_recenter,
     5 i_recenter_flag
      common/real_variables/  xprd, yprd, zprd, box(2,3),  
     1 x(3,max_atoms)
	  common/char_variables/ con_file_path, dcd_file_path
      character*76 con_file_path
      character*76 dcd_file_path
      
c local variables

      integer i,j,itag,itrue

      real*8  buf(8)

      read(21,*) 
      read(21,*) n_timesteps
      read(21,*)
      read(21,*) n_atoms
      read(21,*)
      read(21,*) box(1,1),box(2,1)
      read(21,*) box(1,2),box(2,2)
      read(21,*) box(1,3),box(2,3)
      read(21,*)
    
      xprd = box(2,1) - box(1,1)
      yprd = box(2,2) - box(1,2)
      zprd = box(2,3) - box(1,3)

      do i = 1, n_atoms
         read (21,*) (buf(j),j=1,5)
c         itag = nint(buf(1)) - n_frozen_atoms
c         x(1,itag) = buf(3)*xprd + box(1,1)
c         x(2,itag) = buf(4)*yprd + box(1,2)
c         x(3,itag) = buf(5)*zprd + box(1,3)
         x(1,i) = buf(3)*xprd + box(1,1)
         x(2,i) = buf(4)*yprd + box(1,2)
         x(3,i) = buf(5)*zprd + box(1,3)
      enddo

      return
      end

c -------------------------------------------------------------------------

c read data from in_mkdcd file

      subroutine read_in_mkdcd
      parameter 		(max_atoms = 100000)
      common/int_variables/	n_con_files, n_con_files_to_skip, 
     1 n_frames_per_con_file, n_frames_per_dump_to_dcd,
     2 max_dcd_dumps_per_dcd_file, n_dcd_files_to_skip, n_frames,
     3 n_frames_to_skip, ith_frame, ith_dump_to_dcd, n_timesteps, 
     4 n_atoms, n_frozen_atoms, n_atoms_dump, n_atoms_recenter,
     5 i_recenter_flag
      common/real_variables/  xprd, yprd, zprd, box(2,3),  
     1 x(3,max_atoms)
	  common/char_variables/ con_file_path, dcd_file_path
      character*76 con_file_path
      character*76 dcd_file_path
      
 100  format (a)

      open(22,file='in_mkdcd')
      rewind 22

      read (22,*) n_con_files
      read (22,*) n_con_files_to_skip
      read (22,*) n_frames_per_con_file
      read (22,*) n_frames_per_dump_to_dcd
      read (22,*) max_dcd_dumps_per_dcd_file
      read (22,*) n_dcd_files_to_skip
      read (22,*) n_atoms_dump
      read (22,*) n_atoms_recenter
      read (22,100) con_file_path
      read (22,100) dcd_file_path

      n_frames = n_con_files * n_frames_per_con_file
      n_frames_to_skip = n_con_files_to_skip * n_frames_per_con_file
      ith_frame = n_frames_to_skip
      ith_dump_to_dcd = 0

      close (22)

      return
      end

c -------------------------------------------------------------------------

c recenter box on the average position of the first n_atoms_recenter atoms

      subroutine recenter
      parameter 		(max_atoms = 100000)
      common/int_variables/	n_con_files, n_con_files_to_skip, 
     1 n_frames_per_con_file, n_frames_per_dump_to_dcd,
     2 max_dcd_dumps_per_dcd_file, n_dcd_files_to_skip, n_frames,
     3 n_frames_to_skip, ith_frame, ith_dump_to_dcd, n_timesteps, 
     4 n_atoms, n_frozen_atoms, n_atoms_dump, n_atoms_recenter,
     5 i_recenter_flag
      common/real_variables/  xprd, yprd, zprd, box(2,3),  
     1 x(3,max_atoms)
	  common/char_variables/ con_file_path, dcd_file_path
      character*76 con_file_path
      character*76 dcd_file_path
      
      integer i
      real*8 ave_x, ave_y, ave_z
      
      ave_x = 0.0
      ave_y = 0.0
      ave_z = 0.0
      
      do i = 1, n_atoms_recenter
         ave_x = ave_x + x(1,i)
         ave_y = ave_y + x(2,i)
         ave_z = ave_z + x(3,i)
      enddo
      
      ave_x = ave_x/float(n_atoms_recenter)
      ave_y = ave_y/float(n_atoms_recenter)
      ave_z = ave_z/float(n_atoms_recenter)
      
      do i = 1, n_atoms_dump
         x(1,i) = x(1,i) - ave_x
         x(2,i) = x(2,i) - ave_y
         x(3,i) = x(3,i) - ave_z
      enddo
      
      do i = 1, n_atoms_dump
      	x(1,i) = x(1,i) - xprd*(nint(x(1,i)/xprd))
      	x(2,i) = x(2,i) - xprd*(nint(x(2,i)/yprd))
      	x(3,i) = x(3,i) - xprd*(nint(x(3,i)/zprd))
      enddo
      
      return
      end

c -------------------------------------------------------------------------
