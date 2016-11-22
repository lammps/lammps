! The main program which runs our driver test case potentials
! 
! Copyright (C) 2013, Joshua More and Michele Ceriotti
! 
! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
! 
! The above copyright notice and this permission notice shall be included
! in all copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
! CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
! TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
!
! Currently the potentials implemented are the Lennard-Jones
! potential, the Silvera-Goldman para-hydrogen potential and
! the ideal gas (i.e. no interaction at all)

      PROGRAM DRIVER
         USE LJ
         USE SG
      IMPLICIT NONE

      ! SOCKET VARIABLES
      INTEGER, PARAMETER :: MSGLEN=12   ! length of the headers of the driver/wrapper communication protocol
      INTEGER socket, inet, port        ! socket ID & address of the server
      CHARACTER*1024 :: host

      ! COMMAND LINE PARSING
      CHARACTER*1024 :: cmdbuffer, vops
      INTEGER ccmd, vstyle
      LOGICAL verbose
      INTEGER commas(4), par_count      ! stores the index of commas in the parameter string
      DOUBLE PRECISION vpars(4)         ! array to store the parameters of the potential

      ! SOCKET COMMUNICATION BUFFERS
      CHARACTER*12 :: header
      LOGICAL :: isinit=.false., hasdata=.false.
      INTEGER cbuf
      CHARACTER*2048 :: initbuffer      ! it's unlikely a string this large will ever be passed...
      DOUBLE PRECISION, ALLOCATABLE :: msgbuffer(:)

      ! PARAMETERS OF THE SYSTEM (CELL, ATOM POSITIONS, ...)
      DOUBLE PRECISION sigma, eps, rc, rn, ks ! potential parameters
      INTEGER nat
      DOUBLE PRECISION pot
      DOUBLE PRECISION, ALLOCATABLE :: atoms(:,:), forces(:,:)
      DOUBLE PRECISION cell_h(3,3), cell_ih(3,3), virial(3,3)
      DOUBLE PRECISION volume

      ! NEIGHBOUR LIST ARRAYS
      INTEGER, DIMENSION(:), ALLOCATABLE :: n_list, index_list
      DOUBLE PRECISION init_volume, init_rc ! needed to correctly adjust the cut-off radius for variable cell dynamics
      DOUBLE PRECISION, ALLOCATABLE :: last_atoms(:,:) ! Holds the positions when the neighbour list is created
      DOUBLE PRECISION displacement ! Tracks how far each atom has moved since the last call of nearest_neighbours

      INTEGER i

      ! parse the command line parameters
      ! intialize defaults
      ccmd = 0
      inet = 1
      host = "localhost"//achar(0)
      port = 31415
      verbose = .false.
      par_count = 0
      vstyle = -1

      DO i = 1, IARGC()
         CALL GETARG(i, cmdbuffer)
         IF (cmdbuffer == "-u") THEN ! flag for unix socket
            inet = 0
            ccmd = 0
         ELSEIF (cmdbuffer == "-h") THEN ! read the hostname
            ccmd = 1
         ELSEIF (cmdbuffer == "-p") THEN ! reads the port number
            ccmd = 2
         ELSEIF (cmdbuffer == "-m") THEN ! reads the style of the potential function
            ccmd = 3
         ELSEIF (cmdbuffer == "-o") THEN ! reads the parameters
            ccmd = 4
         ELSEIF (cmdbuffer == "-v") THEN ! flag for verbose standard output
            verbose = .true.
         ELSE
            IF (ccmd == 0) THEN
               WRITE(*,*) " Unrecognized command line argument", ccmd
               WRITE(*,*) " SYNTAX: driver.x [-u] -h hostname -p port -m [gas|lj|sg|harm] -o 'comma_separated_parameters' [-v] "
               WRITE(*,*) ""
               WRITE(*,*) " For LJ potential use -o sigma,epsilon,cutoff "
               WRITE(*,*) " For SG potential use -o cutoff "
               WRITE(*,*) " For 1D harmonic oscillator use -o k "
               WRITE(*,*) " For the ideal gas, no options needed! "
               CALL EXIT(-1)
            ENDIF
            IF (ccmd == 1) THEN
               host = trim(cmdbuffer)//achar(0)
            ELSEIF (ccmd == 2) THEN
               READ(cmdbuffer,*) port
            ELSEIF (ccmd == 3) THEN
               IF (trim(cmdbuffer) == "lj") THEN
                  vstyle = 1
               ELSEIF (trim(cmdbuffer) == "sg") THEN
                  vstyle = 2
               ELSEIF (trim(cmdbuffer) == "harm") THEN
                  vstyle = 3
               ELSEIF (trim(cmdbuffer) == "gas") THEN
                  vstyle = 0  ! ideal gas
               ELSE
                  WRITE(*,*) " Unrecognized potential type ", trim(cmdbuffer)
                  WRITE(*,*) " Use -m [gas|lj|sg|harm] "
                  CALL EXIT(-1)
               ENDIF
            ELSEIF (ccmd == 4) THEN
               par_count = 1
               commas(1) = 0
               DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0) 
                  commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                  READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) vpars(par_count)
                  par_count = par_count + 1
               ENDDO
               READ(cmdbuffer(commas(par_count)+1:),*) vpars(par_count)
            ENDIF
            ccmd = 0
         ENDIF
      ENDDO

      IF (vstyle == -1) THEN
         WRITE(*,*) " Error, type of potential not specified."
         WRITE(*,*) " SYNTAX: driver.x [-u] -h hostname -p port -m [gas|lj|sg|harm] -o 'comma_separated_parameters' [-v] "
         WRITE(*,*) ""
         WRITE(*,*) " For LJ potential use -o sigma,epsilon,cutoff "
         WRITE(*,*) " For SG potential use -o cutoff "
         WRITE(*,*) " For the ideal gas, no options needed! "
         CALL EXIT(-1)
      ELSEIF (vstyle == 0) THEN
         IF (par_count /= 0) THEN
            WRITE(*,*) "Error: no initialization string needed for ideal gas."
            CALL EXIT(-1) 
         ENDIF   
         isinit = .true.
      ELSEIF (vstyle == 1) THEN
         IF (par_count /= 3) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For LJ potential use -o sigma,epsilon,cutoff "
            CALL EXIT(-1) ! Note that if initialization from the wrapper is implemented this exit should be removed.
         ENDIF   
         sigma = vpars(1)
         eps = vpars(2)
         rc = vpars(3)
         rn = rc*1.2
         isinit = .true.
      ELSEIF (vstyle == 2) THEN
         IF (par_count /= 1) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For SG potential use -o cutoff "
            CALL EXIT(-1) ! Note that if initialization from the wrapper is implemented this exit should be removed.
         ENDIF
         rc = vpars(1)
         rn = rc*1.2
         isinit = .true.
      ELSEIF (vstyle == 3) THEN
         IF (par_count /= 1) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For 1D harmonic potential use -o k "
            CALL EXIT(-1) ! Note that if initialization from the wrapper is implemented this exit should be removed.
         ENDIF
         ks = vpars(1)
         isinit = .true.
      ENDIF

      IF (verbose) THEN
         WRITE(*,*) " DRIVER - Connecting to host ", trim(host)
         IF (inet > 0) THEN
            WRITE(*,*) " on port ", port, " using an internet socket."
         ELSE
            WRITE(*,*) " using an UNIX socket."
         ENDIF
      ENDIF

      ! Calls the interface to the C sockets to open a communication channel
      CALL open_socket(socket, inet, port, host)
      nat = -1
      DO WHILE (.true.) ! Loops forever (or until the wrapper ends!)

         ! Reads from the socket one message header
         CALL readbuffer(socket, header, MSGLEN)
         IF (verbose) WRITE(*,*) " Message from server: ", trim(header)

         IF (trim(header) == "STATUS") THEN
            ! The wrapper is inquiring on what we are doing
            IF (.not. isinit) THEN
               CALL writebuffer(socket,"NEEDINIT    ",MSGLEN)  ! Signals that we need initialization data
            ELSEIF (hasdata) THEN
               CALL writebuffer(socket,"HAVEDATA    ",MSGLEN)  ! Signals that we are done computing and can return forces
            ELSE
               CALL writebuffer(socket,"READY       ",MSGLEN)  ! We are idling and eager to compute something
            ENDIF
         ELSEIF (trim(header) == "INIT") THEN     ! The driver is kindly providing a string for initialization
            CALL readbuffer(socket, cbuf, 4)
            CALL readbuffer(socket, initbuffer, cbuf)
            IF (verbose) WRITE(*,*) " Initializing system from wrapper, using ", trim(initbuffer)
            isinit=.true. ! We actually do nothing with this string, thanks anyway. Could be used to pass some information (e.g. the input parameters, or the index of the replica, from the driver
         ELSEIF (trim(header) == "POSDATA") THEN  ! The driver is sending the positions of the atoms. Here is where we do the calculation!

            ! Parses the flow of data from the socket
            CALL readbuffer(socket, cell_h,  9*8)  ! Cell matrix
            CALL readbuffer(socket, cell_ih, 9*8)  ! Inverse of the cell matrix (so we don't have to invert it every time here)

            ! The wrapper uses atomic units for everything, and row major storage.
            ! At this stage one should take care that everything is converted in the
            ! units and storage mode used in the driver.
            cell_h = transpose(cell_h)
            cell_ih = transpose(cell_ih)
            ! We assume an upper triangular cell-vector matrix
            volume = cell_h(1,1)*cell_h(2,2)*cell_h(3,3)

            CALL readbuffer(socket, cbuf, 4)       ! The number of atoms in the cell
            IF (nat < 0) THEN  ! Assumes that the number of atoms does not change throughout a simulation, so only does this once
               nat = cbuf
               IF (verbose) WRITE(*,*) " Allocating buffer and data arrays, with ", nat, " atoms"
               ALLOCATE(msgbuffer(3*nat))
               ALLOCATE(atoms(nat,3))
               ALLOCATE(forces(nat,3))
            ENDIF

            CALL readbuffer(socket, msgbuffer, nat*3*8)
            DO i = 1, nat
               atoms(i,:) = msgbuffer(3*(i-1)+1:3*i)
            ENDDO

            IF (vstyle == 0) THEN   ! ideal gas, so no calculation done
               pot = 0
               forces = 0
               virial = 0
            ELSEIF (vstyle == 3) THEN ! 1D harmonic potential, so only uses the first position variable
               pot = 0.5*ks*atoms(1,1)**2
               forces = 0
               forces(1,1) = -ks*atoms(1,1)
               virial = 0
               virial(1,1) = forces(1,1)*atoms(1,1)
            ELSE
               IF ((allocated(n_list) .neqv. .true.)) THEN
                  IF (verbose) WRITE(*,*) " Allocating neighbour lists."
                  ALLOCATE(n_list(nat*(nat-1)/2))
                  ALLOCATE(index_list(nat))
                  ALLOCATE(last_atoms(nat,3))
                  CALL nearest_neighbours(rn, nat, atoms, cell_h, cell_ih, index_list, n_list)
                  last_atoms = atoms
                  init_volume = volume
                  init_rc = rc
               ENDIF

               ! Checking to see if we need to re-calculate the neighbour list
               rc = init_rc*(volume/init_volume)**(1.0/3.0)
               DO i = 1, nat
                  CALL separation(cell_h, cell_ih, atoms(i,:), last_atoms(i,:), displacement)
                  ! Note that displacement is the square of the distance moved by atom i since the last time the neighbour list was created.
                  IF (4*displacement > (rn-rc)*(rn-rc)) THEN
                     IF (verbose) WRITE(*,*) " Recalculating neighbour lists"
                     CALL nearest_neighbours(rn, nat, atoms, cell_h, cell_ih, index_list, n_list)
                     last_atoms = atoms
                     rn = 1.2*rc
                     EXIT
                  ENDIF
               ENDDO

               IF (vstyle == 1) THEN
                  CALL LJ_getall(rc, sigma, eps, nat, atoms, cell_h, cell_ih, index_list, n_list, pot, forces, virial)
               ELSEIF (vstyle == 2) THEN
                  CALL SG_getall(rc, nat, atoms, cell_h, cell_ih, index_list, n_list, pot, forces, virial)
               ENDIF
               IF (verbose) WRITE(*,*) " Calculated energy is ", pot
            ENDIF
            hasdata = .true. ! Signal that we have data ready to be passed back to the wrapper
         ELSEIF (trim(header) == "GETFORCE") THEN  ! The driver calculation is finished, it's time to send the results back to the wrapper

            ! Data must be re-formatted (and units converted) in the units and shapes used in the wrapper
            DO i = 1, nat
               msgbuffer(3*(i-1)+1:3*i) = forces(i,:)
            ENDDO
            virial = transpose(virial)

            CALL writebuffer(socket,"FORCEREADY  ",MSGLEN)
            CALL writebuffer(socket,pot,8)  ! Writing the potential
            CALL writebuffer(socket,nat,4)  ! Writing the number of atoms
            CALL writebuffer(socket,msgbuffer,3*nat*8) ! Writing the forces
            CALL writebuffer(socket,virial,9*8)  ! Writing the virial tensor, NOT divided by the volume
            cbuf = 7 ! Size of the "extras" string
            CALL writebuffer(socket,cbuf,4) ! This would write out the "extras" string, but in this case we only use a dummy string.
            CALL writebuffer(socket,"nothing",7)

            hasdata = .false.
         ELSE
            WRITE(*,*) " Unexpected header ", header
            CALL EXIT(-1)
         ENDIF
      ENDDO
      IF (nat > 0) DEALLOCATE(atoms, forces, msgbuffer)
      END PROGRAM
