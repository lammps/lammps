MODULE MPI
  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: MPI_COMM_WORLD=0
  INTEGER, PARAMETER :: MPI_SUCCESS=0

  PUBLIC :: MPI_COMM_WORLD, MPI_SUCCESS, &
      mpi_comm_split

CONTAINS
  
  SUBROUTINE mpi_comm_split(comm,color,key,newcomm,ierr)
    INTEGER, INTENT(in)  :: comm,color,key
    INTEGER, INTENT(out) :: newcomm,ierr

    newcomm = comm + 1
    ierr = 0
  END SUBROUTINE mpi_comm_split
END MODULE MPI
