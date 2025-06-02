// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.


#include "colvarmodule.h"
#include "colvarproxy_replicas.h"


colvarproxy_replicas::colvarproxy_replicas()
{
#ifdef COLVARS_MPI
  replicas_mpi_comm = MPI_COMM_NULL;
#endif
}


colvarproxy_replicas::~colvarproxy_replicas() {}


void colvarproxy_replicas::set_replicas_mpi_communicator(replicas_mpi_comm_t comm)
{
  replicas_mpi_comm = comm;
#ifdef COLVARS_MPI
  if (comm != MPI_COMM_NULL) {
    MPI_Comm_rank(comm, &replicas_mpi_rank);
    MPI_Comm_size(comm, &replicas_mpi_num);
    cvm::log("Enabling multiple replicas: this is replica number " +
             cvm::to_str(replica_index() + 1) + " of " + cvm::to_str(num_replicas()) + ".\n");
  }
#endif
}


int colvarproxy_replicas::check_replicas_enabled()
{
#ifdef COLVARS_MPI
  if (replicas_mpi_comm != MPI_COMM_NULL) {
    return num_replicas() > 1 ? COLVARS_OK : COLVARS_ERROR;
  }
  return COLVARS_ERROR;
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy_replicas::replica_index()
{
  return replicas_mpi_rank;
}


int colvarproxy_replicas::num_replicas()
{
  return replicas_mpi_num;
}


void colvarproxy_replicas::replica_comm_barrier()
{
#ifdef COLVARS_MPI
  MPI_Barrier(replicas_mpi_comm);
#endif
}


int colvarproxy_replicas::replica_comm_recv(char *buffer, int buffer_length, int source_rank)
{
#ifdef COLVARS_MPI
  MPI_Status status;
  int retval = MPI_Recv(buffer, buffer_length, MPI_CHAR, source_rank, 0, replicas_mpi_comm, &status);
  if (retval == MPI_SUCCESS) {
    MPI_Get_count(&status, MPI_CHAR, &retval);
  } else {
    retval = 0;
  }
  return retval;
#else
  (void)buffer;
  (void)buffer_length;
  (void)source_rank;
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy_replicas::replica_comm_send(char *buffer, int buffer_length, int destination_rank)
{
#ifdef COLVARS_MPI
  int retval = MPI_Send(buffer, buffer_length, MPI_CHAR, destination_rank, 0, replicas_mpi_comm);
  if (retval == MPI_SUCCESS) {
    retval = buffer_length;
  } else {
    retval = 0;
  }
  return retval;
#else
  (void)buffer;
  (void)buffer_length;
  (void)destination_rank;
  return COLVARS_NOT_IMPLEMENTED;
#endif
}
