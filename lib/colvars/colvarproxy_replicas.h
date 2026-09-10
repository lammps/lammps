// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARPROXY_REPLICAS_H
#define COLVARPROXY_REPLICAS_H


#ifdef COLVARS_MPI
#include <mpi.h>
typedef MPI_Comm replicas_mpi_comm_t;
#else
typedef void * replicas_mpi_comm_t;
#endif


/// \brief Methods for multiple-replica communication
class colvarproxy_replicas {

public:

  /// Constructor
  colvarproxy_replicas();

  /// Destructor
  virtual ~colvarproxy_replicas();

  /// Set the multiple replicas communicator
  virtual void set_replicas_mpi_communicator(replicas_mpi_comm_t comm);

  /// Indicate if multi-replica support is available and active
  virtual int check_replicas_enabled();

  /// Index of this replica
  virtual int replica_index();

  /// Total number of replicas
  virtual int num_replicas();

  /// Synchronize replica with others
  virtual void replica_comm_barrier();

  /// Receive data from other replica
  virtual int replica_comm_recv(char* msg_data, int buf_len, int src_rep);

  /// Send data to other replica
  virtual int replica_comm_send(char* msg_data, int msg_len, int dest_rep);

protected:

  /// MPI communicator containint 1 root proc from each world
  replicas_mpi_comm_t replicas_mpi_comm;

  /// Index (rank) of this replica in the MPI implementation
  int replicas_mpi_rank = 0;

  /// Number of replicas in the MPI implementation
  int replicas_mpi_num = 1;
};

#endif
