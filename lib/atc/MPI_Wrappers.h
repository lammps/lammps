#ifndef MPI_WRAPPERS_H
#define MPI_WRAPPERS_H

#include <iostream>
#include <string>
#include "mpi.h"

namespace MPI_Wrappers {

  typedef struct {double val; int rank; } DOUBLE_RANK;

  int rank(MPI_Comm comm);
  bool rank_zero(MPI_Comm comm);
  int size(MPI_Comm comm);
  bool serial(MPI_Comm comm);
  void broadcast(MPI_Comm comm, double *buf, int count = 1);
  void int_broadcast(MPI_Comm comm, int *buf, int count = 1);
  void allsum(MPI_Comm comm, void *send_buf, double *rec_buf, int count = 1);
  void int_allsum(MPI_Comm comm, void *send_buf, int *rec_buf, int count = 1);
  void int_scansum(MPI_Comm comm, int *send_buf, int *rec_buf, int count = 1);
  void allmax(MPI_Comm comm, double *send_buf, double *rec_buf, int count = 1);
  void int_allmax(MPI_Comm comm, int *send_buf, int *rec_buf, int count = 1);
  void allmin(MPI_Comm comm, double *send_buf, double *rec_buf, int count = 1);
  void int_allmin(MPI_Comm comm, int *send_buf, int *rec_buf, int count = 1);
  int rank_min(MPI_Comm comm, double *send_buf, double *rec_buf, int count);
  void int_recv(MPI_Comm comm, int *recv_buf, int max_size,int iproc);
  void recv(MPI_Comm comm, double *recv_buf, int max_size,int iproc);
  void int_send(MPI_Comm comm, int *send_buf,int send_size);
  void send(MPI_Comm comm, double *send_buf,int send_size);
  void allgatherv(MPI_Comm comm, double *send_buf, int send_count,
                  double *rec_buf, int *rec_counts, int *displacements);
  void gather(MPI_Comm comm, double send, double * recv);
  void int_allgather(MPI_Comm comm, int send, int* recv);
  void logical_or(MPI_Comm comm, void *send_buf, int *rec_buf, int count = 1);
  void barrier(MPI_Comm comm);
  void stop(MPI_Comm comm, std::string msg="");
  void int_scatter(MPI_Comm comm, int *send_buf, int *rec_buf, int count = 1);

//  void sparse_allsum(MPI_Comm comm, SparseMatrix<double> &toShare);

  void print_msg(MPI_Comm comm, std::string msg);
  void print_msg_once(MPI_Comm comm,std::string msg,bool prefix=true,bool endline=true);

}

#endif
