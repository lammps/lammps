#include "MPI_Wrappers.h"
#include "Utility.h"
using ATC_Utility::to_string;
#include "ATC_Error.h"
using ATC::ATC_Error;
using std::cout;
using std::string;
#ifdef ISOLATE_FE
#include "Matrix.h"
using ATC_Matrix::SparseMatrix;
#endif

namespace MPI_Wrappers {

  int rank(MPI_Comm comm)
  {
    int rank;
    MPI_Comm_rank(comm, &rank);
    return rank;
  }
  bool rank_zero(MPI_Comm comm) { return rank(comm)==0;}

  int size(MPI_Comm comm)
  {
    int size;
    MPI_Comm_size(comm, &size);
    return size;
  }
  bool serial(MPI_Comm comm) { return size(comm) == 0; }

  void broadcast(MPI_Comm comm, double *buf, int count)
  {
    int error = MPI_Bcast(buf, count, MPI_DOUBLE, 0, comm);
    if (error != MPI_SUCCESS) throw ATC_Error("error in broadcast "+to_string(error));
  }

  void int_broadcast(MPI_Comm comm, int *buf, int count)
  {
    int error = MPI_Bcast(buf, count, MPI_INT, 0, comm);
    if (error != MPI_SUCCESS) throw ATC_Error("error in broadcast "+to_string(error));
  }

  void allsum(MPI_Comm comm, void *send_buf, double *rec_buf, int count)
  {
    int error = MPI_Allreduce(send_buf, rec_buf, count, MPI_DOUBLE, MPI_SUM,
                              comm);
    if (error != MPI_SUCCESS) throw ATC_Error("error in allsum "+to_string(error));
  }

  void int_allsum(MPI_Comm comm, void *send_buf, int *rec_buf, int count)
  {
    int error = MPI_Allreduce(send_buf, rec_buf, count, MPI_INT, MPI_SUM,
                              comm);
    if (error != MPI_SUCCESS) throw ATC_Error("error in int_allsum "+to_string(error));
  }

  void int_scansum(MPI_Comm comm, int *send_buf, int *rec_buf, int count)
  {
    int error = MPI_Scan(send_buf, rec_buf, count, MPI_INT, MPI_SUM,
                              comm);
    if (error != MPI_SUCCESS) throw ATC_Error("error in int_scansum "+to_string(error));
  }


  void allmax(MPI_Comm comm, double *send_buf, double *rec_buf, int count)
  {
    int error = MPI_Allreduce(send_buf, rec_buf, count, MPI_DOUBLE, MPI_MAX,
                              comm);
    if (error != MPI_SUCCESS) throw ATC_Error("error in allmax "+to_string(error));
  }

  void int_allmax(MPI_Comm comm, int *send_buf, int *rec_buf, int count)
  {
    int error = MPI_Allreduce(send_buf, rec_buf, count, MPI_INT, MPI_MAX,
                              comm);
    if (error != MPI_SUCCESS) throw ATC_Error("error in int_allmax "+to_string(error));
  }

  void allmin(MPI_Comm comm, double *send_buf, double *rec_buf, int count)
  {
    int error = MPI_Allreduce(send_buf, rec_buf, count, MPI_DOUBLE, MPI_MIN,
                              comm);
    if (error != MPI_SUCCESS) throw ATC_Error("error in allmax "+to_string(error));
  }

  void int_allmin(MPI_Comm comm, int *send_buf, int *rec_buf, int count)
  {
    int error = MPI_Allreduce(send_buf, rec_buf, count, MPI_INT, MPI_MIN,
                              comm);
    if (error != MPI_SUCCESS) throw ATC_Error("error in int_allmax "+to_string(error));
  }

  int rank_min(MPI_Comm comm, double *send_buf, double *rec_buf, int count)
  {
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    DOUBLE_RANK in[count],out[count];
    for (int i = 0; i < count; i++) {
      in[i].val = send_buf[i];
      in[i].rank = myRank;
    }
    int error = MPI_Allreduce(in, out, count, MPI_DOUBLE_INT, MPI_MINLOC,
                              comm);
    if (error != MPI_SUCCESS) throw ATC_Error("error in rank_min "+to_string(error));
    for (int i = 0; i < count; i++) {
      rec_buf[i] = out[i].val;
    }
    return out[0].rank;
  }

  void int_recv(MPI_Comm comm, int *recv_buf, int max_size, int iproc)
  {
    MPI_Status status;
    MPI_Request request;
    int tmp, error, recv_size;
    error = MPI_Irecv(recv_buf,max_size,MPI_INT,iproc,0,comm,&request);
    error = error && MPI_Send(&tmp,0,MPI_INT,iproc,0,comm);
    error = error && MPI_Wait(&request,&status);
    error = error && MPI_Get_count(&status,MPI_DOUBLE,&recv_size);
    if (error != MPI_SUCCESS) throw ATC_Error("error in int_recv "+to_string(error));
  }

  void recv(MPI_Comm comm, double *recv_buf, int max_size,int iproc)
  {
    MPI_Status status;
    MPI_Request request;
    int tmp, error, recv_size;
    error = MPI_Irecv(recv_buf,max_size,MPI_DOUBLE,iproc,0,comm,&request);
    error = error && MPI_Send(&tmp,0,MPI_INT,iproc,0,comm);
    error = error && MPI_Wait(&request,&status);
    error = error && MPI_Get_count(&status,MPI_DOUBLE,&recv_size);
    if (error != MPI_SUCCESS) throw ATC_Error("error in recv "+to_string(error));
  }

  void int_send(MPI_Comm comm, int *send_buf,int send_size)
  {
    MPI_Status status;
    int tmp, error;
    error = MPI_Recv(&tmp,0,MPI_INT,0,0,comm,&status);
    error = error && MPI_Rsend(send_buf,send_size,MPI_INT,0,0,comm);
    if (error != MPI_SUCCESS) throw ATC_Error("error in int_send "+to_string(error));
  }

  void send(MPI_Comm comm, double *send_buf,int send_size)
  {
    MPI_Status status;
    int tmp, error;
    error = MPI_Recv(&tmp,0,MPI_INT,0,0,comm,&status);
    error = error && MPI_Rsend(send_buf,send_size,MPI_DOUBLE,0,0,comm);
    if (error != MPI_SUCCESS) throw ATC_Error("error in int_send "+to_string(error));
  }

  void int_scatter(MPI_Comm comm, int *send_buf, int *rec_buf, int count)
  {
    int error;
    int numprocs = size(comm);
    int sizes[numprocs];
    int displacements[numprocs];
    for (int i = 0; i < numprocs; ++i) {
      sizes[i] = 1;
      displacements[i] = i;
    }
    error = MPI_Scatterv(send_buf, sizes, displacements, MPI_INT, rec_buf, count, MPI_INT, 0, comm);
    if (error != MPI_SUCCESS) throw ATC_Error("error in int_scatter "+to_string(error));
  }

  void allgatherv(MPI_Comm comm, double *send_buf, int send_count,
                  double *rec_buf, int *rec_counts, int *displacements)
  {
    int error = MPI_Allgatherv(send_buf, send_count, MPI_DOUBLE,
                               rec_buf, rec_counts, displacements, MPI_DOUBLE,
                               comm);
    if (error != MPI_SUCCESS) throw ATC_Error("error in allgatherv "+to_string(error));
  }

  void gather(MPI_Comm comm, double send, double* recv)
  {
    int send_count = 1;
    int recv_count = 1;
    int root = 0;
    int error = MPI_Gather(&send, send_count, MPI_DOUBLE,
                            recv, recv_count,  MPI_DOUBLE, root, comm);
    if (error != MPI_SUCCESS) throw ATC_Error("error in allgatherv "+to_string(error));
  }

  void int_allgather(MPI_Comm comm, int send, int* recv)
  {
    int send_count = 1;
    int recv_count = 1;
    int error = MPI_Allgather(&send, send_count, MPI_INT,
                              recv, recv_count,  MPI_INT, comm);
    if (error != MPI_SUCCESS) throw ATC_Error("error in allgatherv "+to_string(error));
  }

  void logical_or(MPI_Comm comm, void *send_buf, int *rec_buf, int count)
  {
    int error = MPI_Allreduce(send_buf, rec_buf, count, MPI_INT, MPI_LOR,
                              comm);
    if (error != MPI_SUCCESS) throw ATC_Error("error in logical_or "+to_string(error));
  }

  void barrier(MPI_Comm comm)
  {
    int error = MPI_Barrier(comm);
    if (error != MPI_SUCCESS) throw ATC_Error("error in barrier "+to_string(error));
  }

  void stop(MPI_Comm comm, string msg)
  {
    int error = MPI_Barrier(comm);
    if (error != MPI_SUCCESS) throw ATC_Error("error in barrier "+to_string(error));
    throw ATC_Error("...stopping "+msg);
  }

#ifdef ISOLATE_FE
void sparse_allsum(MPI_Comm comm,SparseMatrix<double> &toShare) const
{
  toShare.compress();
  
  // initialize MPI information
  int nProcs = size(comm);
  int myRank = rank(comm);;

  int error;

  // get numbers of rows, columns, rowsCRS, and 
  // sizes (number of nonzero elements in matrix)
  SparseMatInfo *recInfo = new SparseMatInfo[nProcs];
  SparseMatInfo myInfo;
  myInfo.rows    = toShare.nRows();
  myInfo.cols    = toShare.nCols();
  myInfo.rowsCRS = toShare.nRowsCRS();
  myInfo.size    = toShare.size();

  error = MPI_Allgather(&myInfo, 4, MPI_INT,
                        recInfo, 4, MPI_INT, lammps_->world);
  if (error != MPI_SUCCESS) throw ATC_Error("error in sparse_allsum_numrows "+to_string(error));

  // adjust row sendcounts because recRowsCRS is off by one
  int rowCounts[nProcs];
  int sizeCounts[nProcs];
  // set up total size of receive buffers for Allgatherv calls
  int totalRowsCRS = 0;
  int totalSize = 0;
  // set up array of displacements for Allgatherv calls
  int rowOffsets[nProcs];
  rowOffsets[0] = 0;
  int sizeOffsets[nProcs];
  sizeOffsets[0] = 0;
  for (int i = 0; i < nProcs; i++) {
    // find the total number of entries to share in the mpi calls below
    rowCounts[i] = recInfo[i].rowsCRS + 1;
    sizeCounts[i] = recInfo[i].size;
    totalRowsCRS += rowCounts[i];
    totalSize += recInfo[i].size;
    // these already have their 0th slot filled in
    if (i == 0) continue;
    rowOffsets[i] = rowOffsets[i-1] + rowCounts[i-1];
    sizeOffsets[i] = sizeOffsets[i-1] + sizeCounts[i-1];
  }

  // get actual rows
  INDEX *rec_ia = new INDEX[totalRowsCRS];
  if (toShare.size() == 0) {
    double dummy[0];
    error = MPI_Allgatherv(dummy, 0, MPI_INT,
                           rec_ia, rowCounts, rowOffsets, MPI_INT, lammps_->world);
  }
  else
    error = MPI_Allgatherv(toShare.rows(), rowCounts[myRank], MPI_INT,
                           rec_ia, rowCounts, rowOffsets, MPI_INT, lammps_->world);
  if (error != MPI_SUCCESS)
    throw ATC_Error("error in sparse_allsum_rowarray "+to_string(error));

  // get actual cols
  INDEX *rec_ja = new INDEX[totalSize];
  error = MPI_Allgatherv(toShare.cols(), sizeCounts[myRank], MPI_INT,
                         rec_ja, sizeCounts, sizeOffsets, MPI_INT, lammps_->world);
  if (error != MPI_SUCCESS)
    throw ATC_Error("error in sparse_allsum_colarray "+to_string(error));
     
  // get the array of values
  double *rec_vals = new double[totalSize];
  error = MPI_Allgatherv(toShare.ptr(), sizeCounts[myRank], MPI_DOUBLE,
                         rec_vals, sizeCounts, sizeOffsets, MPI_DOUBLE, lammps_->world);
  if (error != MPI_SUCCESS)
    throw ATC_Error("error in sparse_allsum_valarray "+to_string(error));

  INDEX *rec_ia_proc; 
  INDEX *rec_ja_proc; 
  double *rec_vals_proc;
  for (int i = 0; i < nProcs; i++) {
    if (myRank != i) {
      // deallocated when tempMat is deleted since it wraps them
      rec_ia_proc = new INDEX[rowCounts[i]];
      rec_ja_proc = new INDEX[sizeCounts[i]];
      rec_vals_proc = new double[sizeCounts[i]];
       
      // copy the data passed with MPI into the new spots
      copy(rec_ia + rowOffsets[i], 
           rec_ia + rowOffsets[i] + rowCounts[i],
           rec_ia_proc);
      copy(rec_ja + sizeOffsets[i], 
           rec_ja + sizeOffsets[i] + sizeCounts[i],
           rec_ja_proc);
      copy(rec_vals + sizeOffsets[i], 
           rec_vals + sizeOffsets[i] + sizeCounts[i],
           rec_vals_proc);

      // Does anyone know why we have to declare tempMat here (as well as set it equal to
      // something) to avoid segfaults? there are still segfaults, but they happen at a much 
      // later stage of the game now (and for less benchmarks overall).
      SparseMatrix<double> tempMat =
        SparseMatrix<double>(rec_ia_proc, rec_ja_proc, rec_vals_proc, 
                             recInfo[i].size, recInfo[i].rows, 
                             recInfo[i].cols, recInfo[i].rowsCRS);
      toShare += tempMat;
    }
  }

  delete[] recInfo;
  delete[] rec_ia;
  delete[] rec_ja;
  delete[] rec_vals;
}
#endif

 void print_msg(MPI_Comm comm, string msg) 
  {
    if (serial(comm)) { cout << " ATC: " << msg << "\n"; }
    else { cout << " ATC: P" << rank(comm) << ", " << msg << "\n"; }
  }

  void print_msg_once(MPI_Comm comm, string msg, bool prefix, bool endline) 
  {
    if (rank_zero(comm)) {
      if (prefix) cout << " ATC: ";
      cout << msg;
      if (endline) cout << "\n";
    }
  }
}
