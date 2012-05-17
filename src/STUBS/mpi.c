/* -----------------------------------------------------------------------
   LAMMPS 2003 (July 31) - Molecular Dynamics Simulator
   Sandia National Laboratories, www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------ */

/* Single-processor "stub" versions of MPI routines */

#include "stdlib.h"
#include "string.h"
#include "stdio.h"
#include "stdint.h"
#include <sys/time.h>
#include "mpi.h"

/* lo-level data structure */

struct _mpi_double_int {
  double value;
  int proc;
};
typedef struct _mpi_double_int double_int;

/* ---------------------------------------------------------------------- */
/* MPI Functions */
/* ---------------------------------------------------------------------- */

int MPI_Init(int *argc, char ***argv) {return 0;}

/* ---------------------------------------------------------------------- */

int MPI_Initialized(int *flag)
{
  *flag = 1;
  return 0;
}

/* ---------------------------------------------------------------------- */

/* return "localhost" as name of the processor */

void MPI_Get_processor_name(char *name, int *resultlen)
{
  const char host[] = "localhost";
  int len;

  if (!name || !resultlen) return;

  len = strlen(host);
  memcpy(name,host,len+1);
  *resultlen = len;
  return;
}

/* ---------------------------------------------------------------------- */

int MPI_Comm_rank(MPI_Comm comm, int *me)
{
  *me = 0;
  return 0;
}

/* ---------------------------------------------------------------------- */

int MPI_Comm_size(MPI_Comm comm, int *nprocs)
{
  *nprocs = 1;
  return 0;
}

/* ---------------------------------------------------------------------- */

int MPI_Abort(MPI_Comm comm, int errorcode)
{
  exit(1);
  return 0;
}

/* ---------------------------------------------------------------------- */

int MPI_Finalize() {return 0;}

/* ---------------------------------------------------------------------- */

double MPI_Wtime()
{
  double time;
  struct timeval tv;

  gettimeofday(&tv,NULL);
  time = 1.0 * tv.tv_sec + 1.0e-6 * tv.tv_usec;
  return time;
}

/* ---------------------------------------------------------------------- */

int MPI_Type_size(MPI_Datatype datatype, int *size)
{
  if (datatype == MPI_INT) *size = sizeof(int);
  else if (datatype == MPI_FLOAT) *size = sizeof(float);
  else if (datatype == MPI_DOUBLE) *size = sizeof(double);
  else if (datatype == MPI_CHAR) *size = sizeof(char);
  else if (datatype == MPI_BYTE) *size = sizeof(char);
  else if (datatype == MPI_LONG_LONG) *size = sizeof(uint64_t);
  else if (datatype == MPI_DOUBLE_INT) *size = sizeof(double_int);

  return 0;
}

/* ---------------------------------------------------------------------- */

int MPI_Send(void *buf, int count, MPI_Datatype datatype,
             int dest, int tag, MPI_Comm comm)
{
  printf("MPI Stub WARNING: Should not send message to self\n");
  return 0;
}

/* ---------------------------------------------------------------------- */

int MPI_Isend(void *buf, int count, MPI_Datatype datatype,
              int source, int tag, MPI_Comm comm, MPI_Request *request)
{
  printf("MPI Stub WARNING: Should not send message to self\n");
  return 0;
}

/* ---------------------------------------------------------------------- */

int MPI_Rsend(void *buf, int count, MPI_Datatype datatype,
              int dest, int tag, MPI_Comm comm)
{
  printf("MPI Stub WARNING: Should not rsend message to self\n");
  return 0;
}

/* ---------------------------------------------------------------------- */

int MPI_Recv(void *buf, int count, MPI_Datatype datatype,
             int source, int tag, MPI_Comm comm, MPI_Status *status)
{
  printf("MPI Stub WARNING: Should not recv message from self\n");
  return 0;
}

/* ---------------------------------------------------------------------- */

int MPI_Irecv(void *buf, int count, MPI_Datatype datatype,
              int source, int tag, MPI_Comm comm, MPI_Request *request)
{
  printf("MPI Stub WARNING: Should not recv message from self\n");
  return 0;
}

/* ---------------------------------------------------------------------- */

int MPI_Wait(MPI_Request *request, MPI_Status *status)
{
  printf("MPI Stub WARNING: Should not wait on message from self\n");
  return 0;
}

/* ---------------------------------------------------------------------- */

int MPI_Waitall(int n, MPI_Request *request, MPI_Status *status)
{
  printf("MPI Stub WARNING: Should not wait on message from self\n");
  return 0;
}

/* ---------------------------------------------------------------------- */

int MPI_Waitany(int count, MPI_Request *request, int *index, 
                MPI_Status *status)
{
  printf("MPI Stub WARNING: Should not wait on message from self\n");
  return 0;
}

/* ---------------------------------------------------------------------- */

int MPI_Sendrecv(void *sbuf, int scount, MPI_Datatype sdatatype,
                 int dest, int stag, void *rbuf, int rcount,
                 MPI_Datatype rdatatype, int source, int rtag,
                 MPI_Comm comm, MPI_Status *status)
{
  printf("MPI Stub WARNING: Should not send message to self\n");
  return 0;
}

/* ---------------------------------------------------------------------- */

int MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count)
{
  printf("MPI Stub WARNING: Should not get count of message to self\n");
  return 0;
}

/* ---------------------------------------------------------------------- */

int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *comm_out)
{
  *comm_out = comm;
  return 0;
}

/* ---------------------------------------------------------------------- */

int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *comm_out)
{
  *comm_out = comm;
  return 0;
}

/* ---------------------------------------------------------------------- */

int MPI_Comm_free(MPI_Comm *comm) {return 0;}

/* ---------------------------------------------------------------------- */

int MPI_Cart_create(MPI_Comm comm_old, int ndims, int *dims, int *periods,
                    int reorder, MPI_Comm *comm_cart)
{
  *comm_cart = comm_old;
  return 0;
}

/* ---------------------------------------------------------------------- */

int MPI_Cart_get(MPI_Comm comm, int maxdims, int *dims, int *periods,
                 int *coords)
{
  dims[0] = dims[1] = dims[2] = 1;
  periods[0] = periods[1] = periods[2] = 1;
  coords[0] = coords[1] = coords[2] = 0;
  return 0;
}

/* ---------------------------------------------------------------------- */

int MPI_Cart_shift(MPI_Comm comm, int direction, int displ,
                   int *source, int *dest)
{
  *source = *dest = 0;
  return 0;
}

/* ---------------------------------------------------------------------- */

int MPI_Cart_rank(MPI_Comm comm, int *coords, int *rank)
{
  *rank = 0;
  return 0;
}

/* ---------------------------------------------------------------------- */

int MPI_Barrier(MPI_Comm comm) {return 0;}

/* ---------------------------------------------------------------------- */

int MPI_Bcast(void *buf, int count, MPI_Datatype datatype,
              int root, MPI_Comm comm) {return 0;}

/* ---------------------------------------------------------------------- */

/* copy values from data1 to data2 */

int MPI_Allreduce(void *sendbuf, void *recvbuf, int count,
                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  int n;
  if (datatype == MPI_INT) n = count*sizeof(int);
  else if (datatype == MPI_FLOAT) n = count*sizeof(float);
  else if (datatype == MPI_DOUBLE) n = count*sizeof(double);
  else if (datatype == MPI_CHAR) n = count*sizeof(char);
  else if (datatype == MPI_BYTE) n = count*sizeof(char);
  else if (datatype == MPI_LONG_LONG) n = count*sizeof(uint64_t);
  else if (datatype == MPI_DOUBLE_INT) n = count*sizeof(double_int);

  if (sendbuf == MPI_IN_PLACE || recvbuf == MPI_IN_PLACE) return;
  memcpy(recvbuf,sendbuf,n);
  return 0;
}

/* ---------------------------------------------------------------------- */

/* copy values from data1 to data2 */

int MPI_Reduce(void *sendbuf, void *recvbuf, int count,
		   MPI_Datatype datatype, MPI_Op op,
		   int root, MPI_Comm comm)
{
  int n;
  if (datatype == MPI_INT) n = count*sizeof(int);
  else if (datatype == MPI_FLOAT) n = count*sizeof(float);
  else if (datatype == MPI_DOUBLE) n = count*sizeof(double);
  else if (datatype == MPI_CHAR) n = count*sizeof(char);
  else if (datatype == MPI_BYTE) n = count*sizeof(char);
  else if (datatype == MPI_LONG_LONG) n = count*sizeof(uint64_t);
  else if (datatype == MPI_DOUBLE_INT) n = count*sizeof(double_int);

  if (sendbuf == MPI_IN_PLACE || recvbuf == MPI_IN_PLACE) return;
  memcpy(recvbuf,sendbuf,n);
  return 0;
}


/* ---------------------------------------------------------------------- */

int MPI_Scan(void *sendbuf, void *recvbuf, int count,
             MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  int n;
  if (datatype == MPI_INT) n = count*sizeof(int);
  else if (datatype == MPI_FLOAT) n = count*sizeof(float);
  else if (datatype == MPI_DOUBLE) n = count*sizeof(double);
  else if (datatype == MPI_CHAR) n = count*sizeof(char);
  else if (datatype == MPI_BYTE) n = count*sizeof(char);
  else if (datatype == MPI_LONG_LONG) n = count*sizeof(uint64_t);
  else if (datatype == MPI_DOUBLE_INT) n = count*sizeof(double_int);

  if (sendbuf == MPI_IN_PLACE || recvbuf == MPI_IN_PLACE) return;
  memcpy(recvbuf,sendbuf,n);
  return 0;
}

/* ---------------------------------------------------------------------- */

/* copy values from data1 to data2 */

int MPI_Allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                  void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  MPI_Comm comm)
{
  int n;
  if (sendtype == MPI_INT) n = sendcount*sizeof(int);
  else if (sendtype == MPI_FLOAT) n = sendcount*sizeof(float);
  else if (sendtype == MPI_DOUBLE) n = sendcount*sizeof(double);
  else if (sendtype == MPI_CHAR) n = sendcount*sizeof(char);
  else if (sendtype == MPI_BYTE) n = sendcount*sizeof(char);
  else if (sendtype == MPI_LONG_LONG) n = sendcount*sizeof(uint64_t);
  else if (sendtype == MPI_DOUBLE_INT) n = sendcount*sizeof(double_int);

  if (sendbuf == MPI_IN_PLACE || recvbuf == MPI_IN_PLACE) return;
  memcpy(recvbuf,sendbuf,n);
  return 0;
}

/* ---------------------------------------------------------------------- */

/* copy values from data1 to data2 */

int MPI_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                   void *recvbuf, int *recvcounts, int *displs,
                   MPI_Datatype recvtype, MPI_Comm comm)
{
  int n;
  if (sendtype == MPI_INT) n = sendcount*sizeof(int);
  else if (sendtype == MPI_FLOAT) n = sendcount*sizeof(float);
  else if (sendtype == MPI_DOUBLE) n = sendcount*sizeof(double);
  else if (sendtype == MPI_CHAR) n = sendcount*sizeof(char);
  else if (sendtype == MPI_BYTE) n = sendcount*sizeof(char);
  else if (sendtype == MPI_LONG_LONG) n = sendcount*sizeof(uint64_t);
  else if (sendtype == MPI_DOUBLE_INT) n = sendcount*sizeof(double_int);

  if (sendbuf == MPI_IN_PLACE || recvbuf == MPI_IN_PLACE) return;
  memcpy(recvbuf,sendbuf,n);
  return 0;
}

/* ---------------------------------------------------------------------- */

/* copy values from data1 to data2 */

int MPI_Reduce_scatter(void *sendbuf, void *recvbuf, int *recvcounts,
                       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  int n;
  if (datatype == MPI_INT) n = *recvcounts*sizeof(int);
  else if (datatype == MPI_FLOAT) n = *recvcounts*sizeof(float);
  else if (datatype == MPI_DOUBLE) n = *recvcounts*sizeof(double);
  else if (datatype == MPI_CHAR) n = *recvcounts*sizeof(char);
  else if (datatype == MPI_BYTE) n = *recvcounts*sizeof(char);
  else if (datatype == MPI_LONG_LONG) n = *recvcounts*sizeof(uint64_t);
  else if (datatype == MPI_DOUBLE_INT) n = *recvcounts*sizeof(double_int);

  if (sendbuf == MPI_IN_PLACE || recvbuf == MPI_IN_PLACE) return;
  memcpy(recvbuf,sendbuf,n);
  return 0;
}

/* ---------------------------------------------------------------------- */

/* copy values from data1 to data2 */

int MPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
               void *recvbuf, int recvcount, MPI_Datatype recvtype,
               int root, MPI_Comm comm)
{
  int n;
  if (sendtype == MPI_INT) n = sendcount*sizeof(int);
  else if (sendtype == MPI_FLOAT) n = sendcount*sizeof(float);
  else if (sendtype == MPI_DOUBLE) n = sendcount*sizeof(double);
  else if (sendtype == MPI_CHAR) n = sendcount*sizeof(char);
  else if (sendtype == MPI_BYTE) n = sendcount*sizeof(char);
  else if (sendtype == MPI_LONG_LONG) n = sendcount*sizeof(uint64_t);
  else if (sendtype == MPI_DOUBLE_INT) n = sendcount*sizeof(double_int);

  if (sendbuf == MPI_IN_PLACE || recvbuf == MPI_IN_PLACE) return;
  memcpy(recvbuf,sendbuf,n);
  return 0;
}

/* ---------------------------------------------------------------------- */

/* copy values from data1 to data2 */

int MPI_Gatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		void *recvbuf, int *recvcounts, int *displs,
		MPI_Datatype recvtype, int root, MPI_Comm comm)
{
  int n;
  if (sendtype == MPI_INT) n = sendcount*sizeof(int);
  else if (sendtype == MPI_FLOAT) n = sendcount*sizeof(float);
  else if (sendtype == MPI_DOUBLE) n = sendcount*sizeof(double);
  else if (sendtype == MPI_CHAR) n = sendcount*sizeof(char);
  else if (sendtype == MPI_BYTE) n = sendcount*sizeof(char);
  else if (sendtype == MPI_LONG_LONG) n = sendcount*sizeof(uint64_t);
  else if (sendtype == MPI_DOUBLE_INT) n = sendcount*sizeof(double_int);

  if (sendbuf == MPI_IN_PLACE || recvbuf == MPI_IN_PLACE) return;
  memcpy(recvbuf,sendbuf,n);
  return 0;
}

/* ---------------------------------------------------------------------- */

/* copy values from data1 to data2 */

int MPI_Scatterv(void *sendbuf, int *sendcounts, int *displs,
		 MPI_Datatype sendtype, void *recvbuf, int recvcount,
		 MPI_Datatype recvtype, int root, MPI_Comm comm)
{
  int n;
  if (sendtype == MPI_INT) n = recvcount*sizeof(int);
  else if (sendtype == MPI_FLOAT) n = recvcount*sizeof(float);
  else if (sendtype == MPI_DOUBLE) n = recvcount*sizeof(double);
  else if (sendtype == MPI_CHAR) n = recvcount*sizeof(char);
  else if (sendtype == MPI_BYTE) n = recvcount*sizeof(char);
  else if (sendtype == MPI_LONG_LONG) n = recvcount*sizeof(uint64_t);
  else if (sendtype == MPI_DOUBLE_INT) n = recvcount*sizeof(double_int);

  if (sendbuf == MPI_IN_PLACE || recvbuf == MPI_IN_PLACE) return;
  memcpy(recvbuf,sendbuf,n);
  return 0;
}
