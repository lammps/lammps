/* -*- c++ -*- -----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------ */

#ifndef MPI_STUBS
#define MPI_STUBS

#include <stdlib.h>

/* We compile STUBS with C++ so the symbols embedded
 * the serial shared library will not collide with any
 * corresponding symbols from a real MPI library (which
 * uses C bindings). As a consequence the header *must*
 * enforce compiling with C++ only. */

#ifndef __cplusplus
#error "MPI STUBS must be compiled with a C++ compiler"
#endif

/* Dummy defs for MPI stubs */

#define MPI_COMM_WORLD 0

#define MPI_SUCCESS 0
#define MPI_ERR_ARG -1

#define MPI_INT 1
#define MPI_FLOAT 2
#define MPI_DOUBLE 3
#define MPI_CHAR 4
#define MPI_BYTE 5
#define MPI_LONG 6
#define MPI_LONG_LONG 7
#define MPI_DOUBLE_INT 8

#define MPI_SUM 1
#define MPI_MAX 2
#define MPI_MIN 3
#define MPI_MAXLOC 4
#define MPI_MINLOC 5
#define MPI_LOR 6

#define MPI_UNDEFINED -1
#define MPI_COMM_NULL -1
#define MPI_GROUP_EMPTY -1
#define MPI_GROUP_NULL -1

#define MPI_ANY_SOURCE -1
#define MPI_STATUS_IGNORE NULL

#define MPI_TAG_UB 0

#define MPI_Comm int
#define MPI_Request int
#define MPI_Datatype int
#define MPI_Op int
#define MPI_Fint int
#define MPI_Group int
#define MPI_Offset long

#define MPI_IN_PLACE NULL

#define MPI_MAX_PROCESSOR_NAME 128
#define MPI_MAX_LIBRARY_VERSION_STRING 128

typedef void MPI_User_function(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype);

/* MPI data structs */

struct _MPI_Status {
  int MPI_SOURCE;
};
typedef struct _MPI_Status MPI_Status;

/* Function prototypes for MPI stubs */

int MPI_Init(int *argc, char ***argv);
int MPI_Initialized(int *flag);
int MPI_Finalized(int *flag);
int MPI_Get_library_version(char *version, int *resultlen);
int MPI_Get_processor_name(char *name, int *resultlen);
int MPI_Get_version(int *major, int *minor);

int MPI_Comm_rank(MPI_Comm comm, int *me);
int MPI_Comm_size(MPI_Comm comm, int *nprocs);
int MPI_Abort(MPI_Comm comm, int errorcode);
int MPI_Finalize();
double MPI_Wtime();

int MPI_Type_size(int, int *);
int MPI_Request_free(MPI_Request *request);

int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
              MPI_Request *request);
int MPI_Rsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
             MPI_Status *status);
int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
              MPI_Request *request);
int MPI_Wait(MPI_Request *request, MPI_Status *status);
int MPI_Waitall(int n, MPI_Request *request, MPI_Status *status);
int MPI_Waitany(int count, MPI_Request *request, int *index, MPI_Status *status);
int MPI_Sendrecv(const void *sbuf, int scount, MPI_Datatype sdatatype, int dest, int stag,
                 void *rbuf, int rcount, MPI_Datatype rdatatype, int source, int rtag,
                 MPI_Comm comm, MPI_Status *status);
int MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count);

int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *comm_out);
int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *comm_out);
int MPI_Comm_free(MPI_Comm *comm);
MPI_Fint MPI_Comm_c2f(MPI_Comm comm);
MPI_Comm MPI_Comm_f2c(MPI_Fint comm);
int MPI_Comm_group(MPI_Comm comm, MPI_Group *group);
int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm);
int MPI_Group_incl(MPI_Group group, int n, int *ranks, MPI_Group *newgroup);
int MPI_Group_free(MPI_Group *group);

int MPI_Cart_create(MPI_Comm comm_old, int ndims, int *dims, int *periods, int reorder,
                    MPI_Comm *comm_cart);
int MPI_Cart_get(MPI_Comm comm, int maxdims, int *dims, int *periods, int *coords);
int MPI_Cart_shift(MPI_Comm comm, int direction, int displ, int *source, int *dest);
int MPI_Cart_rank(MPI_Comm comm, int *coords, int *rank);

int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype);
int MPI_Type_commit(MPI_Datatype *datatype);
int MPI_Type_free(MPI_Datatype *datatype);

int MPI_Op_create(MPI_User_function *function, int commute, MPI_Op *op);
int MPI_Op_free(MPI_Op *op);

int MPI_Barrier(MPI_Comm comm);
int MPI_Bcast(void *buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm);
int MPI_Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                  MPI_Comm comm);
int MPI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root,
               MPI_Comm comm);
int MPI_Scan(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
             MPI_Comm comm);
int MPI_Allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
                  MPI_Datatype recvtype, MPI_Comm comm);
int MPI_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                   int *recvcounts, int *displs, MPI_Datatype recvtype, MPI_Comm comm);
int MPI_Reduce_scatter(void *sendbuf, void *recvbuf, int *recvcounts, MPI_Datatype datatype,
                       MPI_Op op, MPI_Comm comm);
int MPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
               MPI_Datatype recvtype, int root, MPI_Comm comm);
int MPI_Gatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int *recvcounts,
                int *displs, MPI_Datatype recvtype, int root, MPI_Comm comm);
int MPI_Scatter(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
                MPI_Datatype recvtype, int root, MPI_Comm comm);
int MPI_Scatterv(void *sendbuf, int *sendcounts, int *displs, MPI_Datatype sendtype, void *recvbuf,
                 int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
int MPI_Alltoall(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
                 MPI_Datatype recvtype, MPI_Comm comm);
int MPI_Alltoallv(void *sendbuf, int *sendcounts, int *sdispls, MPI_Datatype sendtype,
                  void *recvbuf, int *recvcounts, int *rdispls, MPI_Datatype recvtype,
                  MPI_Comm comm);
int MPI_Comm_get_attr(MPI_Comm comm, int comm_keyval, void **attribute_val,
                      int *flag);
/* ---------------------------------------------------------------------- */

#endif
