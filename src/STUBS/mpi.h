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

#ifndef MPI_STUBS
#define MPI_STUBS

/* Dummy defs for MPI stubs */

#define MPI_COMM_WORLD 0

#define MPI_INT 1
#define MPI_FLOAT 2
#define MPI_DOUBLE 3
#define MPI_CHAR 4
#define MPI_BYTE 5
#define MPI_DOUBLE_INT 6

#define MPI_SUM 1
#define MPI_MAX 2
#define MPI_MIN 3
#define MPI_MAXLOC 4
#define MPI_MINLOC 5
#define MPI_LOR 6

#define MPI_ANY_SOURCE -1

#define MPI_Comm int
#define MPI_Request int
#define MPI_Datatype int
#define MPI_Op int

/* MPI data structs */

struct MPI_Status {
  int MPI_SOURCE;
};

/* Function prototypes for MPI stubs */

void MPI_Init(int *argc, char ***argv);
void MPI_Comm_rank(MPI_Comm comm, int *me);
void MPI_Comm_size(MPI_Comm comm, int *nprocs);
void MPI_Abort(MPI_Comm comm, int errorcode);
void MPI_Finalize();
double MPI_Wtime();

void MPI_Send(void *buf, int count, MPI_Datatype datatype,
	      int dest, int tag, MPI_Comm comm);
void MPI_Rsend(void *buf, int count, MPI_Datatype datatype,
	       int dest, int tag, MPI_Comm comm);
void MPI_Recv(void *buf, int count, MPI_Datatype datatype,
	      int source, int tag, MPI_Comm comm, MPI_Status *status);
void MPI_Irecv(void *buf, int count, MPI_Datatype datatype,
	       int source, int tag, MPI_Comm comm, MPI_Request *request);
void MPI_Wait(MPI_Request *request, MPI_Status *status);
void MPI_Waitall(int n, MPI_Request *request, MPI_Status *status);
void MPI_Waitany(int count, MPI_Request *request, int *index, 
		 MPI_Status *status);
void MPI_Sendrecv(void *sbuf, int scount, MPI_Datatype sdatatype,
		  int dest, int stag, void *rbuf, int rcount,
		  MPI_Datatype rdatatype, int source, int rtag,
		  MPI_Comm comm, MPI_Status *status);
void MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count);

void MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *comm_out);
void MPI_Comm_dup(MPI_Comm comm, MPI_Comm *comm_out);
void MPI_Comm_free(MPI_Comm *comm);

void MPI_Cart_create(MPI_Comm comm_old, int ndims, int *dims, int *periods,
		     int reorder, MPI_Comm *comm_cart);
void MPI_Cart_get(MPI_Comm comm, int maxdims, int *dims, int *periods,
		  int *coords);
void MPI_Cart_shift(MPI_Comm comm, int direction, int displ,
		    int *source, int *dest);
void MPI_Cart_rank(MPI_Comm comm, int *coords, int *rank);

void MPI_Barrier(MPI_Comm comm);
void MPI_Bcast(void *buf, int count, MPI_Datatype datatype,
	       int root, MPI_Comm comm);
void MPI_Allreduce(void *sendbuf, void *recvbuf, int count,
		   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
void MPI_Reduce(void *sendbuf, void *recvbuf, int count,
		   MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
void MPI_Scan(void *sendbuf, void *recvbuf, int count,
	      MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
void MPI_Allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		   void *recvbuf, int recvcount, MPI_Datatype recvtype,
		   MPI_Comm comm);
void MPI_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		    void *recvbuf, int *recvcounts, int *displs,
		    MPI_Datatype recvtype, MPI_Comm comm);
void MPI_Reduce_scatter(void *sendbuf, void *recvbuf, int *recvcounts,
			MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
void MPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		void *recvbuf, int recvcount, MPI_Datatype recvtype,
		int root, MPI_Comm comm);
void MPI_Gatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		    void *recvbuf, int *recvcounts, int *displs,
		    MPI_Datatype recvtype, int root, MPI_Comm comm);

#endif
