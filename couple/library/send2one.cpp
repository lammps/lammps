#include "mpi.h"
#include "stdlib.h"
#include "stdio.h"
#include "send2one.h"
#include "memory.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

Send2One::Send2One(MPI_Comm caller_comm)
{
  comm = caller_comm;
  MPI_Comm_rank(comm,&me);
  MPI_Comm_size(comm,&nprocs);

  memory = new Memory(comm);
  error = new Error(comm);

  buf = NULL;
  maxbuf = 0;
}

/* ---------------------------------------------------------------------- */

Send2One::~Send2One()
{
  delete memory;
  delete error;
  memory->sfree(buf);
}

/* ---------------------------------------------------------------------- */

void Send2One::execute()
{
  int nme,nmax,nsize,ping;
  MPI_Status status;
  MPI_Request request;

  // pre-processing before ping loop

  pre();

  // nme = size of data I contribute, in bytes
  // nmax = max size of data on any proc, in bytes
  // reallocate buf if necessary

  nme = size();
  MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,comm);

  if (nmax > maxbuf) {
    maxbuf = nmax;
    memory->sfree(buf);
    buf = (char *) memory->smalloc(maxbuf,"foo:buf");
  }

  // pack my data into buf

  pack(buf);

  // proc 0 pings each proc, receives its data
  // all other procs wait for ping, send their data to proc 0
  // invoke process() to work with data

  if (me == 0) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
	MPI_Irecv(buf,maxbuf,MPI_CHAR,iproc,0,comm,&request);
	MPI_Send(&ping,0,MPI_INT,iproc,0,comm);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_CHAR,&nsize);
      } else nsize = nme;
      
      process(nsize,buf);
    }
    
  } else {
    MPI_Recv(&ping,0,MPI_INT,0,0,comm,&status);
    MPI_Rsend(buf,nme,MPI_CHAR,0,0,comm);
  }

  post();
}
