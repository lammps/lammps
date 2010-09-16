#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "many2one.h"
#include "memory.h"

/* ---------------------------------------------------------------------- */

Many2One::Many2One(MPI_Comm caller_comm)
{
  comm = caller_comm;
  MPI_Comm_rank(comm,&me);
  MPI_Comm_size(comm,&nprocs);

  memory = new Memory(comm);

  if (me == 0) {
    counts = new int[nprocs];
    multicounts = new int[nprocs];
    displs = new int[nprocs];
    multidispls = new int[nprocs];
  } else counts = multicounts = displs = multidispls = NULL;

  idall = NULL;
}

/* ---------------------------------------------------------------------- */

Many2One::~Many2One()
{
  delete memory;

  delete [] counts;
  delete [] multicounts;
  delete [] displs;
  delete [] multidispls;

  memory->sfree(idall);
}

/* ---------------------------------------------------------------------- */

void Many2One::setup(int nsrc_in, int *id, int ndest)
{
  nsrc = nsrc_in;
  MPI_Allreduce(&nsrc,&nall,1,MPI_INT,MPI_SUM,comm);
  MPI_Gather(&nsrc,1,MPI_INT,counts,1,MPI_INT,0,comm);

  if (me == 0) {
    displs[0] = 0;
    for (int i = 1; i < nprocs; i++)
      displs[i] = displs[i-1] + counts[i-1];
  }

  // gather IDs into idall

  idall = NULL;
  if (me == 0)
    idall = (int *) memory->smalloc(nall*sizeof(int),"many2one:idall");
  MPI_Gatherv(id,nsrc,MPI_INT,idall,counts,displs,MPI_INT,0,comm);
}

/* ---------------------------------------------------------------------- */

void Many2One::gather(double *src, int n, double *dest)
{
  int i,j,ii,jj,m;

  if (me == 0)
    for (int i = 0; i < nprocs; i++) {
      multicounts[i] = n*counts[i];
      multidispls[i] = n*displs[i];
    }

  // allgather src into desttmp

  double *desttmp = NULL;
  if (me == 0)
    desttmp = (double *) memory->smalloc(n*nall*sizeof(double),
					 "many2one:idsttmp");
  MPI_Gatherv(src,n*nsrc,MPI_DOUBLE,desttmp,multicounts,multidispls,
	      MPI_DOUBLE,0,comm);

  // use idall to move datums from desttmp to dest

  if (me == 0) {
    if (n == 1)
      for (i = 0; i < nall; i++) {
	j = idall[i] - 1;
	dest[j] = desttmp[i];
      }
    else
      for (i = 0; i < nall; i++) {
	j = idall[i] - 1;
	ii = n*i;
	jj = n*j;
	for (m = 0; m < n; m++)
	  dest[jj++] = desttmp[ii++];
      }
  }

  memory->sfree(desttmp);
}
