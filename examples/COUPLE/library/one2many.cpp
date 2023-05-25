#include <mpi.h>
#include <cstdlib>
#include "one2many.h"
#include "memorylib.h"

#include <map>

/* ---------------------------------------------------------------------- */

One2Many::One2Many(MPI_Comm caller_comm)
{
  comm = caller_comm;
  MPI_Comm_rank(comm,&me);
  MPI_Comm_size(comm,&nprocs);
  
  memory = new MemoryLib(comm);
  hash = new std::map<int,int>();
}

/* ---------------------------------------------------------------------- */

One2Many::~One2Many()
{
  delete memory;
  delete hash;
}

/* ---------------------------------------------------------------------- */

void One2Many::setup(int nsrc_in, int ndest, int *id)
{
  nsrc = nsrc_in;

  // store my local IDs in hash

  hash->clear();
  for (int i = 0; i < ndest; i++)
    hash->insert(std::pair<int,int> (id[i],i));
}

/* ---------------------------------------------------------------------- */

void One2Many::scatter(double *src, int n, double *dest)
{
  int i,j,k,m;
 
  // allocate src on procs that don't have it

  int flag = 0;
  if (src == NULL) {
    src = (double *) memory->smalloc(n*nsrc*sizeof(double),"one2many:src");
    flag = 1;
  }

  // broadcast src from 0 to other procs

  MPI_Bcast(src,n*nsrc,MPI_DOUBLE,0,comm);

  // each proc loops over entire src
  // if I own the global ID, copy src values into dest

  std::map<int,int>::iterator loc;
  for (m = 1; m <= nsrc; m++) {
    loc = hash->find(m);
    if (loc == hash->end()) continue;
    i = n*loc->second;
    j = 3*(m-1);
    if (n == 1) dest[i] = src[j];
    else
      for (k = 0; k < n; k++)
	dest[i++] = src[j++];
  }

  // free locally allocated src

  if (flag) memory->sfree(src);
}
